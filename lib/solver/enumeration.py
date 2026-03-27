"""
Обобщённые помощники для конечного перебора кандидатов-генераторов.

Обобщение lib/numeric/enumeration.py на произвольное число переменных n.
BasisTerm хранит мультииндекс alpha как кортеж длины n вместо фиксированных ax, ay.

Функции build_basis_terms, iter_candidates, build_generator и вспомогательные
переиспользуются из lib/numeric/enumeration.py для n=2 через backward-compat shim.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Iterator

from lib.numeric.types import Generator


@dataclass(frozen=True)
class BasisTerm:
    """
    Базисный терм — одно мономиальное слагаемое генератора для произвольного n.

    Атрибуты
    ---------
    alpha      : мультииндекс степеней переменных, len(alpha) == n
    target_var : индекс производной (0-based, в диапазоне 0..n-1)
    n          : число переменных
    """

    alpha: tuple[int, ...]
    target_var: int
    n: int

    def __post_init__(self) -> None:
        if len(self.alpha) != self.n:
            raise ValueError(
                f"alpha должен иметь длину n={self.n}, получено {len(self.alpha)}"
            )
        if not (0 <= self.target_var < self.n):
            raise ValueError(
                f"target_var={self.target_var} вне диапазона [0, {self.n})"
            )

    def factor_spec(self, variables: list[str] | None = None) -> str:
        """
        Форматирует терм в человекочитаемую строку вида 'x^2*y*dz'.

        Parameters
        ----------
        variables : список имён переменных длины n; по умолчанию ['x','y','z','u','v','w']
        """
        _default_vars = ["x", "y", "z", "u", "v", "w"]
        if variables is None:
            if self.n <= len(_default_vars):
                variables = _default_vars[: self.n]
            else:
                variables = [f"z{i}" for i in range(self.n)]

        factors: list[str] = []
        for i, exp in enumerate(self.alpha):
            if exp == 1:
                factors.append(variables[i])
            elif exp > 1:
                factors.append(f"{variables[i]}^{exp}")

        deriv = "d" + variables[self.target_var]
        factors.append(deriv)
        return "*".join(factors)


def build_basis_terms(
    n: int,
    max_degree: int,
    include_constants: bool = False,
) -> list[BasisTerm]:
    """
    Строит список базисных термов для пространства с n переменными.

    Parameters
    ----------
    n              : число переменных
    max_degree     : максимальная суммарная степень мономиального коэффициента
    include_constants : включить ли константные члены (alpha = (0,...,0))

    Returns
    -------
    Список BasisTerm, охватывающий все комбинации alpha и target_var.
    """
    terms: list[BasisTerm] = []
    # Перебираем все суммарные степени от 0 до max_degree
    for total in range(max_degree + 1):
        for alpha in _multiindices(n, total):
            if not include_constants and all(a == 0 for a in alpha):
                continue
            for target_var in range(n):
                terms.append(BasisTerm(alpha=alpha, target_var=target_var, n=n))
    return terms


def _multiindices(n: int, total: int) -> Iterator[tuple[int, ...]]:
    """Генератор всех мультиндексов длины n с суммой total."""
    if n == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _multiindices(n - 1, total - first):
            yield (first,) + rest


def build_generator(candidate: list[tuple[BasisTerm, float]]) -> Generator:
    """
    Строит Generator dict из списка (терм, коэффициент).

    Parameters
    ----------
    candidate : список пар (BasisTerm, float)

    Returns
    -------
    Generator: dict[int, dict[tuple, float]]
    """
    generator: Generator = {}
    for term, coeff in candidate:
        generator.setdefault(term.target_var, {})
        alpha = term.alpha
        generator[term.target_var][alpha] = (
            generator[term.target_var].get(alpha, 0.0) + coeff
        )
    # Убираем нулевые мономы и пустые направления
    result: Generator = {}
    for target, monomials in generator.items():
        clean = {a: c for a, c in monomials.items() if abs(c) > 1e-15}
        if clean:
            result[target] = clean
    return result


def estimate_candidate_count(
    basis_size: int,
    coeff_count: int,
    min_terms: int,
    max_terms: int,
) -> int:
    """Оценка числа кандидатов для заданного базиса и множества коэффициентов."""
    total = 0
    for term_count in range(min_terms, max_terms + 1):
        total += math.comb(basis_size, term_count) * (coeff_count**term_count)
    return total


def iter_candidates(
    basis_terms: list[BasisTerm],
    coeffs: list[float],
    min_terms: int,
    max_terms: int,
) -> Iterator[list[tuple[BasisTerm, float]]]:
    """
    Итератор по всем кандидатам-генераторам.

    Перебирает все подмножества basis_terms размером min_terms..max_terms
    с каждым из коэффициентов из coeffs.
    """
    for term_count in range(min_terms, max_terms + 1):
        for term_combo in itertools.combinations(basis_terms, term_count):
            for coeff_combo in itertools.product(coeffs, repeat=term_count):
                yield list(zip(term_combo, coeff_combo))


def format_coefficient(value: float) -> str:
    """Форматирует коэффициент: целые числа без дробной части."""
    if abs(value - round(value)) <= 1e-12:
        return str(int(round(value)))
    return f"{value:g}"


def candidate_to_spec(
    candidate: list[tuple[BasisTerm, float]],
    variables: list[str] | None = None,
) -> str:
    """Форматирует кандидата в строку вида '2*x^2*dx + y*dy'."""
    parts = []
    for term, coeff in candidate:
        factor = term.factor_spec(variables=variables)
        if abs(coeff - 1.0) <= 1e-12:
            parts.append(factor)
        else:
            parts.append(f"{format_coefficient(coeff)}*{factor}")
    return " + ".join(parts)


def parse_coefficients(text: str) -> list[float]:
    """
    Разбирает строку коэффициентов вида '1,-1,2'.

    Дедуплицирует по абсолютному значению (с допуском 1e-15).
    """
    coeffs: list[float] = []
    for raw in text.split(","):
        item = raw.strip()
        if not item:
            continue
        value = float(item)
        if abs(value) <= 1e-15:
            continue
        if not any(abs(value - seen) <= 1e-15 for seen in coeffs):
            coeffs.append(value)
    if not coeffs:
        raise ValueError("нужен хотя бы один ненулевой коэффициент в --coeffs")
    return coeffs
