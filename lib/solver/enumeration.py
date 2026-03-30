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
    Построить полный список мономиальных базисных термов для перебора генераторов.

    На вход принимает число переменных ``n``, максимальную степень коэффициента ``max_degree`` и флаг ``include_constants``.
    На выходе возвращает список ``BasisTerm``, покрывающий все допустимые пары ``(alpha, target_var)``.

    >>> len(build_basis_terms(2, 1, include_constants=False))
    4
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
    Собрать разреженный generator из конечного набора мономиальных термов.

    На вход принимает список пар ``(BasisTerm, coeff)``.
    На выходе возвращает generator в формате ``{axis: {alpha: coeff}}`` с удалением нулевых коэффициентов.

    >>> term = BasisTerm(alpha=(1, 0), target_var=0, n=2)
    >>> build_generator([(term, 2.0)])
    {0: {(1, 0): 2.0}}
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
    """
    Оценить мощность конечного пространства перебора генераторов.

    На вход принимает размер базиса ``basis_size``, число коэффициентов ``coeff_count`` и границы числа термов.
    На выходе возвращает суммарное число возможных кандидатов.

    >>> estimate_candidate_count(3, 2, 1, 2)
    18
    """
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
    Итеративно перечислить все кандидаты конечного перебора.

    На вход принимает набор ``basis_terms``, допустимые коэффициенты ``coeffs`` и диапазон числа термов.
    На выходе возвращает итератор по спискам пар ``(BasisTerm, coeff)``.

    >>> basis = build_basis_terms(1, 1, include_constants=False)
    >>> len(list(iter_candidates(basis, [1.0], 1, 1)))
    1
    """
    for term_count in range(min_terms, max_terms + 1):
        for term_combo in itertools.combinations(basis_terms, term_count):
            for coeff_combo in itertools.product(coeffs, repeat=term_count):
                yield list(zip(term_combo, coeff_combo))


def format_coefficient(value: float) -> str:
    """
    Отформатировать коэффициент для строковой спецификации кандидата.

    На вход принимает вещественное число ``value``.
    На выходе возвращает строку без дробной части для почти целых чисел и обычный ``g``-формат иначе.

    >>> format_coefficient(2.0), format_coefficient(0.5)
    ('2', '0.5')
    """
    if abs(value - round(value)) <= 1e-12:
        return str(int(round(value)))
    return f"{value:g}"


def candidate_to_spec(
    candidate: list[tuple[BasisTerm, float]],
    variables: list[str] | None = None,
) -> str:
    """
    Перевести кандидата из списка термов в строковую спецификацию.

    На вход принимает список пар ``(BasisTerm, coeff)`` и необязательные имена переменных.
    На выходе возвращает строку вида ``2*x^2*dx + y*dy``.

    >>> term = BasisTerm(alpha=(1, 0), target_var=1, n=2)
    >>> candidate_to_spec([(term, 1.0)])
    'x*dy'
    """
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
    Разобрать пользовательскую строку коэффициентов для перебора.

    На вход принимает строку ``text`` вида ``1,-1,2``.
    На выходе возвращает список ненулевых коэффициентов с удалением близких дубликатов.

    >>> parse_coefficients("1, -1, 1, 0, 2")
    [1.0, -1.0, 2.0]
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
