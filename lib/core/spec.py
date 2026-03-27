"""
Каноническая модель дифференцирования: DerivationSpec и примитивные кирпичики.

DerivationSpec — единое представление дифференцирования из алгебры Ли W_n,
не зависящее ни от Sage, ни от NumPy. Все формулы теорем сначала выражаются
через DerivationSpec, а затем адаптируются конкретным backend-ом.
"""

from __future__ import annotations

import random as _random
from dataclasses import dataclass, field
from typing import Union

# Скалярный тип коэффициента: int, float или рациональное число Python.
Scalar = Union[int, float]


@dataclass(frozen=True)
class DerivationSpec:
    """
    Каноническое представление дифференцирования из W_n.

    Семантика:
        terms[i][alpha] = коэффициент при z^alpha * ∂/∂z_i

    Пример (x^2*y*∂/∂x + 3*∂/∂y в W_2):
        DerivationSpec(n=2, terms={0: {(2, 1): 1.0}, 1: {(0, 0): 3.0}})

    Инварианты, гарантируемые при создании:
    - target_var (ключ terms) находится в диапазоне 0..n-1
    - len(alpha) == n для каждого мультиидекса alpha
    - нулевые коэффициенты удаляются при нормализации
    """

    n: int
    terms: dict[int, dict[tuple[int, ...], Scalar]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.n <= 0:
            raise ValueError(f"n должно быть положительным, получено n={self.n}")

        # Валидация и нормализация
        normalized: dict[int, dict[tuple[int, ...], Scalar]] = {}
        for target_var, mono_dict in self.terms.items():
            if not (0 <= target_var < self.n):
                raise ValueError(
                    f"target_var={target_var} вне диапазона 0..{self.n - 1}"
                )
            cleaned: dict[tuple[int, ...], Scalar] = {}
            for alpha, coeff in mono_dict.items():
                if len(alpha) != self.n:
                    raise ValueError(
                        f"len(alpha)={len(alpha)} != n={self.n} "
                        f"для alpha={alpha} в компоненте {target_var}"
                    )
                if isinstance(coeff, int) and coeff != 0:
                    cleaned[alpha] = coeff
                elif isinstance(coeff, float) and abs(coeff) > 1e-15:
                    cleaned[alpha] = coeff
            if cleaned:
                normalized[target_var] = cleaned

        object.__setattr__(self, "terms", normalized)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, DerivationSpec):
            return NotImplemented
        if self.n != other.n:
            return False
        if set(self.terms) != set(other.terms):
            return False
        for axis, mono_dict in self.terms.items():
            other_dict = other.terms.get(axis, {})
            if set(mono_dict) != set(other_dict):
                return False
            for alpha, coeff in mono_dict.items():
                other_coeff = other_dict[alpha]
                # Сравнение с допуском для float
                if isinstance(coeff, float) or isinstance(other_coeff, float):
                    if abs(float(coeff) - float(other_coeff)) > 1e-12:
                        return False
                else:
                    if coeff != other_coeff:
                        return False
        return True

    def __hash__(self) -> int:
        # Замороженный dataclass требует хеш; делаем по n и осям
        return hash((self.n, tuple(sorted(self.terms.keys()))))


# ---------------------------------------------------------------------------
# Примитивные кирпичики
# ---------------------------------------------------------------------------


def monomial_spec(
    n: int,
    axis: int,
    alpha: tuple[int, ...],
    coeff: Scalar = 1,
) -> DerivationSpec:
    """
    Мономиальное поле: coeff * z^alpha * ∂/∂z_axis.

    Args:
        n:     размерность пространства
        axis:  индекс компоненты производной (0..n-1)
        alpha: мультиидекс (кортеж длины n)
        coeff: коэффициент (по умолчанию 1)

    Пример: monomial_spec(2, 0, (2, 1)) = x^2*y*∂/∂x
    """
    if len(alpha) != n:
        raise ValueError(f"len(alpha)={len(alpha)} != n={n}")
    return DerivationSpec(n=n, terms={axis: {alpha: coeff}})


def partial_spec(n: int, axis: int) -> DerivationSpec:
    """
    Частная производная ∂/∂z_axis.

    Эквивалентно monomial_spec(n, axis, (0,...,0), 1).

    Пример: partial_spec(2, 0) = ∂/∂x
    """
    zero_alpha = (0,) * n
    return DerivationSpec(n=n, terms={axis: {zero_alpha: 1}})


def combine_specs(*specs: DerivationSpec) -> DerivationSpec:
    """
    Сумма дифференцирований: сложить все specs в одно.

    Коэффициенты одинаковых мономов суммируются. Нулевые результаты удаляются.

    Все specs должны иметь одинаковое n.
    """
    if not specs:
        raise ValueError("combine_specs требует хотя бы одного аргумента")

    n = specs[0].n
    for spec in specs[1:]:
        if spec.n != n:
            raise ValueError(
                f"Все specs должны иметь одинаковое n, "
                f"получено {spec.n} вместо {n}"
            )

    merged: dict[int, dict[tuple[int, ...], Scalar]] = {}
    for spec in specs:
        for axis, mono_dict in spec.terms.items():
            if axis not in merged:
                merged[axis] = {}
            for alpha, coeff in mono_dict.items():
                prev = merged[axis].get(alpha, 0)
                merged[axis][alpha] = prev + coeff  # type: ignore[operator]

    return DerivationSpec(n=n, terms=merged)


def spec_degree(spec: DerivationSpec) -> int:
    """
    Максимальная степень дифференцирования.

    Степень монома z^alpha * ∂/∂z_i равна sum(alpha) - 1
    (стандартное градуирование алгебры W_n).

    Возвращает -1 для нулевого (пустого) дифференцирования.
    """
    max_deg = -1
    for mono_dict in spec.terms.values():
        for alpha in mono_dict:
            deg = sum(alpha) - 1
            if deg > max_deg:
                max_deg = deg
    return max_deg


def random_spec(
    n: int,
    max_degree: int,
    sparsity: float = 0.5,
    seed: int | None = None,
) -> DerivationSpec:
    """
    Случайное разреженное дифференцирование из W_n до степени max_degree.

    Args:
        n:          размерность
        max_degree: максимальная степень дифференцирования (sum(alpha)-1 <= max_degree)
        sparsity:   вероятность включить каждый моном (0..1)
        seed:       seed для воспроизводимости

    Использует только стандартную библиотеку Python (без NumPy).
    """
    import itertools

    rng = _random.Random(seed)
    terms: dict[int, dict[tuple[int, ...], Scalar]] = {}

    # Перебираем все мономы z^alpha с sum(alpha) в диапазоне 0..max_degree+1
    # (степень дифференцирования = sum(alpha) - 1, поэтому берём до max_degree+1)
    for axis in range(n):
        mono_dict: dict[tuple[int, ...], Scalar] = {}
        for total in range(0, max_degree + 2):
            # Все мультиидексы длины n с суммой total
            for alpha in _multiindices(n, total):
                if rng.random() < sparsity:
                    coeff = rng.gauss(0.0, 1.0)
                    if abs(coeff) > 1e-15:
                        mono_dict[alpha] = coeff
        if mono_dict:
            terms[axis] = mono_dict

    if not terms:
        # Гарантируем непустой результат: берём одну случайную частную производную
        axis = rng.randrange(n)
        terms = {axis: {(0,) * n: 1.0}}

    return DerivationSpec(n=n, terms=terms)


# ---------------------------------------------------------------------------
# Скобка Ли и присоединённое действие
# ---------------------------------------------------------------------------


def _poly_diff(poly: dict, i: int) -> dict:
    """Дифференцировать полином по переменной i."""
    result: dict[tuple[int, ...], Scalar] = {}
    for alpha, c in poly.items():
        if alpha[i] > 0:
            new_alpha = tuple(a - (1 if k == i else 0) for k, a in enumerate(alpha))
            new_c = c * alpha[i]  # type: ignore[operator]
            prev = result.get(new_alpha, 0)
            result[new_alpha] = prev + new_c  # type: ignore[operator]
    return {a: c for a, c in result.items() if abs(float(c)) > 1e-15}


def _poly_mul(f: dict, g: dict) -> dict:
    """Перемножить два полинома."""
    result: dict[tuple[int, ...], Scalar] = {}
    for alpha, cf in f.items():
        for beta, cg in g.items():
            gamma = tuple(a + b for a, b in zip(alpha, beta))
            prev = result.get(gamma, 0)
            result[gamma] = prev + cf * cg  # type: ignore[operator]
    return {a: c for a, c in result.items() if abs(float(c)) > 1e-15}


def _spec_apply(spec: DerivationSpec, poly: dict) -> dict:
    """Применить деривацию к полиному: D(poly) = sum_i f_i * ∂(poly)/∂x_i."""
    result: dict[tuple[int, ...], Scalar] = {}
    for axis, f_i in spec.terms.items():
        dg = _poly_diff(poly, axis)
        if not dg:
            continue
        prod = _poly_mul(f_i, dg)
        for alpha, c in prod.items():
            prev = result.get(alpha, 0)
            result[alpha] = prev + c  # type: ignore[operator]
    return {a: c for a, c in result.items() if abs(float(c)) > 1e-15}


def lie_bracket(
    spec1: DerivationSpec,
    spec2: DerivationSpec,
    D_max: int | None = None,
) -> DerivationSpec:
    """
    Вычислить скобку Ли [spec1, spec2].

    Формула: ([D1, D2])_k = D1(g_k) - D2(f_k),
    где f_k, g_k — k-е компоненты spec1, spec2.

    Args:
        spec1: первое дифференцирование
        spec2: второе дифференцирование
        D_max: если задан, отбросить мономы степени > D_max

    Оба spec должны иметь одинаковое n.
    """
    if spec1.n != spec2.n:
        raise ValueError(
            f"lie_bracket требует одинакового n: {spec1.n} != {spec2.n}"
        )
    n = spec1.n
    terms: dict[int, dict[tuple[int, ...], Scalar]] = {}

    for k in range(n):
        g_k = spec2.terms.get(k, {})
        f_k = spec1.terms.get(k, {})

        contrib: dict[tuple[int, ...], Scalar] = {}
        for alpha, c in _spec_apply(spec1, g_k).items():
            prev = contrib.get(alpha, 0)
            contrib[alpha] = prev + c  # type: ignore[operator]
        for alpha, c in _spec_apply(spec2, f_k).items():
            prev = contrib.get(alpha, 0)
            contrib[alpha] = prev - c  # type: ignore[operator]

        clean: dict[tuple[int, ...], Scalar] = {
            a: c for a, c in contrib.items() if abs(float(c)) > 1e-15
        }
        if D_max is not None:
            clean = {a: c for a, c in clean.items() if sum(a) - 1 <= D_max}

        if clean:
            terms[k] = clean

    return DerivationSpec(n=n, terms=terms)


def ad(spec: DerivationSpec):
    """
    Присоединённое действие: ad(D)(E) = [D, E].

    Возвращает функцию D_max=None -> DerivationSpec.

    Пример:
        ad_D = ad(D)
        bracket = ad_D(E)           # [D, E]
        bracket_trunc = ad_D(E, D_max=5)
    """

    def _ad(other: DerivationSpec, D_max: int | None = None) -> DerivationSpec:
        return lie_bracket(spec, other, D_max=D_max)

    return _ad


# ---------------------------------------------------------------------------
# Конвертация в старый numeric формат
# ---------------------------------------------------------------------------


def to_generator(spec: DerivationSpec) -> "dict[int, dict[tuple[int, ...], float]]":
    """
    Конвертация DerivationSpec в формат Generator (dict[int, dict[tuple, float]]).

    Это обратно совместимый формат, используемый в lib/numeric/.
    Коэффициенты приводятся к float.
    """
    return {
        axis: {alpha: float(coeff) for alpha, coeff in mono_dict.items()}
        for axis, mono_dict in spec.terms.items()
    }


# ---------------------------------------------------------------------------
# Вспомогательные функции
# ---------------------------------------------------------------------------


def _multiindices(n: int, total: int):
    """Генератор всех мультиидексов длины n с суммой total."""
    import itertools

    if n == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _multiindices(n - 1, total - first):
            yield (first,) + rest
