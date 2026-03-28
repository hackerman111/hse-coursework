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


@dataclass(frozen=True)
class LeadingTerm:
    """
    Старший член в каноническом spec-представлении.

    Поле `degree` дублируется явно для удобства инспекции и сортировки.
    """

    axis: int
    alpha: tuple[int, ...]
    coeff: Scalar
    degree: int


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

    На вход принимает размерность `n`, номер компоненты `axis`,
    мультииндекс `alpha` и коэффициент `coeff`.
    На выходе возвращает `DerivationSpec` с одним мономом.

    >>> monomial_spec(2, axis=0, alpha=(2, 1), coeff=3)
    DerivationSpec(n=2, terms={0: {(2, 1): 3}})
    """
    if len(alpha) != n:
        raise ValueError(f"len(alpha)={len(alpha)} != n={n}")
    return DerivationSpec(n=n, terms={axis: {alpha: coeff}})


def partial_spec(n: int, axis: int) -> DerivationSpec:
    """
    Частная производная ∂/∂z_axis.

    На вход принимает размерность `n` и номер координаты `axis`.
    На выходе возвращает `DerivationSpec`, соответствующий `∂/∂z_axis`.

    >>> partial_spec(2, axis=1)
    DerivationSpec(n=2, terms={1: {(0, 0): 1}})
    """
    zero_alpha = (0,) * n
    return DerivationSpec(n=n, terms={axis: {zero_alpha: 1}})


def combine_specs(*specs: DerivationSpec) -> DerivationSpec:
    """
    Сумма дифференцирований: сложить все specs в одно.

    На вход принимает один или несколько `DerivationSpec` с одинаковой
    размерностью `n`.
    На выходе возвращает их сумму с объединёнными компонентами.

    >>> combine_specs(partial_spec(2, 0), partial_spec(2, 1))
    DerivationSpec(n=2, terms={0: {(0, 0): 1}, 1: {(0, 0): 1}})
    """
    if not specs:
        raise ValueError("combine_specs требует хотя бы одного аргумента")

    n = specs[0].n
    for spec in specs[1:]:
        if spec.n != n:
            raise ValueError(
                f"Все specs должны иметь одинаковое n, получено {spec.n} вместо {n}"
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

    На вход принимает один `DerivationSpec`.
    На выходе возвращает максимальную степень `sum(alpha) - 1`,
    а для нулевого spec возвращает `-1`.

    >>> spec_degree(monomial_spec(2, axis=0, alpha=(2, 1)))
    2
    """
    max_deg = -1
    for mono_dict in spec.terms.values():
        for alpha in mono_dict:
            deg = sum(alpha) - 1
            if deg > max_deg:
                max_deg = deg
    return max_deg


def zero_spec(n: int) -> DerivationSpec:
    """
    Нулевое дифференцирование в W_n.

    На вход принимает только размерность `n`.
    На выходе возвращает пустой `DerivationSpec` без мономов.

    >>> zero_spec(2)
    DerivationSpec(n=2, terms={})
    """
    return DerivationSpec(n=n)


def scale_spec(spec: DerivationSpec, scalar: Scalar) -> DerivationSpec:
    """
    Умножить дифференцирование на скаляр.

    На вход принимает `DerivationSpec` и число `scalar`.
    На выходе возвращает новый `DerivationSpec` с масштабированными коэффициентами.
    >>> scale_spec(partial_spec(2, 0), -2)
    DerivationSpec(n=2, terms={0: {(0, 0): -2}})
    """
    if scalar == 0:
        return zero_spec(spec.n)
    return DerivationSpec(
        n=spec.n,
        terms={
            axis: {
                alpha: coeff * scalar  # type: ignore[operator]
                for alpha, coeff in mono_dict.items()
            }
            for axis, mono_dict in spec.terms.items()
        },
    )


def normalize_spec(spec: DerivationSpec, eps: float = 1e-15) -> DerivationSpec:
    """
    Нормализовать spec, удаляя коэффициенты с модулем <= eps.

    На вход принимает `DerivationSpec` и порог `eps`.
    На выходе возвращает новый `DerivationSpec` без достаточно малых коэффициентов.

    >>> normalize_spec(monomial_spec(2, axis=0, alpha=(1, 0), coeff=1e-10), eps=1e-5)
    DerivationSpec(n=2, terms={})
    """
    return DerivationSpec(
        n=spec.n,
        terms={
            axis: {
                alpha: coeff
                for alpha, coeff in mono_dict.items()
                if abs(float(coeff)) > eps
            }
            for axis, mono_dict in spec.terms.items()
        },
    )


def homogeneous_component(spec: DerivationSpec, degree: int) -> DerivationSpec:
    """
    Извлечь однородную компоненту фиксированной степени.

    На вход принимает `DerivationSpec` и целую степень `degree`.
    На выходе возвращает только те мономы, у которых `sum(alpha) - 1 == degree`.

    >>> source = combine_specs(partial_spec(2, 0), monomial_spec(2, 1, (2, 0), coeff=5))
    >>> homogeneous_component(source, -1)
    DerivationSpec(n=2, terms={0: {(0, 0): 1}})
    """
    return DerivationSpec(
        n=spec.n,
        terms={
            axis: {
                alpha: coeff
                for alpha, coeff in mono_dict.items()
                if sum(alpha) - 1 == degree
            }
            for axis, mono_dict in spec.terms.items()
        },
    )


def homogeneous_components(spec: DerivationSpec) -> dict[int, DerivationSpec]:
    """
    Разбить spec на однородные компоненты по стандартному градуированию W_n.

    На вход принимает один `DerivationSpec`.
    На выходе возвращает словарь `степень -> однородная компонента`.

    >>> source = combine_specs(partial_spec(2, 0), monomial_spec(2, 1, (2, 0), coeff=5))
    >>> parts = homogeneous_components(source)
    >>> sorted(parts)
    [-1, 1]
    >>> parts[-1] == partial_spec(2, 0)
    True
    >>> parts[1] == monomial_spec(2, 1, (2, 0), coeff=5)
    True
    """
    grouped: dict[int, dict[int, dict[tuple[int, ...], Scalar]]] = {}
    for axis, mono_dict in spec.terms.items():
        for alpha, coeff in mono_dict.items():
            degree = sum(alpha) - 1
            grouped.setdefault(degree, {}).setdefault(axis, {})[alpha] = coeff

    return {
        degree: DerivationSpec(n=spec.n, terms=terms)
        for degree, terms in grouped.items()
    }


def leading_term_spec(spec: DerivationSpec, order: str = "lex") -> LeadingTerm | None:
    """
    Выбрать старший член spec по заданному мономиальному порядку.

    На вход принимает `DerivationSpec` и имя порядка `order`.
    На выходе возвращает `LeadingTerm` для старшего монома или `None` для нуля.

    >>> source = combine_specs(partial_spec(2, 0), monomial_spec(2, 1, (2, 0), coeff=5))
    >>> leading_term_spec(source)
    LeadingTerm(axis=1, alpha=(2, 0), coeff=5, degree=1)
    """
    if order != "lex":
        raise ValueError(
            f"Неподдерживаемый порядок {order!r}; сейчас доступен только 'lex'"
        )

    best: LeadingTerm | None = None
    for axis, mono_dict in spec.terms.items():
        for alpha, coeff in mono_dict.items():
            current = LeadingTerm(
                axis=axis,
                alpha=alpha,
                coeff=coeff,
                degree=sum(alpha) - 1,
            )
            if best is None:
                best = current
                continue
            if current.alpha > best.alpha:
                best = current
                continue
            if current.alpha == best.alpha and current.axis > best.axis:
                best = current

    return best


def random_spec(
    n: int,
    max_degree: int,
    sparsity: float = 0.5,
    seed: int | None = None,
) -> DerivationSpec:
    """
    Случайное разреженное дифференцирование из W_n до степени max_degree.

    На вход принимает размерность `n`, верхнюю границу `max_degree`,
    разреженность `sparsity` и необязательный `seed`.
    На выходе возвращает воспроизводимый случайный `DerivationSpec`.

    >>> random_spec(2, max_degree=1, sparsity=0.0, seed=7)
    DerivationSpec(n=2, terms={1: {(0, 0): 1.0}})
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
