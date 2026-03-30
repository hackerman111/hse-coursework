"""
Re-export shim: функции разбора и форматирования перенесены в lib/io/parse.py.
Конструкторы генераторов (make_*) остаются здесь как legacy API.

Оставлен для обратной совместимости.
"""

from __future__ import annotations

import numpy as np

# Реэкспорт из lib/io/parse.py
from lib.io.parse import (
    VARIABLE_ALIASES,
    describe_generator,
    describe_random_generator,
    generator_max_degree,
    generator_term_count,
    parse_generator_spec,
)

from .indexing import multiindices
from .types import Generator

__all__ = [
    "VARIABLE_ALIASES",
    "describe_generator",
    "describe_random_generator",
    "generator_max_degree",
    "generator_term_count",
    "parse_generator_spec",
    "make_monomial_generator",
    "make_monomial_D",
    "make_general_monomial_generator",
    "make_random_generator",
    "make_random_E",
    "resolve_hypothesis_generator",
]


def make_monomial_generator(
    n: int,
    target_var: int,
    alpha: tuple[int, ...],
    coeff: float = 1.0,
) -> Generator:
    """
    Построить мономиальный генератор ``coeff * z^alpha * d/dz_target_var``.

    На вход принимает размерность ``n``, индекс компоненты ``target_var``, мультииндекс ``alpha`` и коэффициент ``coeff``.
    На выходе возвращает generator в разреженном словарном формате.

    >>> make_monomial_generator(2, target_var=1, alpha=(2, 0), coeff=3.0)
    {1: {(2, 0): 3.0}}
    """
    if len(alpha) != n:
        raise ValueError(f"expected multiindex of length {n}, got {len(alpha)}")
    return {target_var: {tuple(alpha): coeff}}


def make_monomial_D(
    n: int,
    a: int,
    b: int,
    target_var: int = 0,
) -> Generator:
    """
    Построить legacy-генератор вида ``x^a y^b d(target_var)``.

    На вход принимает размерность ``n``, показатели ``a`` и ``b`` и целевую компоненту ``target_var``.
    На выходе возвращает generator, где ненулевые степени по умолчанию сосредоточены в первых двух координатах.

    >>> make_monomial_D(2, 2, 1, target_var=0)
    {0: {(2, 1): 1.0}}
    """
    alpha = [0] * n
    if n >= 1:
        alpha[0] = a
    if n >= 2:
        alpha[1] = b
    return make_monomial_generator(n, target_var=target_var, alpha=tuple(alpha))


def make_general_monomial_generator(
    n: int,
    exponents: list[int],
    target_var: int = 0,
) -> Generator:
    """
    Построить мономиальный генератор с произвольным списком показателей.

    На вход принимает размерность ``n``, список ``exponents`` и индекс компоненты ``target_var``.
    На выходе возвращает generator для монома ``z_0^e0 ... z_{n-1}^e_{n-1} * d/dz_target_var``.

    >>> make_general_monomial_generator(3, [1, 0, 2], target_var=2)
    {2: {(1, 0, 2): 1.0}}
    """
    alpha = [0] * n
    for index, exponent in enumerate(exponents[:n]):
        alpha[index] = exponent
    return make_monomial_generator(n, target_var=target_var, alpha=tuple(alpha))


def make_random_generator(
    n: int,
    max_degree: int,
    rng: np.random.Generator | None = None,
    sparsity: float = 0.5,
) -> Generator:
    """
    Сгенерировать случайный разреженный элемент алгебры ``W_n``.

    На вход принимает размерность ``n``, максимальную степень ``max_degree``, генератор случайных чисел ``rng`` и параметр разреженности ``sparsity``.
    На выходе возвращает generator со случайными коэффициентами у мономов степени не выше ``max_degree``.

    >>> rng = np.random.default_rng(0)
    >>> sorted(make_random_generator(1, 0, rng=rng, sparsity=1.0).keys())
    [0]
    """
    if rng is None:
        rng = np.random.default_rng()

    generator: Generator = {}
    for axis in range(n):
        monomials: dict[tuple[int, ...], float] = {}
        for total_degree in range(0, max_degree + 2):
            for alpha in multiindices(n, total_degree):
                if rng.random() < sparsity:
                    monomials[alpha] = float(rng.standard_normal())
        if monomials:
            generator[axis] = monomials

    return generator


def make_random_E(n: int, max_degree: int, seed: int | None = None) -> Generator:
    """
    Построить legacy-случайный второй генератор для численных экспериментов.

    На вход принимает размерность ``n``, максимальную степень ``max_degree`` и необязательное зерно ``seed``.
    На выходе возвращает воспроизводимый случайный generator.

    >>> make_random_E(1, 0, seed=0) == make_random_E(1, 0, seed=0)
    True
    """
    rng = np.random.default_rng(seed)
    return make_random_generator(n, max_degree=max_degree, rng=rng, sparsity=0.3)


def resolve_hypothesis_generator(
    n: int,
    spec: str | None,
    a: int,
    b: int,
) -> Generator:
    """
    Разрешить описание второго генератора из явной строки или legacy-пары степеней.

    На вход принимает размерность ``n``, необязательную строку ``spec`` и показатели ``a, b``.
    На выходе возвращает generator, полученный либо парсером, либо из ``make_monomial_D``.

    >>> resolve_hypothesis_generator(2, "dx", 3, 4)
    {0: {(0, 0): 1.0}}
    """
    if spec is not None:
        return parse_generator_spec(spec, n)
    return make_monomial_D(n, a, b, target_var=0)
