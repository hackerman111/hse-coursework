"""
Re-export shim: функции разбора и форматирования перенесены в lib/io/parse.py.
Конструкторы генераторов (make_*) остаются здесь как legacy API.

Оставлен для обратной совместимости.
"""

from __future__ import annotations

import numpy as np

# Re-exports из lib/io/parse.py
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
    """Create ``coeff * z^alpha * ∂_{z_target_var}``."""
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
    Legacy helper для ``x^a y^b d(target_var)``.

    Для n > 2 только первые две координаты получают ненулевые степени.
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
    """Create ``z_0^e0 * ... * z_{n-1}^e_{n-1} * ∂_{z_target_var}``."""
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
    """Случайный разреженный генератор в W_n до указанной степени."""
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
    """Legacy helper: случайный второй генератор."""
    rng = np.random.default_rng(seed)
    return make_random_generator(n, max_degree=max_degree, rng=rng, sparsity=0.3)


def resolve_hypothesis_generator(
    n: int,
    spec: str | None,
    a: int,
    b: int,
) -> Generator:
    """Возвращает генератор из строковой спецификации или из (a, b)."""
    if spec is not None:
        return parse_generator_spec(spec, n)
    return make_monomial_D(n, a, b, target_var=0)
