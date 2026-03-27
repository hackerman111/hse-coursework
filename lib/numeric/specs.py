"""
Public builders, parser, and formatting helpers for numeric generators.
"""

from __future__ import annotations

import re

import numpy as np

from .indexing import multiindices
from .types import Generator


VARIABLE_ALIASES = ("x", "y", "z", "u", "v", "w")


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
    Legacy helper for ``x^a y^b d(target_var)``.

    For ``n > 2`` only the first two coordinates receive nonzero exponents.
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
    """
    Create a random sparse generator in ``W_n`` up to the requested degree.
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
    """Legacy helper returning a random second generator."""
    rng = np.random.default_rng(seed)
    return make_random_generator(n, max_degree=max_degree, rng=rng, sparsity=0.3)


def generator_term_count(generator: Generator) -> int:
    """Count nonzero monomials across all components."""
    return sum(len(monomials) for monomials in generator.values())


def generator_max_degree(generator: Generator) -> int:
    """Maximum homogeneous degree appearing in the generator."""
    max_degree = -1
    for monomials in generator.values():
        for alpha in monomials:
            max_degree = max(max_degree, sum(alpha) - 1)
    return max_degree


def describe_generator(
    generator: Generator,
    n: int,
    max_terms: int | None = None,
) -> str:
    """Human-readable description of a sparse generator."""
    variables = [f"z{i}" for i in range(n)]
    parts: list[str] = []
    for axis, monomials in sorted(generator.items()):
        for alpha, coeff in sorted(monomials.items()):
            factors = []
            for index, exponent in enumerate(alpha):
                if exponent == 1:
                    factors.append(variables[index])
                elif exponent > 1:
                    factors.append(f"{variables[index]}^{exponent}")
            monomial = "*".join(factors) if factors else "1"
            coeff_prefix = "" if coeff == 1.0 else f"{coeff:g}·"
            parts.append(f"{coeff_prefix}{monomial}·∂_{variables[axis]}")
    if max_terms is not None and len(parts) > max_terms:
        return " + ".join(parts[:max_terms]) + f" + ... ({len(parts)} terms)"
    return " + ".join(parts) if parts else "0"


def describe_random_generator(generator: Generator, n: int) -> str:
    """Compact summary tailored for randomly generated inputs."""
    return (
        f"random terms={generator_term_count(generator)}, "
        f"max_degree={generator_max_degree(generator)}, "
        f"targets={sorted(generator.keys())}, "
        f"example={describe_generator(generator, n, max_terms=4)}"
    )


def parse_generator_spec(spec: str, n: int) -> Generator:
    """
    Parse strings such as ``x^5*y^5*dx`` or ``dx + x*dy``.
    """
    text = spec.replace(" ", "")
    if not text:
        raise ValueError("пустая спецификация генератора")

    normalized = text if text[0] in "+-" else f"+{text}"
    terms = re.findall(r"[+-][^+-]+", normalized)
    if not terms:
        raise ValueError(f"не удалось разобрать генератор: {spec!r}")

    generator: Generator = {}
    for raw_term in terms:
        sign = -1.0 if raw_term[0] == "-" else 1.0
        body = raw_term[1:]
        factors = [factor for factor in body.split("*") if factor]
        if not factors:
            raise ValueError(f"пустой член в генераторе: {spec!r}")

        coeff = sign
        alpha = [0] * n
        target_var = None

        for factor in factors:
            if re.fullmatch(r"\d+(?:\.\d+)?", factor):
                coeff *= float(factor)
                continue

            deriv_index = _parse_derivative_factor(factor, n)
            if deriv_index is not None:
                if target_var is not None:
                    raise ValueError(
                        f"в одном члене больше одной производной: {raw_term!r}"
                    )
                target_var = deriv_index
                continue

            var_index, exponent = _parse_variable_factor(factor, n)
            if var_index is None:
                raise ValueError(f"неизвестный фактор {factor!r} в {raw_term!r}")
            if exponent > 0:
                alpha[var_index] += exponent

        if target_var is None:
            raise ValueError(f"в члене {raw_term!r} не указана производная")

        generator.setdefault(target_var, {})
        alpha_tuple = tuple(alpha)
        generator[target_var][alpha_tuple] = generator[target_var].get(alpha_tuple, 0.0) + coeff
        if abs(generator[target_var][alpha_tuple]) <= 1e-15:
            del generator[target_var][alpha_tuple]
        if not generator[target_var]:
            del generator[target_var]

    if not generator:
        raise ValueError(f"генератор {spec!r} занулился после сокращения коэффициентов")
    return generator


def _parse_derivative_factor(token: str, n: int) -> int | None:
    if token.startswith("d") and token[1:].isdigit():
        return _validate_axis(int(token[1:]), n, token)

    if token.startswith("dz") and token[2:].isdigit():
        return _validate_axis(int(token[2:]), n, token)

    alias_index = _alias_index(token[1:]) if token.startswith("d") else None
    if alias_index is not None:
        return _validate_axis(alias_index, n, token)

    return None


def _parse_variable_factor(token: str, n: int) -> tuple[int | None, int]:
    exponent = 1
    base = token
    if "^" in token:
        base, exponent_str = token.split("^", 1)
        if not exponent_str.isdigit():
            raise ValueError(f"степень должна быть целой неотрицательной: {token!r}")
        exponent = int(exponent_str)

    if base.startswith("z") and base[1:].isdigit():
        return _validate_axis(int(base[1:]), n, token), exponent

    alias_index = _alias_index(base)
    if alias_index is not None:
        return _validate_axis(alias_index, n, token), exponent

    if base == "1":
        return 0, 0
    return None, exponent


def _alias_index(name: str) -> int | None:
    try:
        return VARIABLE_ALIASES.index(name)
    except ValueError:
        return None


def _validate_axis(index: int, n: int, token: str) -> int:
    if index < 0 or index >= n:
        raise ValueError(f"индекс {index} вне диапазона 0..{n - 1} в {token!r}")
    return index

