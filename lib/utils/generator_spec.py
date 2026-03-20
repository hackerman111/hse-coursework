"""
Разбор и форматирование текстовых спецификаций генераторов.
"""

from __future__ import annotations

import re
from typing import Any

from ..derivation import LieDerivation


VARIABLE_ALIASES = ("x", "y", "z", "u", "v", "w")


def parse_generator_spec(spec: str, n: int) -> dict[int, dict[tuple[int, ...], float]]:
    """
    Разобрать строку вида ``dx + x^2*dy - 3*z0*d1`` в sparse-словарь.
    """
    text = spec.replace(" ", "")
    if not text:
        raise ValueError("пустая спецификация генератора")

    normalized = text if text[0] in "+-" else f"+{text}"
    terms = re.findall(r"[+-][^+-]+", normalized)
    if not terms:
        raise ValueError(f"не удалось разобрать генератор: {spec!r}")

    generator: dict[int, dict[tuple[int, ...], float]] = {}
    for raw_term in terms:
        sign = -1.0 if raw_term[0] == "-" else 1.0
        body = raw_term[1:]
        factors = [factor for factor in body.split("*") if factor]
        if not factors:
            raise ValueError(f"пустой член в генераторе: {spec!r}")

        coefficient = sign
        exponents = [0] * n
        target_var = None

        for factor in factors:
            if re.fullmatch(r"\d+(?:\.\d+)?", factor):
                coefficient *= float(factor)
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
            exponents[var_index] += exponent

        if target_var is None:
            raise ValueError(f"в члене {raw_term!r} не указана производная")

        generator.setdefault(target_var, {})
        key = tuple(exponents)
        generator[target_var][key] = generator[target_var].get(key, 0.0) + coefficient
        if abs(generator[target_var][key]) <= 1e-15:
            del generator[target_var][key]
        if not generator[target_var]:
            del generator[target_var]

    if not generator:
        raise ValueError(f"генератор {spec!r} занулился после сокращения коэффициентов")
    return generator


def mapping_to_derivation(algebra: Any, mapping: dict[int, dict[tuple[int, ...], float]]):
    gens = algebra.gens()
    n = len(gens)
    derivation_mapping = {gen: 0 for gen in gens}

    for component, monomials in mapping.items():
        value = 0
        for exponents, coefficient in monomials.items():
            key = exponents[0] if n == 1 else tuple(exponents)
            value += algebra({key: coefficient})
        derivation_mapping[gens[component]] = value

    return LieDerivation.from_mapping(algebra, derivation_mapping)


def derivation_to_mapping(
    derivation: LieDerivation,
) -> dict[int, dict[tuple[int, ...], Any]]:
    algebra = derivation.codomain()
    mapping: dict[int, dict[tuple[int, ...], Any]] = {}
    for component, gen in enumerate(algebra.gens()):
        poly = algebra(derivation(gen))
        if poly == 0:
            continue
        normalized = {}
        for exp, coeff in poly.dict().items():
            try:
                key = tuple(int(value) for value in exp)
            except TypeError:
                key = (int(exp),)
            normalized[key] = coeff
        mapping[component] = normalized
    return mapping


def describe_generator_mapping(
    mapping: dict[int, dict[tuple[int, ...], Any]],
    n: int,
    *,
    max_terms: int | None = None,
) -> str:
    vars_ = [f"z{i}" for i in range(n)]
    parts = []
    for component, monomials in sorted(mapping.items()):
        for exponents, coefficient in sorted(monomials.items()):
            factors = []
            for axis, power in enumerate(exponents):
                if power == 1:
                    factors.append(vars_[axis])
                elif power > 1:
                    factors.append(f"{vars_[axis]}^{power}")
            monomial = "*".join(factors) if factors else "1"
            coeff_text = "" if coefficient == 1 else f"{coefficient}*"
            parts.append(f"{coeff_text}{monomial}*d{component}")
    if max_terms is not None and len(parts) > max_terms:
        return " + ".join(parts[:max_terms]) + f" + ... ({len(parts)} terms)"
    return " + ".join(parts) if parts else "0"


def generator_term_count(mapping: dict[int, dict[tuple[int, ...], Any]]) -> int:
    return sum(len(monomials) for monomials in mapping.values())


def generator_coeff_degree(mapping: dict[int, dict[tuple[int, ...], Any]]) -> int:
    max_degree = -1
    for monomials in mapping.values():
        for exponents in monomials:
            max_degree = max(max_degree, sum(exponents))
    return max_degree


def _parse_derivative_factor(token: str, n: int):
    if token.startswith("d") and token[1:].isdigit():
        return _validate_axis(int(token[1:]), n, token)

    if token.startswith("dz") and token[2:].isdigit():
        return _validate_axis(int(token[2:]), n, token)

    alias_index = _alias_index(token[1:]) if token.startswith("d") else None
    if alias_index is not None:
        return _validate_axis(alias_index, n, token)

    return None


def _parse_variable_factor(token: str, n: int):
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


def _alias_index(name: str):
    try:
        return VARIABLE_ALIASES.index(name)
    except ValueError:
        return None


def _validate_axis(index: int, n: int, token: str) -> int:
    if index < 0 or index >= n:
        raise ValueError(f"индекс {index} вне диапазона 0..{n - 1} в {token!r}")
    return index
