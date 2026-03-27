"""
Разбор и форматирование генераторов алгебры Ли W_n.

Это библиотечный модуль — все функции парсинга/форматирования из CLI
перенесены сюда. CLI-скрипты должны использовать этот модуль напрямую.

Поддерживаемый формат строк:
    "x^5*y^5*dx"         = x^5 y^5 d/dx
    "dx + x*dy"          = d/dx + x*d/dy
    "3*x^2*y*dx - y^3*dy" = 3x^2 y d/dx - y^3 d/dy

Переменные: x,y,z,u,v,w (по позиции, индексы 0..5).
Производные: dx, dy, dz, ... или d0, d1, ... или dz0, dz1, ...
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from lib.core.spec import DerivationSpec

if TYPE_CHECKING:
    from lib.numeric.types import Generator


# Алиасы переменных: x→0, y→1, z→2, u→3, v→4, w→5
VARIABLE_ALIASES = ("x", "y", "z", "u", "v", "w")


def generator_term_count(generator: "Generator") -> int:
    """Количество ненулевых мономов во всех компонентах генератора."""
    return sum(len(monomials) for monomials in generator.values())


def generator_max_degree(generator: "Generator") -> int:
    """Максимальная однородная степень в генераторе."""
    max_degree = -1
    for monomials in generator.values():
        for alpha in monomials:
            max_degree = max(max_degree, sum(alpha) - 1)
    return max_degree


def describe_generator(
    generator: "Generator",
    n: int,
    max_terms: int | None = None,
) -> str:
    """Человекочитаемое описание разреженного генератора."""
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


def describe_random_generator(generator: "Generator", n: int) -> str:
    """Компактное описание случайного генератора."""
    return (
        f"random terms={generator_term_count(generator)}, "
        f"max_degree={generator_max_degree(generator)}, "
        f"targets={sorted(generator.keys())}, "
        f"example={describe_generator(generator, n, max_terms=4)}"
    )


def parse_generator_spec(spec: str, n: int) -> "Generator":
    """
    Разбор строк вида ``x^5*y^5*dx`` или ``dx + x*dy``.

    Returns:
        Generator: словарь вида {axis: {alpha: coeff}}
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
        generator[target_var][alpha_tuple] = (
            generator[target_var].get(alpha_tuple, 0.0) + coeff
        )
        if abs(generator[target_var][alpha_tuple]) <= 1e-15:
            del generator[target_var][alpha_tuple]
        if not generator[target_var]:
            del generator[target_var]

    if not generator:
        raise ValueError(
            f"генератор {spec!r} занулился после сокращения коэффициентов"
        )
    return generator


def spec_from_string(text: str, n: int) -> DerivationSpec:
    """
    Разбор строки вида ``x^5*y^5*dx`` в DerivationSpec.

    Обёртка над parse_generator_spec, возвращающая каноническую модель.
    """
    generator = parse_generator_spec(text, n)
    return DerivationSpec(n=n, terms=generator)


def spec_to_string(
    spec: DerivationSpec,
    variables: list[str] | None = None,
) -> str:
    """
    Каноническое строковое представление DerivationSpec.

    Args:
        spec:      дифференцирование
        variables: имена переменных (по умолчанию x,y,z,u,v,w или z0,z1,...)

    Returns:
        Строка вида ``x^2*y*dx + 3*dy`` (совместима с spec_from_string).
    """
    n = spec.n
    if variables is None:
        if n <= len(VARIABLE_ALIASES):
            variables = list(VARIABLE_ALIASES[:n])
        else:
            variables = [f"z{i}" for i in range(n)]

    deriv_names = [f"d{v}" for v in variables]

    parts: list[str] = []
    for axis in sorted(spec.terms.keys()):
        mono_dict = spec.terms[axis]
        for alpha in sorted(mono_dict.keys()):
            coeff = mono_dict[alpha]
            # Форматируем коэффициент
            if isinstance(coeff, float):
                coeff_str = f"{coeff:g}"
            else:
                coeff_str = str(coeff)

            # Форматируем моном z^alpha
            monomial_factors: list[str] = []
            for i, exp in enumerate(alpha):
                if exp == 1:
                    monomial_factors.append(variables[i])
                elif exp > 1:
                    monomial_factors.append(f"{variables[i]}^{exp}")

            monomial = "*".join(monomial_factors)
            deriv = deriv_names[axis]

            if monomial:
                if coeff_str == "1":
                    parts.append(f"{monomial}*{deriv}")
                elif coeff_str == "-1":
                    parts.append(f"-{monomial}*{deriv}")
                else:
                    parts.append(f"{coeff_str}*{monomial}*{deriv}")
            else:
                if coeff_str == "1":
                    parts.append(deriv)
                elif coeff_str == "-1":
                    parts.append(f"-{deriv}")
                else:
                    parts.append(f"{coeff_str}*{deriv}")

    if not parts:
        return "0"

    # Собираем строку: первый член без знака, остальные со знаком
    result_parts: list[str] = []
    for i, part in enumerate(parts):
        if i == 0:
            result_parts.append(part)
        else:
            if part.startswith("-"):
                result_parts.append(f" - {part[1:]}")
            else:
                result_parts.append(f" + {part}")

    return "".join(result_parts)


# ---------------------------------------------------------------------------
# Вспомогательные функции разбора (внутренние)
# ---------------------------------------------------------------------------


def _parse_derivative_factor(token: str, n: int) -> int | None:
    """Распознать токен как производную d0, dz1, dx и т.п."""
    if token.startswith("d") and token[1:].isdigit():
        return _validate_axis(int(token[1:]), n, token)

    if token.startswith("dz") and token[2:].isdigit():
        return _validate_axis(int(token[2:]), n, token)

    alias_index = _alias_index(token[1:]) if token.startswith("d") else None
    if alias_index is not None:
        return _validate_axis(alias_index, n, token)

    return None


def _parse_variable_factor(token: str, n: int) -> tuple[int | None, int]:
    """Распознать токен как переменную x^k, z0^3 и т.п."""
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
