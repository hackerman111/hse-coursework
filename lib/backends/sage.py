"""
Sage backend: конвертация DerivationSpec ↔ LieDerivation (SageMath).

ВАЖНО: этот модуль является единственным местом в библиотеке, где используется
import SageMath. Весь остальной код (lib/core/, lib/io/, lib/generators/)
работает без Sage.

Импорт этого модуля требует установленного SageMath.
"""

from __future__ import annotations

from typing import Any

from lib.core.spec import DerivationSpec, monomial_spec, combine_specs


def to_sage(spec: DerivationSpec, algebra: Any) -> Any:
    """
    Конвертация DerivationSpec в LieDerivation (SageMath).

    Args:
        spec:    каноническое представление дифференцирования
        algebra: кольцо многочленов SageMath (PolynomialRing)

    Returns:
        LieDerivation — объект из lib/derivation/
    """
    from lib.derivation.base import LieDerivation

    gens = algebra.gens()
    n = len(gens)

    if spec.n != n:
        raise ValueError(
            f"spec.n={spec.n} не совпадает с числом генераторов алгебры {n}"
        )

    # Строим отображение: gen_i → образ (многочлен из algebra)
    images = [algebra.zero() for _ in range(n)]

    for axis, mono_dict in spec.terms.items():
        for alpha, coeff in mono_dict.items():
            # Строим моном z^alpha
            monomial = algebra.one()
            for i, exp in enumerate(alpha):
                if exp > 0:
                    monomial *= gens[i] ** exp
            images[axis] += algebra(coeff) * monomial

    return LieDerivation.from_mapping(algebra, images)


def from_sage(derivation: Any) -> DerivationSpec:
    """
    Конвертация LieDerivation (SageMath) в DerivationSpec.

    Извлекает коэффициенты многочленов из образов дифференцирования
    и строит каноническое представление.

    Args:
        derivation: LieDerivation — дифференцирование из lib/derivation/

    Returns:
        DerivationSpec — каноническое представление
    """
    algebra = derivation.codomain()
    gens = algebra.gens()
    n = len(gens)

    terms: dict[int, dict[tuple[int, ...], float]] = {}

    for axis, gen in enumerate(gens):
        poly = derivation(gen)
        if poly == 0:
            continue

        mono_dict: dict[tuple[int, ...], float] = {}

        try:
            # Итерируем по мономам многочлена
            for monomial, coeff in poly.dict().items():
                if hasattr(monomial, "__iter__"):
                    alpha = tuple(int(e) for e in monomial)
                else:
                    # Одна переменная: monomial — целое число (степень)
                    alpha = (int(monomial),)
                c = float(coeff)
                if abs(c) > 1e-15:
                    mono_dict[alpha] = c
        except AttributeError:
            # poly — константа
            c = float(poly)
            if abs(c) > 1e-15:
                mono_dict[(0,) * n] = c

        if mono_dict:
            terms[axis] = mono_dict

    return DerivationSpec(n=n, terms=terms)
