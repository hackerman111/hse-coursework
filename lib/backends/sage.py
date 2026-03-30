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


def _representative_for_export(derivation: Any, element: Any) -> Any:
    """
    Получить полиномиального представителя элемента для экспорта в DerivationSpec.
    """
    polynomial_ops = getattr(derivation, "_polynomial_ops", None)
    if polynomial_ops is not None and hasattr(polynomial_ops, "representative"):
        return polynomial_ops.representative(element)
    if hasattr(element, "lift"):
        return element.lift()
    return element


def to_sage(spec: DerivationSpec, algebra: Any) -> Any:
    """
    Конвертировать каноническое spec-представление в объект SageMath-деривации.

    На вход принимает ``DerivationSpec`` и полиномиальную алгебру Sage ``algebra`` той же размерности.
    На выходе возвращает объект ``LieDerivation``, согласованный с данным ``spec``.

    >>> from sage.all import PolynomialRing, QQ
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> derivation = to_sage(monomial_spec(2, 0, (1, 0)), ring)
    >>> type(derivation).__name__
    'LieDerivation'
    """
    from lib.derivation.base import LieDerivation

    gens = algebra.gens()
    n = len(gens)

    if spec.n != n:
        raise ValueError(
            f"spec.n={spec.n} не совпадает с числом генераторов алгебры {n}"
        )

    # Строим отображение: gen_i → образ в виде многочлена из algebra
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
    Экспортировать SageMath-деривацию в каноническое ``DerivationSpec``.

    На вход принимает объект ``LieDerivation``.
    На выходе возвращает ``DerivationSpec``, восстановленный по образам координатных функций.

    >>> from sage.all import PolynomialRing, QQ
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> derivation = to_sage(combine_specs(monomial_spec(2, 0, (1, 0)), monomial_spec(2, 1, (0, 0))), ring)
    >>> from_sage(derivation)
    DerivationSpec(n=2, terms={0: {(1, 0): 1.0}, 1: {(0, 0): 1.0}})
    """
    algebra = derivation.codomain()
    gens = algebra.gens()
    n = len(gens)

    terms: dict[int, dict[tuple[int, ...], float]] = {}

    for axis, gen in enumerate(gens):
        poly = _representative_for_export(derivation, derivation(gen))
        if poly == 0:
            continue

        mono_dict: dict[tuple[int, ...], float] = {}

        try:
            # Итерируем по мономам многочлена
            for monomial, coeff in poly.dict().items():
                if hasattr(monomial, "__iter__"):
                    alpha = tuple(int(e) for e in monomial)
                else:
                    # Для одной переменной monomial хранится как целая степень
                    alpha = (int(monomial),)
                c = float(coeff)
                if abs(c) > 1e-15:
                    mono_dict[alpha] = c
        except AttributeError:
            # Если poly является константой, обрабатываем её отдельно
            c = float(poly)
            if abs(c) > 1e-15:
                mono_dict[(0,) * n] = c

        if mono_dict:
            terms[axis] = mono_dict

    return DerivationSpec(n=n, terms=terms)
