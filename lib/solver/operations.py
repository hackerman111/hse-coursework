"""
Вспомогательные операции solver.
"""

from __future__ import annotations

from typing import Dict, Iterable, Tuple

from ..derivation import LieDerivation


def reduce_derivation(
    derivation: LieDerivation,
    basis: Dict[Tuple[int, object], LieDerivation],
) -> LieDerivation:
    """
    Редуцировать дифференцирование по уже построенному базису старших членов.

    На вход принимает ``LieDerivation`` и словарь базиса ``(axis, monomial) -> LieDerivation``.
    На выходе возвращает остаток после последовательного вычитания базисных элементов с тем же старшим членом.

    >>> from sage.all import PolynomialRing, QQ
    >>> from lib.derivation import LieDerivation
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> d = LieDerivation.from_mapping(ring, [ring.gen(0), 0])
    >>> reduce_derivation(d, {})
    LieDerivation([x, 0])
    """
    while True:
        leading_term = derivation.leading_term
        if leading_term is None:
            return derivation

        comp_idx, monomial, coeff = leading_term
        basis_key = (comp_idx, monomial)
        if basis_key not in basis:
            return derivation

        derivation = derivation - (basis[basis_key] * coeff)


def generate_commutators(
    new_derivation: LieDerivation,
    processed: Iterable[LieDerivation],
) -> list[LieDerivation]:
    """
    Построить все скобки Ли нового элемента с уже обработанными.

    На вход принимает новое дифференцирование ``new_derivation`` и итерируемый набор ``processed``.
    На выходе возвращает список скобок ``[new_derivation, old]``.

    >>> from sage.all import PolynomialRing, QQ
    >>> from lib.derivation import LieDerivation
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> d = LieDerivation.from_mapping(ring, [1, 0])
    >>> generate_commutators(d, [d])
    [LieDerivation([0, 0])]
    """
    return [new_derivation.bracket(old_derivation) for old_derivation in processed]
