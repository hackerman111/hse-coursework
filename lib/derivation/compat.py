"""
Совместимость со старым API.
"""

from __future__ import annotations

from typing import Any

from .base import LieDerivation


def Derivation(algebra: Any, gen_mapping: Any = None) -> LieDerivation:
    """
    Legacy-обёртка над ``LieDerivation.from_mapping`` для совместимости со старым API.

    На вход принимает Sage-алгебру ``algebra`` и отображение образов генераторов ``gen_mapping``.
    На выходе возвращает объект ``LieDerivation``.

    >>> from sage.all import PolynomialRing, QQ
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> type(Derivation(ring, {ring.gen(0): 1, ring.gen(1): 0})).__name__
    'LieDerivation'
    """
    return LieDerivation.from_mapping(algebra, gen_mapping)
