"""
Проверка полноты набора генераторов.
"""

from __future__ import annotations

from typing import Iterable

from ..solver import LieBasisSolver


def check(generators: Iterable, max_iter: int = 1000) -> bool:
    """
    Проверить, порождают ли генераторы всю алгебру Ли.

    На вход принимает итерируемый набор `LieDerivation` и лимит `max_iter`.
    На выходе возвращает `True`, если замыкание по скобке Ли достигло цели.

    >>> from sage.all import PolynomialRing, QQ
    >>> from lib.generator import get_Beldiev
    >>> ring = PolynomialRing(QQ, "z1, z2")
    >>> check(get_Beldiev(ring), max_iter=200)
    True
    """
    solver = LieBasisSolver(list(generators), max_iter)
    return solver.run()
