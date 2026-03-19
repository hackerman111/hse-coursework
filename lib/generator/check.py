"""
Проверка полноты набора генераторов.
"""

from __future__ import annotations

from typing import Iterable

from ..solver import LieBasisSolver


def check(generators: Iterable, max_iter: int = 1000) -> bool:
    solver = LieBasisSolver(list(generators), max_iter)
    return solver.run()
