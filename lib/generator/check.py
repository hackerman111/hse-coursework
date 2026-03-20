"""
Проверка полноты набора генераторов.
"""

from __future__ import annotations

from typing import Iterable

from ..solver import DiffDegreeLieGenerationChecker, LieBasisSolver


def check(generators: Iterable, max_iter: int = 1000) -> bool:
    solver = LieBasisSolver(list(generators), max_iter)
    return solver.run()


def check_diff_degree(
    generators: Iterable,
    degree: int,
    *,
    work_degree: int | None = None,
    max_iter: int = 5000,
    logger=None,
    verbose: bool = False,
    log_every_s: float | None = None,
):
    solver = DiffDegreeLieGenerationChecker(
        list(generators),
        degree,
        work_degree=work_degree,
        max_iter=max_iter,
        logger=logger,
        verbose=verbose,
        log_every_s=log_every_s,
    )
    return solver.run()
