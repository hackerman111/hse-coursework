"""
Критерий порождения: все дифференцирования степени <= 2.

Перечисляет целевые мономы (degree -1, 0, 1, 2) для W_n и строит
stop_condition callback для LieBasisSolver.
"""

from __future__ import annotations

from itertools import combinations_with_replacement
from typing import Any, Callable, Dict, List, Tuple


def enumerate_targets_up_to_degree2(algebra) -> List[Tuple[int, Any]]:
    """
    Перечислить все целевые (comp_idx, monomial) пары для степеней -1, 0, 1, 2.

    Степень дифференцирования z^alpha * ∂/∂z_i равна |alpha| - 1.
    """
    gens = algebra.gens()
    n = len(gens)
    targets: List[Tuple[int, Any]] = []

    for axis in range(n):
        # degree -1: ∂/∂z_i (monomial = 1)
        targets.append((axis, algebra.one()))

        # degree 0: z_j · ∂/∂z_i for all j
        for j in range(n):
            targets.append((axis, gens[j]))

        # degree 1: z_j*z_k · ∂/∂z_i for j <= k
        for j, k in combinations_with_replacement(range(n), 2):
            targets.append((axis, gens[j] * gens[k]))

        # degree 2: z_j*z_k*z_l · ∂/∂z_i for j <= k <= l
        for j, k, l in combinations_with_replacement(range(n), 3):
            targets.append((axis, gens[j] * gens[k] * gens[l]))

    return targets


def make_degree2_stop_condition(
    algebra,
) -> Tuple[Callable, Dict[Tuple[int, Any], bool]]:
    """
    Построить stop_condition callback и словарь статуса целей.

    Returns:
        (stop_condition, targets_status):
            stop_condition — callable(solver) -> bool
            targets_status — dict {(axis, monomial): bool} для отслеживания
    """
    targets = enumerate_targets_up_to_degree2(algebra)
    targets_status: Dict[Tuple[int, Any], bool] = {t: False for t in targets}

    def stop_condition(solver) -> bool:
        for target_key in targets_status:
            if not targets_status[target_key] and target_key in solver.basis:
                targets_status[target_key] = True
        return all(targets_status.values())

    return stop_condition, targets_status


def format_target(axis: int, monomial, algebra) -> str:
    """Человекочитаемое описание целевого дифференцирования."""
    gens = algebra.gens()
    mono_str = str(monomial) if monomial != 1 else ""
    var_name = str(gens[axis])
    if mono_str:
        return f"{mono_str}·∂/∂{var_name}"
    return f"∂/∂{var_name}"
