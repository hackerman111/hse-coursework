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
    Перечислить все целевые мономиальные направления до второй степени включительно.

    На вход принимает полиномиальную алгебру Sage ``algebra``.
    На выходе возвращает список пар ``(axis, monomial)``, соответствующих всем базисным элементам степеней ``-1, 0, 1, 2``.

    >>> from sage.all import PolynomialRing, QQ
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> len(enumerate_targets_up_to_degree2(ring))
    20
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
    Построить критерий остановки solver-а по покрытию всех целей до степени 2.

    На вход принимает полиномиальную алгебру Sage ``algebra``.
    На выходе возвращает пару ``(stop_condition, targets_status)`` для отслеживания найденных базисных элементов.

    >>> from sage.all import PolynomialRing, QQ
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> stop_condition, targets_status = make_degree2_stop_condition(ring)
    >>> callable(stop_condition) and len(targets_status) == 20
    True
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
    """
    Отформатировать целевое дифференцирование в удобочитаемую строку.

    На вход принимает индекс компоненты ``axis``, моном ``monomial`` и алгебру Sage ``algebra``.
    На выходе возвращает строку вида ``x*y·∂/∂x``.

    >>> from sage.all import PolynomialRing, QQ
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> format_target(0, ring.gens()[1], ring)
    'y·∂/∂x'
    """
    gens = algebra.gens()
    mono_str = str(monomial) if monomial != 1 else ""
    var_name = str(gens[axis])
    if mono_str:
        return f"{mono_str}·∂/∂{var_name}"
    return f"∂/∂{var_name}"
