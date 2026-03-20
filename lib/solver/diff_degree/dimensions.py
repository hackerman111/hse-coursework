"""
Размерности для градуировки по степени коэффициентов.
"""

from __future__ import annotations

from math import comb


def dim_Wn_coeff_k(n: int, degree: int) -> int:
    """
    dim пространства полиномиальных векторных полей с коэффициентами
    ровно степени ``degree``.
    """
    if degree < 0:
        return 0
    return n * comb(n + degree - 1, degree)


def dim_Wn_coeff_leq_d(n: int, degree: int) -> int:
    """
    dim пространства полиномиальных векторных полей с коэффициентами
    степени не выше ``degree``.
    """
    if degree < 0:
        return 0
    return n * comb(n + degree, n)
