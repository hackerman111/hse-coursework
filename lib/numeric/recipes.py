"""
Named generator families used by numeric experiments.
"""

from __future__ import annotations

from .types import Generator


def make_beldiev_generators(n: int) -> list[Generator]:
    """
    Beldiev's explicit 2-generator family for ``W_n``.
    """
    u: Generator = {(n - 1): {(0,) * n: 1.0}}
    v: Generator = {}

    for j in range(n - 1):
        target_var = n - 2 - j
        power = 4 * (n - j)
        alpha = [0] * n
        for index in range(n - 1, n - 2 - j, -1):
            alpha[index] = power
        v.setdefault(target_var, {})
        v[target_var][tuple(alpha)] = 1.0

    v.setdefault(n - 1, {})
    v[n - 1][tuple([4] * n)] = 1.0
    return [u, v]

