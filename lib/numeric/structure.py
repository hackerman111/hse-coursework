"""
Public wrapper around precomputed structure constants.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .bracket import precompute_structure_constants


StructureConstantArrays = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


@dataclass
class StructureConstantCache:
    """
    Precompute all structure constants needed for one truncated run.

    The cache is reusable across many checker instances with the same ``n`` and
    ``D_max``.
    """

    n: int
    D_max: int
    cache: dict[tuple[int, int], StructureConstantArrays] = field(
        default_factory=dict,
        init=False,
    )

    def __post_init__(self) -> None:
        for p in range(-1, self.D_max + 1):
            for q in range(p, self.D_max + 1):
                r = p + q
                if -1 <= r <= self.D_max:
                    self.cache[(p, q)] = precompute_structure_constants(self.n, p, q)

    def get(self, p: int, q: int) -> StructureConstantArrays | None:
        return self.cache.get((p, q))

