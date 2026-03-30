"""
Conversions between sparse generator specs and graded numeric components.
"""

from __future__ import annotations

import numpy as np

from .indexing import dim_Wn_k, multiindex_lookup
from .types import Generator


def decompose_generator(
    generator: Generator,
    n: int,
    D_max: int,
) -> dict[int, np.ndarray]:
    """
    Разложить разреженный генератор по однородным компонентам стандартной градуировки.

    На вход принимает generator формата ``{axis: {alpha: coeff}}``, размерность ``n`` и порог усечения ``D_max``.
    На выходе возвращает словарь ``degree -> np.ndarray`` с координатами проекций в ``W_n^(degree)``.

    >>> parts = decompose_generator({0: {(0, 0): 1.0}, 1: {(2, 0): 2.0}}, 2, 2)
    >>> sorted(parts)
    [-1, 1]
    """
    components: dict[int, np.ndarray] = {}

    for axis, monomials in generator.items():
        for alpha, coeff in monomials.items():
            degree = sum(alpha) - 1
            if degree < -1 or degree > D_max:
                continue
            if degree not in components:
                components[degree] = np.zeros(dim_Wn_k(n, degree), dtype=np.float64)
            lookup = multiindex_lookup(n, sum(alpha))
            index = axis * len(lookup) + lookup[alpha]
            components[degree][index] += coeff

    return components
