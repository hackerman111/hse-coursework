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
    Split a sparse generator into homogeneous numeric components.

    ``generator`` has the form ``{axis: {alpha: coeff}}`` and the result maps
    each degree ``k`` to a coordinate vector in ``W_n^(k)``.

    WARNING: The returned components are degree-projections of a single element,
    NOT independent elements of a Lie subalgebra.  Do NOT bracket components
    from different degrees as if they were separate subalgebra members — this
    produces phantom brackets.  Use ``bracket_full_elements`` from ``bracket.py``
    to correctly compute Lie brackets of non-homogeneous elements.
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

