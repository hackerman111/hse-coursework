"""
Генераторы по теореме Андриста.
"""

from __future__ import annotations

from typing import Any

from ..derivation import LieDerivation
from .shared import multiply_all, partial_derivative, require_dimension, zero_mapping


def get_Andristy(algebra: Any):
    """
    Возвращает генераторы (U, V, W) согласно статье Андриста.
    """
    gens, dimension = require_dimension(algebra, "Теорема Андриста")

    u_generator = partial_derivative(algebra, gens, dimension - 1)

    v_mapping = zero_mapping(gens)
    v_mapping[gens[dimension - 1]] = 1
    for index in range(dimension - 2, -1, -1):
        next_var = gens[index + 1]
        tail_vars = gens[index + 2 :]
        v_mapping[gens[index]] = multiply_all(tail_vars) * (next_var ** 3)
    v_generator = LieDerivation.from_mapping(algebra, v_mapping)

    w_mapping = zero_mapping(gens)
    w_mapping[gens[dimension - 1]] = (multiply_all(gens[: dimension - 1]) ** 2) * gens[
        dimension - 1
    ]
    w_generator = LieDerivation.from_mapping(algebra, w_mapping)

    return [u_generator, v_generator, w_generator]
