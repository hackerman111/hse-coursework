"""
Генераторы по теореме Бельдиева.
"""

from __future__ import annotations

from typing import Any

from ..derivation import LieDerivation
from .shared import multiply_all, partial_derivative, require_dimension, zero_mapping


def get_Beldiev(algebra: Any):
    """
    Возвращает генераторы (U, V) для алгебры Ли векторных полей на K^n
    согласно статье Бельдиева.
    """
    gens, dimension = require_dimension(algebra, "Теорема Бельдиева")

    u_generator = partial_derivative(algebra, gens, dimension - 1)
    v_mapping = zero_mapping(gens)

    for index in range(dimension - 1):
        tail_vars = gens[index + 1 :]
        power = 4 * (index + 2)
        v_mapping[gens[index]] = multiply_all(tail_vars) ** power

    v_mapping[gens[dimension - 1]] = multiply_all(gens) ** 4
    v_generator = LieDerivation.from_mapping(algebra, v_mapping)
    return [u_generator, v_generator]
