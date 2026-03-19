"""
Общие утилиты для построения генераторов.
"""

from __future__ import annotations

from functools import reduce
import operator
from typing import Any

from ..derivation import LieDerivation


def require_dimension(algebra: Any, theorem_name: str, minimum: int = 2):
    gens = algebra.gens()
    dimension = len(gens)
    if dimension < minimum:
        raise ValueError(f"{theorem_name} требует n >= {minimum}")
    return gens, dimension


def zero_mapping(gens: tuple[Any, ...] | list[Any]) -> dict[Any, Any]:
    return {gen: 0 for gen in gens}


def partial_derivative(algebra: Any, gens: tuple[Any, ...], target_index: int) -> LieDerivation:
    mapping = zero_mapping(gens)
    mapping[gens[target_index]] = 1
    return LieDerivation.from_mapping(algebra, mapping)


def multiply_all(values: tuple[Any, ...] | list[Any]) -> Any:
    return reduce(operator.mul, values, 1)
