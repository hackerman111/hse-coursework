"""
Совместимость со старым API.
"""

from __future__ import annotations

from typing import Any

from .base import LieDerivation


def Derivation(algebra: Any, gen_mapping: Any = None) -> LieDerivation:
    """
    Старый конструктор поверх LieDerivation.from_mapping.
    """
    return LieDerivation.from_mapping(algebra, gen_mapping)
