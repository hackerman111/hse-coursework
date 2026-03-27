"""
Shared public type aliases for the numeric Lie-generation layer.
"""

from __future__ import annotations

from typing import Dict, Tuple, TypeAlias


Monomial: TypeAlias = Tuple[int, ...]
Generator: TypeAlias = Dict[int, Dict[Monomial, float]]

