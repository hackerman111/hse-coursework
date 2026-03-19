"""
Фабрика для выбора корректной обертки над Sage-дифференцированием.
"""

from __future__ import annotations

from typing import Any

from .base import LieDerivation
from .quotient import QuotientLieDerivation


class LieDerivationFactory:
    @staticmethod
    def create(sage_derivation: Any) -> LieDerivation:
        algebra = sage_derivation.codomain()
        if hasattr(algebra, "cover"):
            return QuotientLieDerivation(sage_derivation)
        return LieDerivation(sage_derivation)
