"""
Поддержка дифференцирований в фактор-кольцах.

QuotientLieDerivation отличается от LieDerivation только тем, как извлекается
многочлен из элемента факторкольца: нужно вызвать .lift() перед lm()/degree().
Это выражается переопределением одного метода _extract_poly(), а не дублированием
leading_term и degree.
"""

from __future__ import annotations

from typing import Any

from .base import LieDerivation
from .polynomial_ops import QuotientPolynomialOps


class QuotientLieDerivation(LieDerivation):
    """
    Специализированная версия дифференцирования для работы в фактор-кольцах (R/I).
    """

    def __init__(self, sage_derivation: Any) -> None:
        super().__init__(sage_derivation, polynomial_ops=QuotientPolynomialOps())
        if not hasattr(self._algebra, "cover"):
            raise ValueError(
                "QuotientLieDerivation требует фактор-алгебру (с методом cover())"
            )

    def __repr__(self) -> str:
        return f"QuotientLieD({self._d})"
