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


class QuotientLieDerivation(LieDerivation):
    """
    Специализированная версия дифференцирования для работы в фактор-кольцах (R/I).
    """

    def __init__(self, sage_derivation: Any) -> None:
        super().__init__(sage_derivation)
        if not hasattr(self._algebra, "cover"):
            raise ValueError(
                "QuotientLieDerivation требует фактор-алгебру (с методом cover())"
            )

        if hasattr(self._algebra, "cover_ring"):
            self._base_ring = self._algebra.cover_ring()
        else:
            cover = self._algebra.cover()
            self._base_ring = cover.domain() if hasattr(cover, "domain") else cover

    def _extract_poly(self, element: Any) -> Any:
        """
        Подъём элемента из факторкольца в накрывающее кольцо.

        Позволяет унаследованным leading_term и degree работать с lm()/degree()
        исходного кольца многочленов, а не факторкольца.
        """
        if hasattr(element, "lift"):
            return element.lift()
        return element

    def __repr__(self) -> str:
        return f"QuotientLieD({self._d})"
