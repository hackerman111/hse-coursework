"""
Поддержка дифференцирований в фактор-кольцах.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple

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

    def _lift_poly(self, element: Any) -> Any:
        if hasattr(element, "lift"):
            return element.lift()
        return element

    @property
    def leading_term(self) -> Optional[Tuple[int, Any, Any]]:
        best_term = None

        for index, gen in enumerate(self._algebra.gens()):
            value = self(gen)
            if value == 0:
                continue

            lifted = self._lift_poly(value)
            try:
                lm = lifted.lm()
                lc = self._algebra(lifted.lc())
            except AttributeError:
                lm = self._base_ring(1)
                lc = value

            if best_term is None:
                best_term = (index, lm, lc)
                continue

            current_lm = best_term[1]
            if lm > current_lm or (lm == current_lm and index > best_term[0]):
                best_term = (index, lm, lc)

        return best_term

    def degree(self) -> int:
        max_deg = -1

        for gen in self._algebra.gens():
            value = self(gen)
            if value == 0:
                continue

            lifted = self._lift_poly(value)
            try:
                degree = lifted.degree()
            except AttributeError:
                degree = 0

            if degree > max_deg:
                max_deg = degree

        return max_deg

    def __repr__(self) -> str:
        return f"QuotientLieD({self._d})"
