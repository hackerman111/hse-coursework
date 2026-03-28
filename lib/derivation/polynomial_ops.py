"""
Helpers for inspecting polynomial coefficients of symbolic derivations.
"""

from __future__ import annotations

from typing import Any, Protocol


class PolynomialOps(Protocol):
    def representative(self, element: Any) -> Any:
        """Return a polynomial representative suitable for inspection."""

    def leading_term(self, element: Any) -> tuple[Any, Any] | None:
        """Return `(lm, lc)` or `None` for zero."""

    def degree(self, element: Any) -> int:
        """Return the polynomial degree or `-1` for zero."""


class BasePolynomialOps:
    """Inspection helpers for ordinary polynomial rings."""

    def representative(self, element: Any) -> Any:
        return element

    def leading_term(self, element: Any) -> tuple[Any, Any] | None:
        poly = self.representative(element)
        if poly == 0:
            return None

        try:
            return poly.lm(), poly.lc()
        except AttributeError:
            return 1, poly

    def degree(self, element: Any) -> int:
        poly = self.representative(element)
        if poly == 0:
            return -1

        try:
            return poly.degree()
        except AttributeError:
            return 0


class QuotientPolynomialOps(BasePolynomialOps):
    """Inspection helpers for quotient rings via lifted representatives."""

    def representative(self, element: Any) -> Any:
        if hasattr(element, "lift"):
            return element.lift()
        return element
