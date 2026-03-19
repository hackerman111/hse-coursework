"""
Базовая обертка над дифференцированием SageMath.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple, Union

from .constructors import (
    from_jacobian as create_from_jacobian,
    from_linear as create_from_linear,
    from_mapping as create_from_mapping,
    from_weitzenbock as create_from_weitzenbock,
)


class LieDerivation:
    """
    Обертка над дифференцированием SageMath для удобства работы.
    """

    def __init__(self, sage_derivation: Any) -> None:
        self._d = sage_derivation
        self._algebra = sage_derivation.codomain()

    def __call__(self, arg: Any) -> Any:
        return self._d(arg)

    def __add__(self, other: Union["LieDerivation", Any]) -> "LieDerivation":
        from .factory import LieDerivationFactory

        if isinstance(other, LieDerivation):
            return LieDerivationFactory.create(self._d + other._d)
        return LieDerivationFactory.create(self._d + other)

    def __sub__(self, other: Union["LieDerivation", Any]) -> "LieDerivation":
        from .factory import LieDerivationFactory

        if isinstance(other, LieDerivation):
            return LieDerivationFactory.create(self._d - other._d)
        return LieDerivationFactory.create(self._d - other)

    def __mul__(self, other: Union["LieDerivation", Any]) -> "LieDerivation":
        from .factory import LieDerivationFactory

        if isinstance(other, LieDerivation):
            raise TypeError(
                "Нельзя умножать два дифференцирования напрямую. "
                "Используйте скобку Ли (bracket)."
            )
        return LieDerivationFactory.create(self._d * other)

    def __rmul__(self, other: Any) -> "LieDerivation":
        from .factory import LieDerivationFactory

        return LieDerivationFactory.create(other * self._d)

    def codomain(self) -> Any:
        return self._algebra

    def bracket(self, other: "LieDerivation") -> "LieDerivation":
        """
        [D1, D2] = D1(D2) - D2(D1)
        """
        from .factory import LieDerivationFactory

        if not isinstance(other, LieDerivation):
            raise ValueError("Можно вычислять скобку только с другим LieDerivation")

        gens = self._algebra.gens()
        images = [self(other(gen)) - other(self(gen)) for gen in gens]
        return LieDerivationFactory.create(self._algebra.derivation(images))

    @property
    def leading_term(self) -> Optional[Tuple[int, Any, Any]]:
        """
        Возвращает старший член дифференцирования (term).
        """
        best_term = None

        for index, gen in enumerate(self._algebra.gens()):
            poly = self(gen)
            if poly == 0:
                continue

            try:
                lm = poly.lm()
                lc = poly.lc()
            except AttributeError:
                lm = self._algebra(1)
                lc = poly

            if best_term is None:
                best_term = (index, lm, lc)
                continue

            current_lm = best_term[1]
            if lm > current_lm or (lm == current_lm and index > best_term[0]):
                best_term = (index, lm, lc)

        return best_term

    def is_triangular(self) -> bool:
        gens = self._algebra.gens()

        for index, gen in enumerate(gens):
            allowed = set(gens[index + 1 :])
            value = self(gen)

            if hasattr(value, "is_constant") and value.is_constant():
                continue
            if isinstance(value, (int, float)):
                continue

            try:
                actual_vars = set(value.variables())
            except AttributeError:
                continue

            if not actual_vars.issubset(allowed):
                return False

        return True

    def degree(self) -> int:
        """
        Возвращает максимальную степень представителей коэффициентов.
        """
        max_deg = -1
        for gen in self._algebra.gens():
            poly = self(gen)
            if poly == 0:
                continue

            try:
                degree = poly.degree()
            except AttributeError:
                degree = 0

            if degree > max_deg:
                max_deg = degree

        return max_deg

    @staticmethod
    def from_mapping(algebra: Any, gen_mapping: Any = None) -> "LieDerivation":
        return create_from_mapping(algebra, gen_mapping)

    @staticmethod
    def from_linear(algebra: Any, matrix: Any) -> "LieDerivation":
        return create_from_linear(algebra, matrix)

    @staticmethod
    def from_weitzenbock(algebra: Any, matrix: Any) -> "LieDerivation":
        return create_from_weitzenbock(algebra, matrix)

    @staticmethod
    def from_jacobian(algebra: Any, polynomials: Any) -> "LieDerivation":
        return create_from_jacobian(algebra, polynomials)

    def __repr__(self) -> str:
        return f"LieD({self._d})"
