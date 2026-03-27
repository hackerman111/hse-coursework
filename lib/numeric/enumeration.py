"""
Помощники для конечного перебора кандидатов-генераторов (n=2, backward compat).

Универсальная версия для произвольного n находится в lib/solver/enumeration.py.
Общие утилиты (parse_coefficients, estimate_candidate_count, format_coefficient,
iter_candidates, candidate_to_spec) реэкспортируются оттуда.

Здесь сохраняются n=2-специфичные BasisTerm(ax, ay) и build_basis_terms(max_x_deg, max_y_deg)
для обратной совместимости с numeric_enumerate.py и тестами.
"""

from __future__ import annotations

from dataclasses import dataclass

from .types import Generator

# --- Реэкспорт общих утилит из lib/solver/enumeration.py ----------------

from lib.solver.enumeration import (  # noqa: F401
    estimate_candidate_count,
    format_coefficient,
    iter_candidates,
    parse_coefficients,
)


# --- n=2-специфичные компоненты (backward compat) -----------------------


@dataclass(frozen=True)
class BasisTerm:
    """
    Базисный терм для n=2 (backward compat).

    Используется в numeric_enumerate.py и связанных тестах.
    Для произвольного n используйте lib.solver.enumeration.BasisTerm.
    """

    target_var: int
    ax: int
    ay: int

    def factor_spec(self) -> str:
        factors: list[str] = []
        if self.ax == 1:
            factors.append("x")
        elif self.ax > 1:
            factors.append(f"x^{self.ax}")

        if self.ay == 1:
            factors.append("y")
        elif self.ay > 1:
            factors.append(f"y^{self.ay}")

        factors.append("dx" if self.target_var == 0 else "dy")
        return "*".join(factors)


def build_basis_terms(
    max_x_deg: int,
    max_y_deg: int,
    include_constants: bool,
) -> list[BasisTerm]:
    """Строит список базисных термов для n=2 (backward compat)."""
    terms: list[BasisTerm] = []
    for target_var in (0, 1):
        for ax in range(max_x_deg + 1):
            for ay in range(max_y_deg + 1):
                if not include_constants and ax == 0 and ay == 0:
                    continue
                terms.append(BasisTerm(target_var=target_var, ax=ax, ay=ay))
    return terms


def build_generator(candidate: list[tuple[BasisTerm, float]]) -> Generator:
    """Строит Generator dict из списка (BasisTerm, коэффициент) для n=2."""
    generator: Generator = {}
    for term, coeff in candidate:
        generator.setdefault(term.target_var, {})
        alpha = (term.ax, term.ay)
        generator[term.target_var][alpha] = generator[term.target_var].get(alpha, 0.0) + coeff
    return {target: monomials for target, monomials in generator.items() if monomials}


def candidate_to_spec(candidate: list[tuple[BasisTerm, float]]) -> str:
    """Форматирует кандидата в строку (n=2, backward compat)."""
    parts = []
    for term, coeff in candidate:
        factor = term.factor_spec()
        if abs(coeff - 1.0) <= 1e-12:
            parts.append(factor)
        else:
            parts.append(f"{format_coefficient(coeff)}*{factor}")
    return " + ".join(parts)
