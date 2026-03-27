"""
Reusable helpers for finite numeric generator enumeration.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from .types import Generator


@dataclass(frozen=True)
class BasisTerm:
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
    terms: list[BasisTerm] = []
    for target_var in (0, 1):
        for ax in range(max_x_deg + 1):
            for ay in range(max_y_deg + 1):
                if not include_constants and ax == 0 and ay == 0:
                    continue
                terms.append(BasisTerm(target_var=target_var, ax=ax, ay=ay))
    return terms


def parse_coefficients(text: str) -> list[float]:
    coeffs: list[float] = []
    for raw in text.split(","):
        item = raw.strip()
        if not item:
            continue
        value = float(item)
        if abs(value) <= 1e-15:
            continue
        if not any(abs(value - seen) <= 1e-15 for seen in coeffs):
            coeffs.append(value)
    if not coeffs:
        raise ValueError("need at least one nonzero coefficient in --coeffs")
    return coeffs


def estimate_candidate_count(
    basis_size: int,
    coeff_count: int,
    min_terms: int,
    max_terms: int,
) -> int:
    total = 0
    for term_count in range(min_terms, max_terms + 1):
        total += math.comb(basis_size, term_count) * (coeff_count**term_count)
    return total


def build_generator(candidate: list[tuple[BasisTerm, float]]) -> Generator:
    generator: Generator = {}
    for term, coeff in candidate:
        generator.setdefault(term.target_var, {})
        alpha = (term.ax, term.ay)
        generator[term.target_var][alpha] = generator[term.target_var].get(alpha, 0.0) + coeff
    return {target: monomials for target, monomials in generator.items() if monomials}


def format_coefficient(value: float) -> str:
    if abs(value - round(value)) <= 1e-12:
        return str(int(round(value)))
    return f"{value:g}"


def candidate_to_spec(candidate: list[tuple[BasisTerm, float]]) -> str:
    parts = []
    for term, coeff in candidate:
        factor = term.factor_spec()
        if abs(coeff - 1.0) <= 1e-12:
            parts.append(factor)
        else:
            parts.append(f"{format_coefficient(coeff)}*{factor}")
    return " + ".join(parts)


def iter_candidates(
    basis_terms: list[BasisTerm],
    coeffs: list[float],
    min_terms: int,
    max_terms: int,
):
    for term_count in range(min_terms, max_terms + 1):
        for term_combo in itertools.combinations(basis_terms, term_count):
            for coeff_combo in itertools.product(coeffs, repeat=term_count):
                yield list(zip(term_combo, coeff_combo))

