"""
Вспомогательные операции solver.
"""

from __future__ import annotations

from typing import Dict, Iterable, Tuple

from ..derivation import LieDerivation


def reduce_derivation(
    derivation: LieDerivation,
    basis: Dict[Tuple[int, object], LieDerivation],
) -> LieDerivation:
    while True:
        leading_term = derivation.leading_term
        if leading_term is None:
            return derivation

        comp_idx, monomial, coeff = leading_term
        basis_key = (comp_idx, monomial)
        if basis_key not in basis:
            return derivation

        derivation = derivation - (basis[basis_key] * coeff)


def generate_commutators(
    new_derivation: LieDerivation,
    processed: Iterable[LieDerivation],
) -> list[LieDerivation]:
    return [new_derivation.bracket(old_derivation) for old_derivation in processed]
