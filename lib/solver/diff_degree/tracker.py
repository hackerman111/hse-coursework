"""
Exact-базис и мономиальный трекер для символьной проверки.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Any, DefaultDict

from ...derivation import LieDerivation
from ..operations import reduce_derivation


MonomialKey = tuple[tuple[int, ...], int]


@dataclass(frozen=True)
class TrackedDerivation:
    id: int
    derivation: LieDerivation
    degree: int
    source: str


@dataclass(frozen=True)
class InsertResult:
    tracked: TrackedDerivation | None
    monomial_added: bool
    partial_added: bool


class SymbolicSubspaceTracker:
    """
    Exact-аналоги REF_INSERT и MONOMIAL_TRACK поверх ``LieDerivation``.
    """

    def __init__(self, algebra: Any, max_degree: int) -> None:
        self.algebra = algebra
        self.gens = algebra.gens()
        self.n = len(self.gens)
        self.max_degree = max_degree

        self.basis: dict[tuple[int, Any], LieDerivation] = {}
        self.elements: list[TrackedDerivation] = []
        self.by_degree: DefaultDict[int, list[TrackedDerivation]] = defaultdict(list)
        self.monomials: set[MonomialKey] = set()
        self.partial_axes: set[int] = set()

        self._next_id = 0

    @property
    def total_rank(self) -> int:
        return len(self.elements)

    def degree_rank(self, degree: int) -> int:
        return len(self.by_degree.get(degree, ()))

    def rank_upto(self, degree: int) -> int:
        return sum(len(items) for current, items in self.by_degree.items() if current <= degree)

    def elements_of_degree(self, degree: int) -> tuple[TrackedDerivation, ...]:
        return tuple(self.by_degree.get(degree, ()))

    def has_monomial(self, exponents: tuple[int, ...], component: int) -> bool:
        return (tuple(exponents), component) in self.monomials

    def insert(self, derivation: LieDerivation, *, source: str) -> InsertResult:
        truncated = self.truncate(derivation)
        if truncated is None:
            return InsertResult(None, False, False)

        monomial_added, partial_added = self.observe(truncated)

        reduced = reduce_derivation(truncated, self.basis)
        leading_term = reduced.leading_term
        if leading_term is None:
            return InsertResult(None, monomial_added, partial_added)

        component, monomial, coefficient = leading_term
        normalized = reduced * (1 / coefficient)
        tracked = TrackedDerivation(
            id=self._next_id,
            derivation=normalized,
            degree=normalized.degree(),
            source=source,
        )
        self._next_id += 1

        self.basis[(component, monomial)] = normalized
        self.elements.append(tracked)
        self.by_degree[tracked.degree].append(tracked)
        return InsertResult(tracked, monomial_added, partial_added)

    def observe(self, derivation: LieDerivation) -> tuple[bool, bool]:
        term = self.extract_monomial(derivation)
        if term is None:
            return False, False

        exponents, component, _ = term
        key = (exponents, component)
        monomial_added = key not in self.monomials
        if monomial_added:
            self.monomials.add(key)

        partial_added = False
        if sum(exponents) == 0 and component not in self.partial_axes:
            self.partial_axes.add(component)
            partial_added = True

        return monomial_added, partial_added

    def truncate(self, derivation: LieDerivation) -> LieDerivation | None:
        mapping = {}
        all_zero = True
        for gen in self.gens:
            poly = self._truncate_polynomial(derivation(gen))
            mapping[gen] = poly
            if poly != 0:
                all_zero = False
        if all_zero:
            return None
        return LieDerivation.from_mapping(self.algebra, mapping)

    def make_monomial_derivation(
        self,
        component: int,
        exponents: tuple[int, ...],
        *,
        coefficient: Any = 1,
    ) -> LieDerivation:
        mapping = {gen: 0 for gen in self.gens}
        mapping[self.gens[component]] = self.algebra(
            {self._ring_key(tuple(exponents)): coefficient}
        )
        return LieDerivation.from_mapping(self.algebra, mapping)

    def extract_monomial(
        self,
        derivation: LieDerivation,
    ) -> tuple[tuple[int, ...], int, Any] | None:
        found: tuple[tuple[int, ...], int, Any] | None = None
        for index, gen in enumerate(self.gens):
            poly = self.algebra(derivation(gen))
            if poly == 0:
                continue
            terms = poly.dict()
            if len(terms) != 1:
                return None
            if found is not None:
                return None
            (exponents, coefficient), = terms.items()
            found = (self._normalize_exponents(exponents), index, coefficient)
        return found

    def _truncate_polynomial(self, poly: Any):
        poly = self.algebra(poly)
        if poly == 0:
            return poly

        kept = {
            self._ring_key(self._normalize_exponents(exponents)): coefficient
            for exponents, coefficient in poly.dict().items()
            if sum(self._normalize_exponents(exponents)) <= self.max_degree
        }
        if not kept:
            return self.algebra.zero()
        return self.algebra(kept)

    def _normalize_exponents(self, exponents: Any) -> tuple[int, ...]:
        try:
            return tuple(int(value) for value in exponents)
        except TypeError:
            return (int(exponents),)

    def _ring_key(self, exponents: tuple[int, ...]):
        if self.n == 1:
            return exponents[0]
        return tuple(exponents)
