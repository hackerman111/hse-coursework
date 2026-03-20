"""
Результаты символьной проверки порождения W_n.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DiffDegreeLayerStatus:
    degree: int
    rank: int
    dimension: int

    @property
    def is_full(self) -> bool:
        return self.rank >= self.dimension

    def __str__(self) -> str:
        mark = "full" if self.is_full else "open"
        return f"deg {self.degree:2d}: {self.rank}/{self.dimension} ({mark})"


@dataclass(frozen=True)
class DiffDegreeResult:
    success: bool
    phase: str
    sufficient_condition: bool
    iterations: int
    elapsed_s: float
    target_rank: int
    target_dim: int
    total_rank: int
    work_degree: int
    partial_axes: tuple[int, ...]
    monomial_count: int
    degree_one_rank: int
    layers: tuple[DiffDegreeLayerStatus, ...]

    def summary(self) -> str:
        headline = "PASS" if self.success else "FAIL"
        reason = "H2" if self.sufficient_condition else "rank"
        lines = [
            "=== Diff-Degree Symbolic Check ===",
            f"status = {headline}, phase = {self.phase}, reason = {reason}",
            f"elapsed = {self.elapsed_s:.3f}s, iterations = {self.iterations}",
            (
                f"target rank <= {self.layers[-1].degree if self.layers else 0} = "
                f"{self.target_rank}/{self.target_dim}"
            ),
            f"total rank <= {self.work_degree} = {self.total_rank}",
            f"partials = {self.partial_axes}, degree-1 rank = {self.degree_one_rank}",
            f"tracked monomials = {self.monomial_count}",
        ]
        for layer in self.layers:
            lines.append(str(layer))
        return "\n".join(lines)
