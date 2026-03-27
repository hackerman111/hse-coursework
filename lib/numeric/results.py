"""
Result objects for numeric Lie-generation checks.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class DegreeStatus:
    """Status of one homogeneous degree layer."""

    degree: int
    rank: int
    target_dim: int

    @property
    def is_full(self) -> bool:
        return self.rank >= self.target_dim

    def __repr__(self) -> str:
        tag = "FULL" if self.is_full else f"{self.rank}/{self.target_dim}"
        return f"Degree {self.degree:3d}: {tag}"


@dataclass
class NumericResult:
    """Result of a truncated numeric Lie-generation run."""

    n: int
    d: int
    D_max: int
    success: bool
    degrees: list[DegreeStatus] = field(default_factory=list)
    iterations: int = 0
    elapsed_s: float = 0.0

    def degree_status(self, degree: int) -> DegreeStatus | None:
        """Return the stored status for one degree if it was checked."""
        for status in self.degrees:
            if status.degree == degree:
                return status
        return None

    def covers_degree(self, degree: int) -> bool:
        """Whether the run checked all layers up to the requested degree."""
        return degree <= self.d and self.degree_status(degree) is not None

    def missing_degrees(self, max_degree: int | None = None) -> list[DegreeStatus]:
        """
        Return incomplete target degrees.

        When ``max_degree`` is given, only degrees up to that bound are reported.
        """
        if max_degree is None:
            return [status for status in self.degrees if not status.is_full]
        return [
            status
            for status in self.degrees
            if status.degree <= max_degree and not status.is_full
        ]

    def is_full_up_to(self, degree: int) -> bool:
        """Whether every checked target degree up to ``degree`` is full."""
        if not self.covers_degree(degree):
            return False
        return not self.missing_degrees(max_degree=degree)

    def summary(self) -> str:
        lines = [
            "=== Numeric Lie Generation Check ===",
            f"n = {self.n}, d = {self.d}, D_max = {self.D_max}",
            f"Итераций: {self.iterations}, Время: {self.elapsed_s:.3f}s",
            (
                "Результат: PASS — все степени ≤ d полны"
                if self.success
                else "Результат: FAIL — есть неполные степени"
            ),
            "",
        ]
        for status in self.degrees:
            lines.append(str(status))
        return "\n".join(lines)

