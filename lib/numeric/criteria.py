"""
Project-level criteria built on top of truncated numeric generation checks.
"""

from __future__ import annotations

from dataclasses import dataclass

from .results import NumericResult


FULL_LIE_GENERATION_DEGREE = 2


@dataclass(frozen=True)
class LieGenerationCriterionResult:
    """Decision for a project-level full-generation criterion."""

    name: str
    required_degree: int
    satisfied: bool
    missing_degrees: tuple[int, ...]
    reason: str

    def summary(self) -> str:
        status = "PASS" if self.satisfied else "FAIL"
        return (
            f"{status}: {self.name} "
            f"(required_degree={self.required_degree}, missing={list(self.missing_degrees)})"
        )


def evaluate_full_lie_generation(
    result: NumericResult,
    required_degree: int = FULL_LIE_GENERATION_DEGREE,
) -> LieGenerationCriterionResult:
    """
    Project criterion: if all derivations up to degree 2 are generated, treat the
    set as generating the whole Lie algebra.
    """
    if not result.covers_degree(required_degree):
        return LieGenerationCriterionResult(
            name="full-lie-from-degree-bound",
            required_degree=required_degree,
            satisfied=False,
            missing_degrees=(),
            reason=(
                f"numeric result only checks up to degree {result.d}, "
                f"but degree {required_degree} is required"
            ),
        )

    missing = tuple(status.degree for status in result.missing_degrees(max_degree=required_degree))
    if missing:
        return LieGenerationCriterionResult(
            name="full-lie-from-degree-bound",
            required_degree=required_degree,
            satisfied=False,
            missing_degrees=missing,
            reason=(
                "criterion failed: not all derivations up to the required degree "
                "were generated"
            ),
        )

    return LieGenerationCriterionResult(
        name="full-lie-from-degree-bound",
        required_degree=required_degree,
        satisfied=True,
        missing_degrees=(),
        reason=(
            "criterion satisfied: all derivations up to the required degree "
            "were generated"
        ),
    )


def satisfies_full_lie_generation(
    result: NumericResult,
    required_degree: int = FULL_LIE_GENERATION_DEGREE,
) -> bool:
    """Convenience boolean wrapper around ``evaluate_full_lie_generation``."""
    return evaluate_full_lie_generation(
        result,
        required_degree=required_degree,
    ).satisfied

