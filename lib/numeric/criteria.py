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
    Оценить проектный критерий полной порождённости по усечённому numeric-результату.

    На вход принимает ``NumericResult`` и требуемую степень ``required_degree``.
    На выходе возвращает структурированное решение ``LieGenerationCriterionResult``.

    >>> result = NumericResult(n=2, d=2, D_max=2, success=True, degrees=[], iterations=0, elapsed_s=0.0)
    >>> evaluate_full_lie_generation(result).satisfied
    False
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
    """
    Проверить проектный критерий полной порождённости в булевом виде.

    На вход принимает ``NumericResult`` и необязательную требуемую степень.
    На выходе возвращает ``True`` тогда и только тогда, когда критерий выполнен.

    >>> result = NumericResult(n=2, d=1, D_max=1, success=True, degrees=[], iterations=0, elapsed_s=0.0)
    >>> satisfies_full_lie_generation(result)
    False
    """
    return evaluate_full_lie_generation(
        result,
        required_degree=required_degree,
    ).satisfied
