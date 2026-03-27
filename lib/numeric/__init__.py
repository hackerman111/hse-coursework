"""
Численный модуль для проверки порождения W_n до степени d.

Не зависит от SageMath — работает на чистом numpy.
"""

from .criteria import (
    FULL_LIE_GENERATION_DEGREE,
    LieGenerationCriterionResult,
    evaluate_full_lie_generation,
    satisfies_full_lie_generation,
)
from .recipes import make_beldiev_generators
from .results import DegreeStatus, NumericResult
from .solver import NumericLieGenerationChecker
from .specs import (
    describe_generator,
    describe_random_generator,
    generator_max_degree,
    generator_term_count,
    make_general_monomial_generator,
    make_monomial_D,
    make_monomial_generator,
    make_random_E,
    make_random_generator,
    parse_generator_spec,
)
from .structure import StructureConstantCache
from .workflows import (
    FullLieGenerationCheck,
    HypothesisRunResult,
    HypothesisTrialResult,
    check_full_lie_generation,
    check_numeric_generation,
    resolve_hypothesis_generator,
    run_beldiev,
    run_hypothesis,
)

__all__ = [
    "DegreeStatus",
    "FULL_LIE_GENERATION_DEGREE",
    "FullLieGenerationCheck",
    "HypothesisRunResult",
    "HypothesisTrialResult",
    "LieGenerationCriterionResult",
    "NumericLieGenerationChecker",
    "NumericResult",
    "StructureConstantCache",
    "check_full_lie_generation",
    "check_numeric_generation",
    "describe_generator",
    "describe_random_generator",
    "evaluate_full_lie_generation",
    "generator_max_degree",
    "generator_term_count",
    "make_beldiev_generators",
    "make_general_monomial_generator",
    "make_monomial_D",
    "make_monomial_generator",
    "make_random_E",
    "make_random_generator",
    "parse_generator_spec",
    "resolve_hypothesis_generator",
    "run_beldiev",
    "run_hypothesis",
    "satisfies_full_lie_generation",
]
