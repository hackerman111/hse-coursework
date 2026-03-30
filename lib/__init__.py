"""
Публичная точка входа библиотеки.

Не импортирует Sage-зависимые модули заранее, чтобы чисто численные
подпакеты (`lib.numeric`) можно было использовать без побочных эффектов.

Публичный API (без SageMath):
    from lib import monomial_spec, partial_spec, combine_specs, spec_degree, random_spec
    from lib import spec_from_string, spec_to_string
    from lib import beldiev_specs, andrist_specs
    from lib import check_generation
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from importlib import import_module
from typing import Dict, Tuple

# --- Публичный интерфейс: не требует SageMath -------------------------------

from lib.core.spec import (  # noqa: F401
    DerivationSpec,
    LeadingTerm,
    combine_specs,
    homogeneous_component,
    homogeneous_components,
    leading_term_spec,
    monomial_spec,
    normalize_spec,
    partial_spec,
    random_spec,
    scale_spec,
    spec_degree,
    zero_spec,
)
from lib.core.results import SymbolicResult  # noqa: F401
from lib.backends.numeric import (  # noqa: F401
    from_numeric_components,
    to_generator,
    to_numeric_components,
)
from lib.numeric import (  # noqa: F401
    DegreeStatus,
    FULL_LIE_GENERATION_DEGREE,
    FullLieGenerationCheck,
    HypothesisRunResult,
    HypothesisTrialResult,
    LieGenerationCriterionResult,
    NumericLieGenerationChecker,
    NumericResult,
    StructureConstantCache,
    check_full_lie_generation,
    check_numeric_generation,
    describe_generator,
    describe_random_generator,
    evaluate_full_lie_generation,
    generator_max_degree,
    generator_term_count,
    make_beldiev_generators,
    make_general_monomial_generator,
    make_monomial_D,
    make_monomial_generator,
    make_random_E,
    make_random_generator,
    parse_generator_spec,
    resolve_hypothesis_generator,
    run_beldiev,
    run_hypothesis,
    satisfies_full_lie_generation,
)
from lib.backends.spec_ops import ad, lie_bracket  # noqa: F401
from lib.io.parse import spec_from_string, spec_to_string  # noqa: F401
from lib.generators.beldiev import beldiev_specs  # noqa: F401
from lib.generators.andrist import andrist_specs  # noqa: F401
from lib.solver import check_generation  # noqa: F401

# --- Лениво загружаемые символьные зависимости (требуют SageMath) ----------

_LAZY_EXPORTS: Dict[str, Tuple[str, str]] = {
    "Derivation": ("lib.derivation", "Derivation"),
    "LieDerivation": ("lib.derivation", "LieDerivation"),
    "LieDerivationFactory": ("lib.derivation", "LieDerivationFactory"),
    "QuotientLieDerivation": ("lib.derivation", "QuotientLieDerivation"),
    "LibraryLogger": ("lib.utils", "LibraryLogger"),
    "LieBasisSolver": ("lib.solver", "LieBasisSolver"),
    "check": ("lib.generator", "check"),
    "check_symbolic_degree2_criterion": ("lib.symbolic", "check_symbolic_degree2_criterion"),
    "check_symbolic_generation": ("lib.symbolic", "check_symbolic_generation"),
    "enumerate_targets_up_to_degree2": ("lib.symbolic", "enumerate_targets_up_to_degree2"),
    "format_target": ("lib.symbolic", "format_target"),
    "get_Beldiev": ("lib.generator", "get_Beldiev"),
    "get_Andristy": ("lib.generator", "get_Andristy"),
    "make_degree2_stop_condition": ("lib.symbolic", "make_degree2_stop_condition"),
    "make_symbolic_algebra": ("lib.symbolic", "make_symbolic_algebra"),
    "PolynomialRing": ("sage.all", "PolynomialRing"),
    "spec_to_symbolic_derivation": ("lib.symbolic", "spec_to_symbolic_derivation"),
    "SymbolicCriterionResult": ("lib.symbolic", "SymbolicCriterionResult"),
    "example_util": ("lib.utils", "example_util"),
}

__all__ = [
    # Канонический spec-интерфейс
    "DerivationSpec",
    "LeadingTerm",
    "monomial_spec",
    "partial_spec",
    "combine_specs",
    "spec_degree",
    "random_spec",
    "zero_spec",
    "scale_spec",
    "normalize_spec",
    "homogeneous_component",
    "homogeneous_components",
    "leading_term_spec",
    # Численный интерфейс
    "DegreeStatus",
    "FULL_LIE_GENERATION_DEGREE",
    "FullLieGenerationCheck",
    "HypothesisRunResult",
    "HypothesisTrialResult",
    "LieGenerationCriterionResult",
    "NumericLieGenerationChecker",
    "NumericResult",
    "StructureConstantCache",
    "SymbolicResult",
    "to_generator",
    "to_numeric_components",
    "from_numeric_components",
    "check_full_lie_generation",
    "check_numeric_generation",
    "describe_generator",
    "describe_random_generator",
    "evaluate_full_lie_generation",
    "generator_max_degree",
    "generator_term_count",
    "lie_bracket",
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
    "ad",
    # Ввод и вывод
    "spec_from_string",
    "spec_to_string",
    # Готовые семейства генераторов
    "beldiev_specs",
    "andrist_specs",
    # Решатель
    "check_generation",
    # Символьная часть (ленивая загрузка, требует SageMath)
    "ABC",
    "Derivation",
    "LibraryLogger",
    "LieBasisSolver",
    "LieDerivation",
    "LieDerivationFactory",
    "PolynomialRing",
    "QuotientLieDerivation",
    "SymbolicCriterionResult",
    "abstractmethod",
    "check",
    "check_symbolic_degree2_criterion",
    "check_symbolic_generation",
    "enumerate_targets_up_to_degree2",
    "example_util",
    "format_target",
    "get_Andristy",
    "get_Beldiev",
    "make_degree2_stop_condition",
    "make_symbolic_algebra",
    "spec_to_symbolic_derivation",
]


def __getattr__(name: str):
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        value = getattr(import_module(module_name), attr_name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
