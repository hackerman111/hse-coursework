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

# --- Публичный API: не требует SageMath ------------------------------------

from lib.core.spec import (  # noqa: F401
    DerivationSpec,
    monomial_spec,
    partial_spec,
    combine_specs,
    spec_degree,
    random_spec,
    to_generator,
    lie_bracket,
    ad,
)
from lib.io.parse import spec_from_string, spec_to_string  # noqa: F401
from lib.generators.beldiev import beldiev_specs  # noqa: F401
from lib.generators.andrist import andrist_specs  # noqa: F401
from lib.solver import check_generation  # noqa: F401

# --- Lazy-загружаемые символьные зависимости (требуют SageMath) -----------

_LAZY_EXPORTS: Dict[str, Tuple[str, str]] = {
    "Derivation": ("lib.derivation", "Derivation"),
    "LieDerivation": ("lib.derivation", "LieDerivation"),
    "LieDerivationFactory": ("lib.derivation", "LieDerivationFactory"),
    "QuotientLieDerivation": ("lib.derivation", "QuotientLieDerivation"),
    "LibraryLogger": ("lib.utils", "LibraryLogger"),
    "LieBasisSolver": ("lib.solver", "LieBasisSolver"),
    "check": ("lib.generator", "check"),
    "get_Beldiev": ("lib.generator", "get_Beldiev"),
    "get_Andristy": ("lib.generator", "get_Andristy"),
    "PolynomialRing": ("sage.all", "PolynomialRing"),
    "example_util": ("lib.utils", "example_util"),
}

__all__ = [
    # Canonical spec API
    "DerivationSpec",
    "monomial_spec",
    "partial_spec",
    "combine_specs",
    "spec_degree",
    "random_spec",
    "to_generator",
    "lie_bracket",
    "ad",
    # I/O
    "spec_from_string",
    "spec_to_string",
    # Generator recipes
    "beldiev_specs",
    "andrist_specs",
    # Solver
    "check_generation",
    # Symbolic (lazy, require SageMath)
    "ABC",
    "Derivation",
    "LibraryLogger",
    "LieBasisSolver",
    "LieDerivation",
    "LieDerivationFactory",
    "PolynomialRing",
    "QuotientLieDerivation",
    "abstractmethod",
    "check",
    "example_util",
    "get_Andristy",
    "get_Beldiev",
]


def __getattr__(name: str):
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        value = getattr(import_module(module_name), attr_name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
