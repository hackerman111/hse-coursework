"""
Публичная точка входа библиотеки.

Не импортирует Sage-зависимые модули заранее, чтобы чисто численные
подпакеты (`lib.numeric`) можно было использовать без побочных эффектов.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from importlib import import_module
from typing import Dict, Tuple


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
