"""
Библиотека дифференцирования алгебр
"""

from .derivation import (
    Derivation,
    LieDerivation,
    LieDerivationFactory,
    QuotientLieDerivation,
)
from .solver import LieBasisSolver
from .generator import check, get_Beldiev, get_Andristy
from .utils import *
from sage.all import PolynomialRing
from abc import ABC, abstractmethod

__all__ = [
    "ABC",
    "Derivation",
    "LieBasisSolver",
    "LieDerivation",
    "LieDerivationFactory",
    "PolynomialRing",
    "QuotientLieDerivation",
    "abstractmethod",
    "check",
    "get_Andristy",
    "get_Beldiev",
]
