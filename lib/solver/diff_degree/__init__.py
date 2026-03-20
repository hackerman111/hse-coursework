"""
Публичный API diff-degree checker.
"""

from .checker import DiffDegreeLieGenerationChecker
from .dimensions import dim_Wn_coeff_k, dim_Wn_coeff_leq_d
from .result import DiffDegreeLayerStatus, DiffDegreeResult

__all__ = [
    "DiffDegreeLayerStatus",
    "DiffDegreeLieGenerationChecker",
    "DiffDegreeResult",
    "dim_Wn_coeff_k",
    "dim_Wn_coeff_leq_d",
]
