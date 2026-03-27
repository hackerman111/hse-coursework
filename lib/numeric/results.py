"""
Re-export shim: каноническое определение перенесено в lib/core/results.py.

Оставлен для обратной совместимости — все существующие импорты продолжают работать.
"""

from lib.core.results import DegreeStatus, NumericResult, SymbolicResult

__all__ = ["DegreeStatus", "NumericResult", "SymbolicResult"]
