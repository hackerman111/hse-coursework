"""
Численный модуль для проверки порождения W_n до степени d.

Не зависит от SageMath — работает на чистом numpy.
"""

from .solver import NumericLieGenerationChecker, NumericResult

__all__ = ["NumericLieGenerationChecker", "NumericResult"]
