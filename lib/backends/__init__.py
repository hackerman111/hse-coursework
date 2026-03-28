"""
Backend-адаптеры: конвертация DerivationSpec в конкретные вычислительные представления.

Доступные адаптеры:
- lib.backends.numeric  — NumPy-представление (без зависимости от Sage)
- lib.backends.sage     — SageMath-представление (требует Sage)
"""

from .numeric import from_numeric_components, to_generator, to_numeric_components
from .spec_ops import ad, lie_bracket

__all__ = [
    "ad",
    "from_numeric_components",
    "lie_bracket",
    "to_generator",
    "to_numeric_components",
]
