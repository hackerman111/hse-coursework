"""
Готовые математические рецепты: генераторы теорем Бельдиева и Андриста.
"""

from .beldiev import beldiev_specs
from .andrist import andrist_specs

__all__ = ["beldiev_specs", "andrist_specs"]
