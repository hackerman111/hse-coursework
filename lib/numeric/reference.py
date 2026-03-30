"""
Small reference checks tied to documented numeric examples.
"""

from __future__ import annotations

from .indexing import dim_Wn_leq_d
from ..utils import LibraryLogger


NUMERIC_MD_DIMENSION_EXAMPLES = [
    (2, 5, 56),
    (2, 10, 156),
    (3, 5, 252),
    (3, 10, 1092),
    (5, 5, 2310),
    (5, 10, 21840),
]


def check_dimension_examples(logger: LibraryLogger | None = None) -> bool:
    """
    Сверить табличные размерности из документации с аналитической формулой.

    На вход принимает необязательный ``LibraryLogger``.
    На выходе возвращает ``True``, если все эталонные значения совпали с ``dim_Wn_leq_d``.

    >>> check_dimension_examples()
    True
    """
    logger = logger or LibraryLogger()
    logger.line("Проверка таблицы размерностей из Numeric.md:", indent=0)
    all_ok = True
    for n, d, expected in NUMERIC_MD_DIMENSION_EXAMPLES:
        actual = dim_Wn_leq_d(n, d)
        ok = "✓" if actual == expected else "✗"
        if actual != expected:
            all_ok = False
        logger.info(f"{ok} dim W_{n}^(≤{d}) = {actual} (ожидалось {expected})")
    logger.line()
    return all_ok
