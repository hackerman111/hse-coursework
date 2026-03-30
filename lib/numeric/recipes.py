"""
Named generator families used by numeric experiments.

Делегирует в lib/generators/ — единственный источник формул.
"""

from __future__ import annotations

from lib.backends.numeric import to_generator
from lib.generators.beldiev import beldiev_specs
from .types import Generator


def make_beldiev_generators(n: int) -> list[Generator]:
    """
    Построить численные генераторы семейства Бельдиева для ``W_n``.

    На вход принимает размерность ``n``.
    На выходе возвращает список из двух generator-объектов, полученных из канонических ``DerivationSpec``.

    >>> len(make_beldiev_generators(2))
    2
    """
    return [to_generator(s) for s in beldiev_specs(n)]
