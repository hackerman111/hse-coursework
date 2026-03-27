"""
Named generator families used by numeric experiments.

Делегирует в lib/generators/ — единственный источник формул.
"""

from __future__ import annotations

from lib.core.spec import to_generator
from lib.generators.beldiev import beldiev_specs
from .types import Generator


def make_beldiev_generators(n: int) -> list[Generator]:
    """
    Beldiev's explicit 2-generator family for ``W_n``.

    Делегирует в beldiev_specs() — каноническое определение формул.
    """
    return [to_generator(s) for s in beldiev_specs(n)]
