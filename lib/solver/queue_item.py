"""
Элементы очереди для solver.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(order=True)
class QueueItem:
    degree: int
    unique_counter: int
    derivation: Any = field(compare=False)
