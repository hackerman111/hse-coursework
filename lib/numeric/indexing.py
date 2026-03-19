"""
Индексация мономов и базисных элементов W_n^{(k)}.

Все мультииндексы — кортежи (tuple) длины n.
Результаты кешируются для повторного использования.
"""

from __future__ import annotations

from functools import lru_cache
from math import comb
from typing import Dict, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Перечисление мультииндексов
# ---------------------------------------------------------------------------

@lru_cache(maxsize=512)
def multiindices(n: int, s: int) -> Tuple[Tuple[int, ...], ...]:
    """
    Все мультииндексы α ∈ Z≥0^n с |α| = s, упорядоченные лексикографически
    (по убыванию). Возвращает кортеж кортежей (неизменяемый, кешируемый).

    >>> multiindices(2, 2)
    ((2, 0), (1, 1), (0, 2))
    """
    if n == 1:
        return ((s,),)
    result: List[Tuple[int, ...]] = []
    for first in range(s, -1, -1):
        for tail in multiindices(n - 1, s - first):
            result.append((first,) + tail)
    return tuple(result)


@lru_cache(maxsize=512)
def multiindex_lookup(n: int, s: int) -> Dict[Tuple[int, ...], int]:
    """
    Обратный словарь: мультииндекс → его порядковый номер в `multiindices(n, s)`.
    """
    return {alpha: idx for idx, alpha in enumerate(multiindices(n, s))}


# ---------------------------------------------------------------------------
# Размерности
# ---------------------------------------------------------------------------

def num_monomials(n: int, s: int) -> int:
    """Число мономов степени s в n переменных = C(s + n - 1, n - 1)."""
    if s < 0:
        return 0
    return comb(s + n - 1, n - 1)


def dim_Wn_k(n: int, k: int) -> int:
    """
    dim W_n^{(k)} = n · C(k + n, n - 1) для k ≥ -1.
    Однородная компонента степени k алгебры W_n.
    """
    if k < -1:
        return 0
    return n * comb(k + n, n - 1)


def dim_Wn_leq_d(n: int, d: int) -> int:
    """
    dim W_n^{(≤d)} = n · C(d + n + 1, n).
    Суммарная размерность компонент степеней от -1 до d.
    """
    if d < -1:
        return 0
    return n * comb(d + n + 1, n)


# ---------------------------------------------------------------------------
# Индексация базиса W_n^{(k)}
# ---------------------------------------------------------------------------

def basis_index(alpha: Tuple[int, ...], i: int, n: int, k: int) -> int:
    """
    Индекс базисного элемента z^α ∂_{z_i} внутри W_n^{(k)}.

    Параметры
    ---------
    alpha : мультииндекс с |α| = k + 1
    i     : номер переменной (0-based)
    n     : число переменных
    k     : степень компоненты

    Возвращает
    ----------
    Индекс в диапазоне [0, dim_Wn_k(n, k)).
    """
    s = k + 1  # |α| должен быть s
    lookup = multiindex_lookup(n, s)
    m = len(lookup)  # = num_monomials(n, s)
    return i * m + lookup[alpha]


def decode_basis_index(idx: int, n: int, k: int) -> Tuple[Tuple[int, ...], int]:
    """
    Обратная операция: по индексу вернуть (α, i).
    """
    s = k + 1
    monoms = multiindices(n, s)
    m = len(monoms)
    i = idx // m
    alpha = monoms[idx % m]
    return alpha, i


# ---------------------------------------------------------------------------
# Вспомогательные numpy-утилиты
# ---------------------------------------------------------------------------

def unit_vector(n: int) -> np.ndarray:
    """Стандартные единичные вектора e_0, …, e_{n-1} как строки."""
    return np.eye(n, dtype=np.int64)


_UNIT_CACHE: Dict[int, np.ndarray] = {}


def e_vec(i: int, n: int) -> np.ndarray:
    """Единичный вектор e_i длины n (int64, кешированный)."""
    if n not in _UNIT_CACHE:
        _UNIT_CACHE[n] = np.eye(n, dtype=np.int64)
    return _UNIT_CACHE[n][i]
