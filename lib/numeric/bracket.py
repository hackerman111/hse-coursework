"""
Вычисление скобки Ли в W_n через разреженные структурные константы.

Формула:
  [z^α ∂_{z_i}, z^β ∂_{z_j}] = β_i · z^{α+β-e_i} ∂_{z_j}
                                - α_j · z^{α+β-e_j} ∂_{z_i}

Ключевое свойство: каждая пара базисных элементов даёт ≤2 ненулевых
компонент → тензор экстремально разрежен.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Dict, List, Tuple

import numpy as np

from .indexing import (
    multiindices,
    multiindex_lookup,
    dim_Wn_k,
    num_monomials,
)


# Тип одной ненулевой структурной константы: (a, b, c, λ)
# Здесь a — индекс в W_n^{(p)}, b — индекс в W_n^{(q)}, c — индекс в W_n^{(r)}, λ — коэффициент
StructConst = Tuple[int, int, int, int]


@lru_cache(maxsize=256)
def precompute_structure_constants(
    n: int, p: int, q: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Предвычислить разреженное представление структурных констант для ``[W_n^(p), W_n^(q)]``.

    На вход принимает размерность ``n`` и степени ``p, q``.
    На выходе возвращает четыре numpy-массива ``(row_a, row_b, col_c, vals)``, кодирующие ненулевые коэффициенты скобки.

    >>> tuple(arr.shape[0] for arr in precompute_structure_constants(2, -1, -1))
    (0, 0, 0, 0)
    """
    r = p + q
    s_p = p + 1  # |α| для W_n^{(p)}
    s_q = q + 1  # |β| для W_n^{(q)}
    s_r = r + 1  # |γ| для W_n^{(r)}

    monoms_p = multiindices(n, s_p)
    monoms_q = multiindices(n, s_q)
    lookup_r = multiindex_lookup(n, s_r)

    m_p = len(monoms_p)  # число мономов степени s_p
    m_q = len(monoms_q)
    m_r = len(lookup_r)

    rows_a: List[int] = []
    rows_b: List[int] = []
    cols_c: List[int] = []
    vals: List[int] = []

    alpha_arr = np.array(monoms_p, dtype=np.int64)  # размер (m_p, n)
    beta_arr = np.array(monoms_q, dtype=np.int64)  # размер (m_q, n)

    for ai, alpha in enumerate(monoms_p):
        alpha_np = alpha_arr[ai]
        for bi, beta in enumerate(monoms_q):
            beta_np = beta_arr[bi]
            for i in range(n):
                a_idx = i * m_p + ai  # индекс z^α ∂_{z_i} в W_n^{(p)}
                for j in range(n):
                    b_idx = j * m_q + bi  # индекс z^β ∂_{z_j} в W_n^{(q)}

                    # Первый член: β_i · z^{α+β-e_i} ∂_{z_j}
                    if beta_np[i] > 0:
                        gamma = tuple(
                            alpha_np[k] + beta_np[k] - (1 if k == i else 0)
                            for k in range(n)
                        )
                        c_idx = j * m_r + lookup_r[gamma]
                        rows_a.append(a_idx)
                        rows_b.append(b_idx)
                        cols_c.append(c_idx)
                        vals.append(int(beta_np[i]))

                    # Второй член: -α_j · z^{α+β-e_j} ∂_{z_i}
                    if alpha_np[j] > 0:
                        gamma = tuple(
                            alpha_np[k] + beta_np[k] - (1 if k == j else 0)
                            for k in range(n)
                        )
                        c_idx = i * m_r + lookup_r[gamma]
                        rows_a.append(a_idx)
                        rows_b.append(b_idx)
                        cols_c.append(c_idx)
                        vals.append(-int(alpha_np[j]))

    return (
        np.array(rows_a, dtype=np.int64),
        np.array(rows_b, dtype=np.int64),
        np.array(cols_c, dtype=np.int64),
        np.array(vals, dtype=np.float64),
    )


def bracket_vectors(
    u: np.ndarray,
    p: int,
    v: np.ndarray,
    q: int,
    n: int,
    sc: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> np.ndarray:
    """
    Вычислить скобку Ли двух однородных координатных векторов.

    На вход принимает векторы ``u`` и ``v`` компонент степеней ``p`` и ``q``, размерность ``n`` и необязательный кеш структурных констант ``sc``.
    На выходе возвращает координатный вектор скобки в ``W_n^(p+q)``.

    >>> u = np.array([1.0, 0.0])
    >>> v = np.array([0.0, 1.0])
    >>> bracket_vectors(u, -1, v, -1, 2).tolist()
    [0.0, 0.0]
    """
    r = p + q
    if sc is None:
        sc = precompute_structure_constants(n, p, q)

    rows_a, rows_b, cols_c, vals = sc
    N_r = dim_Wn_k(n, r)
    result = np.zeros(N_r, dtype=np.float64)

    if len(rows_a) == 0:
        return result

    # Векторизованное накопление: result[c] += u[a] * v[b] * λ
    contributions = u[rows_a] * v[rows_b] * vals
    np.add.at(result, cols_c, contributions)
    return result


def bracket_full_elements(
    f: Dict[int, np.ndarray],
    g: Dict[int, np.ndarray],
    n: int,
    structure_cache: "StructureConstantCache",
) -> Dict[int, np.ndarray]:
    """
    Вычислить скобку Ли двух полных неоднородных элементов.

    На вход принимает словари компонент ``f`` и ``g``, размерность ``n`` и ``structure_cache``.
    На выходе возвращает словарь ``degree -> np.ndarray`` для результирующего элемента, включая только ненулевые степени.

    >>> from .structure import StructureConstantCache
    >>> cache = StructureConstantCache(n=2, D_max=0)
    >>> bracket_full_elements({-1: np.array([1.0, 0.0])}, {-1: np.array([0.0, 1.0])}, 2, cache)
    {}
    """
    result: Dict[int, np.ndarray] = {}

    for p, f_p in f.items():
        for q, g_q in g.items():
            r = p + q
            sc = structure_cache.get(p, q) if p <= q else structure_cache.get(q, p)
            if sc is None:
                continue

            if p <= q:
                contribution = bracket_vectors(f_p, p, g_q, q, n, sc=sc)
            else:
                # При использовании sc для (q, p) применяем антисимметрию скобки
                contribution = -bracket_vectors(g_q, q, f_p, p, n, sc=sc)

            N_r = dim_Wn_k(n, r)
            if r not in result:
                result[r] = np.zeros(N_r, dtype=np.float64)
            result[r] += contribution

    # Отбрасываем нулевые компоненты
    return {k: v for k, v in result.items() if np.linalg.norm(v) > 1e-15}


def batch_bracket(
    Bp: np.ndarray,
    p: int,
    Bq: np.ndarray,
    q: int,
    n: int,
    sc: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> np.ndarray:
    """
    Пакетно вычислить скобки между всеми строками двух матриц базиса.

    На вход принимает матрицы ``Bp`` и ``Bq`` с базисными строками степеней ``p`` и ``q``, размерность ``n`` и необязательный кеш ``sc``.
    На выходе возвращает матрицу, каждая строка которой равна одной скобке ``[Bp_s, Bq_t]``.

    >>> Bp = np.array([[1.0, 0.0]])
    >>> Bq = np.array([[0.0, 1.0]])
    >>> batch_bracket(Bp, -1, Bq, -1, 2).shape
    (0, 2)
    """
    r = p + q
    if sc is None:
        sc = precompute_structure_constants(n, p, q)

    rows_a, rows_b, cols_c, vals = sc
    N_r = dim_Wn_k(n, r)
    rp = Bp.shape[0]
    rq = Bq.shape[0]

    if len(rows_a) == 0 or rp == 0 or rq == 0:
        return np.empty((0, N_r), dtype=np.float64)

    # Для каждой пары (s, t) вычисляем вклад в результирующую строку
    # Bp[:, rows_a] имеет размер (rp, nnz), Bq[:, rows_b] имеет размер (rq, nnz)
    A = Bp[:, rows_a]  # размер (rp, nnz)
    B = Bq[:, rows_b]  # размер (rq, nnz)

    # Вместо полного тензорного произведения работаем построчно ради экономии памяти
    results = np.zeros((rp * rq, N_r), dtype=np.float64)
    for s in range(rp):
        a_row = A[s]  # (nnz,)
        for t in range(rq):
            contrib = a_row * B[t] * vals  # (nnz,)
            row_idx = s * rq + t
            np.add.at(results[row_idx], cols_c, contrib)

    return results
