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
from typing import List, Tuple

import numpy as np

from .indexing import (
    multiindices,
    multiindex_lookup,
    dim_Wn_k,
    num_monomials,
)


# Тип одной ненулевой структурной константы: (a, b, c, λ)
# a — индекс в W_n^{(p)}, b — индекс в W_n^{(q)}, c — индекс в W_n^{(r)}, λ — коэффициент
StructConst = Tuple[int, int, int, int]


@lru_cache(maxsize=256)
def precompute_structure_constants(
    n: int, p: int, q: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Предвычислить разреженный тензор структурных констант C^{(p,q)}.

    Возвращает 4 массива (row_a, row_b, col_c, vals) такие что
      [e_a^{(p)}, e_b^{(q)}]_c = vals[k]   для k-го ненулевого элемента.

    Результат кешируется.
    """
    r = p + q
    s_p = p + 1  # |α| для W_n^{(p)}
    s_q = q + 1  # |β| для W_n^{(q)}
    s_r = r + 1  # |γ| для W_n^{(r)}

    monoms_p = multiindices(n, s_p)
    monoms_q = multiindices(n, s_q)
    lookup_r = multiindex_lookup(n, s_r)

    m_p = len(monoms_p)  # num_monomials(n, s_p)
    m_q = len(monoms_q)
    m_r = len(lookup_r)

    rows_a: List[int] = []
    rows_b: List[int] = []
    cols_c: List[int] = []
    vals: List[int] = []

    alpha_arr = np.array(monoms_p, dtype=np.int64)  # shape (m_p, n)
    beta_arr = np.array(monoms_q, dtype=np.int64)    # shape (m_q, n)

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
    Вычислить [u, v] для u ∈ W_n^{(p)}, v ∈ W_n^{(q)}.

    u, v — одномерные numpy-массивы (координаты в базисе компоненты).
    Возвращает вектор в W_n^{(p+q)}.
    """
    r = p + q
    if sc is None:
        sc = precompute_structure_constants(n, p, q)

    rows_a, rows_b, cols_c, vals = sc
    N_r = dim_Wn_k(n, r)
    result = np.zeros(N_r, dtype=np.float64)

    if len(rows_a) == 0:
        return result

    # Vectorised: result[c] += u[a] * v[b] * λ
    contributions = u[rows_a] * v[rows_b] * vals
    np.add.at(result, cols_c, contributions)
    return result


def batch_bracket(
    Bp: np.ndarray,
    p: int,
    Bq: np.ndarray,
    q: int,
    n: int,
    sc: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None = None,
) -> np.ndarray:
    """
    Вычислить все скобки [row_s(Bp), row_t(Bq)].

    Bp : shape (rp, Np) — строки базиса L^{(p)}
    Bq : shape (rq, Nq) — строки базиса L^{(q)}

    Возвращает матрицу shape (rp * rq, Nr) со строками = [Bp_s, Bq_t].
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

    # Для каждой пары (s, t) вычисляем contributions
    # Bp[:, rows_a] shape (rp, nnz), Bq[:, rows_b] shape (rq, nnz)
    A = Bp[:, rows_a]  # (rp, nnz)
    B = Bq[:, rows_b]  # (rq, nnz)

    # outer: (rp, rq, nnz), масштабировано на vals
    # result[s, t, c] = sum over matching k : A[s, k] * B[t, k] * vals[k]
    # Вместо полного outer-product работаем построчно для экономии памяти
    results = np.zeros((rp * rq, N_r), dtype=np.float64)
    for s in range(rp):
        a_row = A[s]  # (nnz,)
        for t in range(rq):
            contrib = a_row * B[t] * vals  # (nnz,)
            row_idx = s * rq + t
            np.add.at(results[row_idx], cols_c, contrib)

    return results
