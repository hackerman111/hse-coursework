"""
Генераторы по теореме Андриста — единственное каноническое определение.

Формулы Андриста для W_n (3 генератора):

    U = ∂/∂z_{n-1}

    V:
        - компонента ∂/∂z_{n-1}: константа 1 (т.е. ∂/∂z_{n-1})
        - компонента ∂/∂z_k (k = 0..n-2):
            z_{k+1}^3 * z_{k+2} * ... * z_{n-1} * ∂/∂z_k

    W:
        - компонента ∂/∂z_{n-1}: z_0^2 * z_1^2 * ... * z_{n-2}^2 * z_{n-1}

Источник: статья Андриста. Ранее определялось только в lib/generator/andrist.py (Sage).
"""

from __future__ import annotations

from lib.core.spec import DerivationSpec, combine_specs, monomial_spec, partial_spec


def andrist_specs(n: int) -> list[DerivationSpec]:
    """
    Генераторы Андриста [U, V, W] для W_n как список DerivationSpec.

    Args:
        n: размерность пространства (n >= 2)

    Returns:
        Список из трёх DerivationSpec: [U, V, W].
    """
    if n < 2:
        raise ValueError(f"Теорема Андриста требует n >= 2, получено n={n}")

    # U = ∂/∂z_{n-1}
    u = partial_spec(n, axis=n - 1)

    # V: собирается из n слагаемых
    v_parts: list[DerivationSpec] = []

    # Константная компонента V_{n-1}: 1 * ∂/∂z_{n-1}
    v_parts.append(monomial_spec(n, axis=n - 1, alpha=(0,) * n, coeff=1))

    # Полиномиальные компоненты V_k для k = 0..n-2:
    # image = z_{k+1}^3 * z_{k+2} * ... * z_{n-1} * ∂/∂z_k
    for k in range(n - 1):
        alpha = [0] * n
        alpha[k + 1] = 3
        for j in range(k + 2, n):
            alpha[j] = 1
        v_parts.append(monomial_spec(n, axis=k, alpha=tuple(alpha), coeff=1))

    v = combine_specs(*v_parts)

    # W: z_0^2 * ... * z_{n-2}^2 * z_{n-1} * ∂/∂z_{n-1}
    alpha_w = [0] * n
    for j in range(n - 1):
        alpha_w[j] = 2
    alpha_w[n - 1] = 1
    w = monomial_spec(n, axis=n - 1, alpha=tuple(alpha_w), coeff=1)

    return [u, v, w]
