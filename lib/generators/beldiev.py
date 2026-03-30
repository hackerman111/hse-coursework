"""
Генераторы по теореме Бельдиева — единственное каноническое определение.

Формулы Бельдиева для W_n (2 генератора):

    U = ∂/∂z_{n-1}

    V_k (компонента ∂/∂z_k, k = 0..n-2):
        V(z^alpha * ∂/∂z_k) где alpha_j = 4*(k+2) для j = k+1..n-1, остальные 0.

    V_{n-1}:
        z_0^4 * z_1^4 * ... * z_{n-1}^4 * ∂/∂z_{n-1}

Источник: статья Бельдиева. Ранее дублировалось в:
- lib/generator/beldiev.py (Sage-версия)
- lib/numeric/recipes.py (numeric-версия)

Теперь оба слоя используют этот модуль.
"""

from __future__ import annotations

from lib.core.spec import DerivationSpec, combine_specs, monomial_spec, partial_spec


def beldiev_specs(n: int) -> list[DerivationSpec]:
    """
    Генераторы Бельдиева [U, V] для W_n как список DerivationSpec.

    На вход принимает размерность `n >= 2`.
    На выходе возвращает список `[U, V]` в каноническом spec-формате.

    >>> from lib.core.spec import spec_degree
    >>> specs = beldiev_specs(2)
    >>> len(specs), specs[0] == partial_spec(2, axis=1), spec_degree(specs[1])
    (2, True, 7)
    """
    if n < 2:
        raise ValueError(f"Теорема Бельдиева требует n >= 2, получено n={n}")

    # Генератор U равен ∂/∂z_{n-1}
    u = partial_spec(n, axis=n - 1)

    # Генератор V собирается из (n - 1) + 1 слагаемых
    v_parts: list[DerivationSpec] = []

    # Компоненты V_k для k = 0..n-2:
    # образ имеет вид (z_{k+1} * ... * z_{n-1})^{4*(k+2)} * ∂/∂z_k
    for k in range(n - 1):
        power = 4 * (k + 2)
        alpha = [0] * n
        for j in range(k + 1, n):
            alpha[j] = power
        v_parts.append(monomial_spec(n, axis=k, alpha=tuple(alpha), coeff=1))

    # Компонента V_{n-1}: (z_0 * ... * z_{n-1})^4 * ∂/∂z_{n-1}
    alpha_last = (4,) * n
    v_parts.append(monomial_spec(n, axis=n - 1, alpha=alpha_last, coeff=1))

    v = combine_specs(*v_parts)
    return [u, v]
