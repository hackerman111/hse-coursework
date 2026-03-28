"""
Numeric backend: конвертация DerivationSpec в NumPy-представление.

Не зависит от SageMath. Требует NumPy.
"""

from __future__ import annotations

import numpy as np

from lib.core.spec import DerivationSpec


def to_generator(spec: DerivationSpec) -> "dict[int, dict[tuple[int, ...], float]]":
    """
    Конвертация DerivationSpec в формат Generator (dict[int, dict[tuple, float]]).

    На вход принимает `DerivationSpec`.
    На выходе возвращает словарь старого numeric-формата с коэффициентами `float`.

    >>> from lib.core.spec import partial_spec
    >>> to_generator(partial_spec(2, 0))
    {0: {(0, 0): 1.0}}
    """
    return {
        axis: {alpha: float(coeff) for alpha, coeff in mono_dict.items()}
        for axis, mono_dict in spec.terms.items()
    }


def to_numeric_components(
    spec: DerivationSpec,
    D_max: int,
) -> dict[int, np.ndarray]:
    """
    Конвертация DerivationSpec в разложение по однородным компонентам.

    На вход принимает `DerivationSpec` и максимальную степень `D_max`
    для numeric backend.
    На выходе возвращает словарь `степень -> numpy-вектор координат`.

    >>> from lib.core.spec import partial_spec
    >>> components = to_numeric_components(partial_spec(2, 0), D_max=2)
    >>> sorted(components)
    [-1]
    >>> components[-1].tolist()
    [1.0, 0.0]
    """
    from lib.numeric.components import decompose_generator

    generator = to_generator(spec)
    return decompose_generator(generator, spec.n, D_max)


def from_numeric_components(
    n: int,
    components: dict[int, np.ndarray],
) -> DerivationSpec:
    """
    Обратное преобразование из однородных numeric-компонент в DerivationSpec.

    На вход принимает размерность `n` и словарь `степень -> numpy-вектор`.
    На выходе возвращает восстановленный `DerivationSpec`.

    >>> from_numeric_components(2, {-1: np.array([1.0, 0.0])})
    DerivationSpec(n=2, terms={0: {(0, 0): 1.0}})
    """
    from lib.numeric.indexing import decode_basis_index, dim_Wn_k

    terms: dict[int, dict[tuple[int, ...], float]] = {}

    for degree, vector in components.items():
        flat = np.asarray(vector, dtype=np.float64).reshape(-1)
        expected_dim = dim_Wn_k(n, degree)
        if flat.shape != (expected_dim,):
            raise ValueError(
                f"Компонента степени {degree} имеет размер {flat.shape[0]}, "
                f"ожидалось {expected_dim}"
            )

        for index, coeff in enumerate(flat):
            if abs(float(coeff)) <= 1e-15:
                continue
            alpha, axis = decode_basis_index(index, n, degree)
            axis_terms = terms.setdefault(axis, {})
            axis_terms[alpha] = axis_terms.get(alpha, 0.0) + float(coeff)

    return DerivationSpec(n=n, terms=terms)
