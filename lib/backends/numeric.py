"""
Numeric backend: конвертация DerivationSpec в NumPy-представление.

Не зависит от SageMath. Требует NumPy.
"""

from __future__ import annotations

import numpy as np

from lib.core.spec import DerivationSpec, to_generator


def to_numeric_components(
    spec: DerivationSpec,
    D_max: int,
) -> dict[int, np.ndarray]:
    """
    Конвертация DerivationSpec в разложение по однородным компонентам.

    Возвращает dict: степень → координатный вектор в W_n^(k).
    Используется как входной формат для NumericLieGenerationChecker.

    Args:
        spec:  дифференцирование
        D_max: максимальная степень (включительно)

    Returns:
        {degree: np.ndarray} — разложение по однородным слоям
    """
    from lib.numeric.components import decompose_generator

    generator = to_generator(spec)
    return decompose_generator(generator, spec.n, D_max)
