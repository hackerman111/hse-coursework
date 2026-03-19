"""
Фабричные функции для создания дифференцирований.
"""

from __future__ import annotations

from typing import Any

from sage.all import Matrix


def _images_from_mapping(algebra: Any, gen_mapping: Any = None) -> list[Any]:
    if gen_mapping is None:
        gen_mapping = {}

    if isinstance(gen_mapping, (list, tuple)):
        images = list(gen_mapping)
        expected = len(algebra.gens())
        if len(images) != expected:
            raise ValueError(
                f"Ожидалось {expected} образов генераторов, получено {len(images)}"
            )
        return images

    if isinstance(gen_mapping, dict):
        gens = algebra.gens()
        return [gen_mapping.get(gen, 0) for gen in gens]

    raise TypeError(
        f"gen_mapping должен быть словарём или списком, получен {type(gen_mapping)}"
    )


def from_mapping(algebra: Any, gen_mapping: Any = None):
    from .factory import LieDerivationFactory

    images = _images_from_mapping(algebra, gen_mapping)
    return LieDerivationFactory.create(algebra.derivation(images))


def from_linear(algebra: Any, matrix: Any):
    from .factory import LieDerivationFactory

    gens = algebra.gens()
    size = len(gens)

    if matrix.nrows() != size or matrix.ncols() != size:
        raise ValueError(
            f"Матрица должна быть {size}x{size}, "
            f"получено {matrix.nrows()}x{matrix.ncols()}"
        )

    images = []
    for row_index in range(size):
        value = 0
        for col_index, gen in enumerate(gens):
            value += matrix[row_index, col_index] * gen
        images.append(value)

    return LieDerivationFactory.create(algebra.derivation(images))


def from_weitzenbock(algebra: Any, matrix: Any):
    if not matrix.is_nilpotent():
        raise ValueError(
            "Матрица должна быть нильпотентной для дифференцирования Вайтценбёка"
        )
    return from_linear(algebra, matrix)


def from_jacobian(algebra: Any, polynomials: Any):
    from .factory import LieDerivationFactory

    gens = algebra.gens()
    size = len(gens)
    if len(polynomials) != size - 1:
        raise ValueError(
            f"Требуется n-1 полиномов ({size - 1}), получено {len(polynomials)}."
        )

    full_jacobian = []
    for polynomial in polynomials:
        full_jacobian.append([polynomial.derivative(gen) for gen in gens])

    images = []
    for omitted_index in range(size):
        sub_matrix_data = []
        for row in full_jacobian:
            sub_matrix_data.append(row[:omitted_index] + row[omitted_index + 1 :])

        determinant = Matrix(sub_matrix_data).determinant()
        sign = (-1) ** (size + omitted_index + 1)
        images.append(sign * determinant)

    return LieDerivationFactory.create(algebra.derivation(images))
