"""
Фабричные функции для создания дифференцирований.
"""

from __future__ import annotations

from typing import Any


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
    """
    Построить дифференцирование по явному отображению генераторов.

    На вход принимает Sage-алгебру ``algebra`` и либо словарь, либо список образов ``gen_mapping``.
    На выходе возвращает объект ``LieDerivation``, созданный фабрикой.

    >>> from sage.all import PolynomialRing, QQ
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> type(from_mapping(ring, {ring.gen(0): 1, ring.gen(1): 0})).__name__
    'LieDerivation'
    """
    from .factory import LieDerivationFactory

    images = _images_from_mapping(algebra, gen_mapping)
    return LieDerivationFactory.create(algebra.derivation(images))


def from_linear(algebra: Any, matrix: Any):
    """
    Построить линейное дифференцирование по матрице коэффициентов.

    На вход принимает Sage-алгебру ``algebra`` и квадратную матрицу ``matrix`` размера ``n x n``.
    На выходе возвращает ``LieDerivation``, действие которого на генераторах задаётся этой матрицей.

    >>> from sage.all import PolynomialRing, QQ, matrix
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> derivation = from_linear(ring, matrix(QQ, [[1, 0], [0, 0]]))
    >>> type(derivation).__name__
    'LieDerivation'
    """
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
    """
    Построить дифференцирование Вайтценбёка по нильпотентной матрице.

    На вход принимает Sage-алгебру ``algebra`` и нильпотентную матрицу ``matrix``.
    На выходе возвращает линейное ``LieDerivation``, индуцированное этой матрицей.

    >>> from sage.all import PolynomialRing, QQ, matrix
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> derivation = from_weitzenbock(ring, matrix(QQ, [[0, 1], [0, 0]]))
    >>> type(derivation).__name__
    'LieDerivation'
    """
    if not matrix.is_nilpotent():
        raise ValueError(
            "Матрица должна быть нильпотентной для дифференцирования Вайтценбёка"
        )
    return from_linear(algebra, matrix)


def from_jacobian(algebra: Any, polynomials: Any):
    """
    Построить дивергентно-свободное дифференцирование по формуле Якобиана.

    На вход принимает Sage-алгебру ``algebra`` и список из ``n-1`` полиномов ``polynomials``.
    На выходе возвращает ``LieDerivation``, чьи компоненты выражаются через миноры матрицы Якоби.

    >>> from sage.all import PolynomialRing, QQ
    >>> ring = PolynomialRing(QQ, "x,y")
    >>> derivation = from_jacobian(ring, [ring.gen(0)])
    >>> type(derivation).__name__
    'LieDerivation'
    """
    from .factory import LieDerivationFactory
    from sage.all import Matrix

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
