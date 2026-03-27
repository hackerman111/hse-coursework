"""
Ядро библиотеки: каноническая модель DerivationSpec и кирпичики.
"""

from .spec import (
    DerivationSpec,
    Scalar,
    monomial_spec,
    partial_spec,
    combine_specs,
    spec_degree,
    random_spec,
    to_generator,
)

__all__ = [
    "DerivationSpec",
    "Scalar",
    "monomial_spec",
    "partial_spec",
    "combine_specs",
    "spec_degree",
    "random_spec",
    "to_generator",
]
