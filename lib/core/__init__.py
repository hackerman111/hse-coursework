"""
Ядро библиотеки: каноническая модель DerivationSpec и кирпичики.
"""

from .spec import (
    DerivationSpec,
    LeadingTerm,
    Scalar,
    combine_specs,
    homogeneous_component,
    homogeneous_components,
    leading_term_spec,
    monomial_spec,
    normalize_spec,
    partial_spec,
    random_spec,
    scale_spec,
    spec_degree,
    zero_spec,
)

__all__ = [
    "DerivationSpec",
    "LeadingTerm",
    "Scalar",
    "monomial_spec",
    "partial_spec",
    "combine_specs",
    "spec_degree",
    "random_spec",
    "zero_spec",
    "scale_spec",
    "normalize_spec",
    "homogeneous_component",
    "homogeneous_components",
    "leading_term_spec",
]
