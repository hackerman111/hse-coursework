"""
Проверки на наличие doctest у публичных функций API.
"""

from __future__ import annotations

import inspect

import lib


EXCLUDED_PUBLIC_NAMES = {
    "ABC",
    "Derivation",
    "DerivationSpec",
    "LeadingTerm",
    "LibraryLogger",
    "LieBasisSolver",
    "LieDerivation",
    "LieDerivationFactory",
    "PolynomialRing",
    "QuotientLieDerivation",
    "abstractmethod",
}


def _public_api_functions():
    for name in lib.__all__:
        if name in EXCLUDED_PUBLIC_NAMES:
            continue
        obj = getattr(lib, name)
        if inspect.isfunction(obj):
            yield name, obj


def test_every_public_api_function_has_russian_doctest():
    missing = []

    for name, func in _public_api_functions():
        doc = inspect.getdoc(func) or ""
        has_doctest = ">>>" in doc
        has_input_line = "На вход" in doc
        has_output_line = "На выход" in doc
        if not (has_doctest and has_input_line and has_output_line):
            missing.append(name)

    assert not missing, (
        "Для публичных API-функций нужны docstring с doctest и краткими "
        f"описаниями входа/выхода. Нет для: {', '.join(sorted(missing))}"
    )
