"""Тесты public API верхнего уровня и symbolic library wrappers."""

from __future__ import annotations

import lib
import lib.numeric as numeric
import pytest

from lib.core.spec import partial_spec


def test_top_level_lib_reexports_numeric_public_api():
    for name in numeric.__all__:
        assert getattr(lib, name) is getattr(numeric, name)


def test_symbolic_module_exposes_library_entry_points():
    import lib.symbolic as symbolic

    expected = {
        "SymbolicCriterionResult",
        "check_symbolic_degree2_criterion",
        "check_symbolic_generation",
        "format_target",
        "make_degree2_stop_condition",
        "make_symbolic_algebra",
        "spec_to_symbolic_derivation",
    }

    for name in expected:
        assert hasattr(symbolic, name)


def test_check_symbolic_generation_converts_specs_before_solver(monkeypatch):
    import lib.symbolic as symbolic

    specs = [partial_spec(2, axis=0), partial_spec(2, axis=1)]
    calls: dict[str, object] = {}

    class FakeSolver:
        def __init__(self, generators, max_iter=1000, stop_condition=None):
            calls["generators"] = list(generators)
            calls["max_iter"] = max_iter
            calls["stop_condition"] = stop_condition
            self.iter_count = 4
            self.basis = {"lt0": object(), "lt1": object()}
            self.targets_found = {0: True, 1: True}

        def run(self) -> bool:
            return True

    monkeypatch.setattr(symbolic, "LieBasisSolver", FakeSolver)
    monkeypatch.setattr(symbolic, "make_symbolic_algebra", lambda n: f"ALG:{n}")
    monkeypatch.setattr(
        symbolic,
        "spec_to_symbolic_derivation",
        lambda spec, algebra: f"lie:{spec.n}:{algebra}",
    )

    result = symbolic.check_symbolic_generation(specs, max_iter=17)

    assert result.success
    assert result.iterations == 4
    assert result.basis_size == 2
    assert calls["max_iter"] == 17
    assert calls["stop_condition"] is None
    assert calls["generators"] == ["lie:2:ALG:2", "lie:2:ALG:2"]


def test_check_symbolic_degree2_criterion_formats_missing_targets(monkeypatch):
    import lib.symbolic as symbolic

    specs = [partial_spec(2, axis=0), partial_spec(2, axis=1)]
    sentinel_stop = object()

    class FakeSolver:
        def __init__(self, generators, max_iter=1000, stop_condition=None):
            self.iter_count = 6
            self.basis = {}
            self.targets_found = {}
            self.stop_condition = stop_condition

        def run(self) -> bool:
            return False

    monkeypatch.setattr(symbolic, "LieBasisSolver", FakeSolver)
    monkeypatch.setattr(symbolic, "make_symbolic_algebra", lambda n: f"ALG:{n}")
    monkeypatch.setattr(
        symbolic,
        "spec_to_symbolic_derivation",
        lambda spec, algebra: f"lie:{spec.n}:{algebra}",
    )
    monkeypatch.setattr(
        symbolic,
        "make_degree2_stop_condition",
        lambda algebra: (
            sentinel_stop,
            {
                (0, "1"): True,
                (1, "z0"): False,
            },
        ),
    )
    monkeypatch.setattr(
        symbolic,
        "format_target",
        lambda axis, monomial, algebra: f"T({axis},{monomial},{algebra})",
    )

    result = symbolic.check_symbolic_degree2_criterion(specs, max_iter=23)

    assert not result.success
    assert result.iterations == 6
    assert result.basis_size == 0
    assert result.found_targets == ("T(0,1,ALG:2)",)
    assert result.missing_targets == ("T(1,z0,ALG:2)",)
