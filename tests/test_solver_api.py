"""
Тесты для lib/solver/__init__.py: unified check_generation() entry point
и standalone примитивов из lib/numeric/solver.py.
"""

import pytest
import numpy as np

from lib.core.spec import partial_spec, monomial_spec, combine_specs
from lib.core.results import NumericResult
from lib.generators.beldiev import beldiev_specs
from lib.solver import check_generation
from lib.numeric.solver import try_add_to_basis, check_degree_fullness
from lib.numeric.subspace import OrthonormalBasis


class TestCheckGenerationWithSpecs:
    """check_generation с list[DerivationSpec]."""

    def test_beldiev_n2_d2_success(self):
        specs = beldiev_specs(2)
        result = check_generation(specs, mode="numeric", d=2)
        assert isinstance(result, NumericResult)
        assert result.success

    def test_beldiev_n2_d5_success(self):
        specs = beldiev_specs(2)
        result = check_generation(specs, mode="numeric", d=5)
        assert result.success

    def test_spec_detects_failure(self):
        # Одиночный partial spec не порождает всю алгебру
        specs = [partial_spec(2, axis=0)]
        result = check_generation(specs, mode="numeric", d=2)
        assert not result.success

    def test_result_has_degrees(self):
        specs = beldiev_specs(2)
        result = check_generation(specs, mode="numeric", d=2)
        assert len(result.degrees) == 4  # -1, 0, 1, 2
        assert result.degrees[0].degree == -1

    def test_auto_mode_with_d_uses_numeric(self):
        specs = beldiev_specs(2)
        result = check_generation(specs, mode="auto", d=2)
        assert isinstance(result, NumericResult)


class TestCheckGenerationWithGeneratorDicts:
    """check_generation с list[Generator] (dict формат)."""

    def test_dict_generators_require_n(self):
        gen = {0: {(0, 0): 1.0}}
        with pytest.raises(ValueError, match="n="):
            check_generation([gen], mode="numeric", d=2)

    def test_dict_generators_with_n(self):
        from lib.backends.numeric import to_generator
        specs = beldiev_specs(2)
        gens = [to_generator(s) for s in specs]
        result = check_generation(gens, n=2, mode="numeric", d=2)
        assert isinstance(result, NumericResult)
        assert result.success

    def test_dict_partial_fails(self):
        gen = {0: {(0, 0): 1.0}}
        result = check_generation([gen], n=2, mode="numeric", d=2)
        assert not result.success


class TestCheckGenerationAutoMode:
    """Проверка логики auto-dispatching."""

    def test_auto_with_d_goes_numeric(self):
        specs = beldiev_specs(2)
        result = check_generation(specs, d=2)
        assert isinstance(result, NumericResult)

    def test_auto_spec_no_d_goes_symbolic(self):
        # DerivationSpec без d → symbolic (requires SageMath), пропускаем
        pytest.skip("auto mode с DerivationSpec без d идёт в symbolic mode (нужен SageMath)")

    def test_empty_generators_raises(self):
        with pytest.raises(ValueError, match="хотя бы один"):
            check_generation([])

    def test_unknown_mode_raises(self):
        specs = beldiev_specs(2)
        with pytest.raises(ValueError, match="Неизвестный"):
            check_generation(specs, mode="unknown")


class TestCheckGenerationParameters:
    """Параметры D_max и verbose."""

    def test_custom_D_max(self):
        specs = beldiev_specs(2)
        result = check_generation(specs, mode="numeric", d=2, D_max=4)
        assert isinstance(result, NumericResult)

    def test_verbose_false_no_crash(self):
        specs = beldiev_specs(2)
        result = check_generation(specs, mode="numeric", d=2, verbose=False)
        assert result is not None

    def test_result_n_and_d_match(self):
        specs = beldiev_specs(2)
        result = check_generation(specs, mode="numeric", d=3)
        assert result.n == 2
        assert result.d == 3


class TestTryAddToBasis:
    """Standalone примитив try_add_to_basis."""

    def test_adds_new_vector(self):
        basis = OrthonormalBasis(3)
        v = np.array([[1.0, 0.0, 0.0]])
        added = try_add_to_basis(basis, v)
        assert added == 1
        assert basis.rank == 1

    def test_skips_zero_vector(self):
        basis = OrthonormalBasis(3)
        v = np.array([[0.0, 0.0, 0.0]])
        added = try_add_to_basis(basis, v)
        assert added == 0
        assert basis.rank == 0

    def test_skips_dependent_vector(self):
        basis = OrthonormalBasis(3)
        v = np.array([[1.0, 0.0, 0.0]])
        try_add_to_basis(basis, v)
        added = try_add_to_basis(basis, v * 2)  # коллинеарный
        assert added == 0
        assert basis.rank == 1

    def test_adds_independent_vectors(self):
        basis = OrthonormalBasis(3)
        v1 = np.array([[1.0, 0.0, 0.0]])
        v2 = np.array([[0.0, 1.0, 0.0]])
        try_add_to_basis(basis, v1)
        added = try_add_to_basis(basis, v2)
        assert added == 1
        assert basis.rank == 2

    def test_empty_matrix(self):
        basis = OrthonormalBasis(3)
        v = np.zeros((0, 3))
        added = try_add_to_basis(basis, v)
        assert added == 0

    def test_normalizes_before_adding(self):
        basis = OrthonormalBasis(3)
        v = np.array([[100.0, 0.0, 0.0]])
        added = try_add_to_basis(basis, v)
        assert added == 1


class TestCheckDegreeFullness:
    """Standalone примитив check_degree_fullness."""

    def test_empty_bases_returns_zeros(self):
        bases = {}
        statuses = check_degree_fullness(bases, n=2, d=0)
        # degrees -1 and 0
        assert len(statuses) == 2
        for s in statuses:
            assert s.rank == 0
            assert not s.is_full

    def test_full_basis_marked_full(self):
        from lib.numeric.indexing import dim_Wn_k
        bases = {}
        for k in range(-1, 2):
            N_k = dim_Wn_k(2, k)
            b = OrthonormalBasis(N_k)
            # Заполним базис единичной матрицей
            eye = np.eye(N_k)
            b.try_add_many(eye)
            bases[k] = b
        statuses = check_degree_fullness(bases, n=2, d=1)
        for s in statuses:
            assert s.is_full

    def test_correct_degree_range(self):
        bases = {}
        statuses = check_degree_fullness(bases, n=2, d=3)
        degrees = [s.degree for s in statuses]
        assert degrees == [-1, 0, 1, 2, 3]

    def test_partial_basis(self):
        from lib.numeric.indexing import dim_Wn_k
        bases = {}
        k = 0
        N_k = dim_Wn_k(2, 0)  # = 4 для n=2, k=0
        b = OrthonormalBasis(N_k)
        v = np.array([[1.0, 0.0, 0.0, 0.0]])  # один вектор в 4-мерном пространстве
        b.try_add_many(v)
        bases[k] = b
        statuses = check_degree_fullness(bases, n=2, d=0)
        deg0 = next(s for s in statuses if s.degree == 0)
        assert deg0.rank == 1
        assert not deg0.is_full
