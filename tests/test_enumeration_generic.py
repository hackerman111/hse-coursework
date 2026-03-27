"""
Тесты для lib/solver/enumeration.py: обобщённый BasisTerm для произвольного n.
"""

import pytest

from lib.solver.enumeration import (
    BasisTerm,
    build_basis_terms,
    build_generator,
    candidate_to_spec,
    estimate_candidate_count,
    format_coefficient,
    iter_candidates,
    parse_coefficients,
)


class TestBasisTermGeneric:
    """Обобщённый BasisTerm для произвольного n."""

    def test_n2_basic(self):
        term = BasisTerm(alpha=(2, 1), target_var=0, n=2)
        assert term.alpha == (2, 1)
        assert term.target_var == 0
        assert term.n == 2

    def test_n3_basic(self):
        term = BasisTerm(alpha=(1, 0, 2), target_var=2, n=3)
        assert term.n == 3

    def test_alpha_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="длину n="):
            BasisTerm(alpha=(1, 0), target_var=0, n=3)

    def test_target_var_out_of_range_raises(self):
        with pytest.raises(ValueError, match="target_var"):
            BasisTerm(alpha=(1, 0), target_var=2, n=2)

    def test_factor_spec_n2(self):
        term = BasisTerm(alpha=(2, 1), target_var=0, n=2)
        s = term.factor_spec()
        assert "x^2" in s
        assert "y" in s
        assert "dx" in s

    def test_factor_spec_n3(self):
        term = BasisTerm(alpha=(1, 0, 2), target_var=2, n=3)
        s = term.factor_spec()
        assert "x" in s
        assert "z^2" in s
        assert "dz" in s

    def test_factor_spec_constant(self):
        term = BasisTerm(alpha=(0, 0), target_var=1, n=2)
        s = term.factor_spec()
        assert "dy" in s
        assert "x" not in s

    def test_factor_spec_custom_variables(self):
        term = BasisTerm(alpha=(1, 0), target_var=0, n=2)
        s = term.factor_spec(variables=["a", "b"])
        assert "a" in s
        assert "da" in s

    def test_hashable_and_frozen(self):
        t1 = BasisTerm(alpha=(1, 0), target_var=0, n=2)
        t2 = BasisTerm(alpha=(1, 0), target_var=0, n=2)
        assert t1 == t2
        assert hash(t1) == hash(t2)
        s = {t1, t2}
        assert len(s) == 1


class TestBuildBasisTermsGeneric:
    """build_basis_terms для произвольного n."""

    def test_n2_max_degree_0_no_constants(self):
        # max_degree=0, нет константных членов → нет термов
        terms = build_basis_terms(n=2, max_degree=0, include_constants=False)
        assert len(terms) == 0

    def test_n2_max_degree_0_with_constants(self):
        # max_degree=0, с константными членами → 2 варианта: dx, dy
        terms = build_basis_terms(n=2, max_degree=0, include_constants=True)
        assert len(terms) == 2

    def test_n2_max_degree_1_with_constants(self):
        # Мономы: (0,0), (1,0), (0,1) — 3 варианта * 2 направления = 6
        terms = build_basis_terms(n=2, max_degree=1, include_constants=True)
        assert len(terms) == 6

    def test_n3_max_degree_1_with_constants(self):
        # Мономы для n=3, degree≤1: (0,0,0),(1,0,0),(0,1,0),(0,0,1) — 4 * 3 = 12
        terms = build_basis_terms(n=3, max_degree=1, include_constants=True)
        assert len(terms) == 12

    def test_all_terms_have_correct_n(self):
        terms = build_basis_terms(n=3, max_degree=2, include_constants=True)
        for t in terms:
            assert t.n == 3
            assert len(t.alpha) == 3

    def test_target_var_range(self):
        terms = build_basis_terms(n=3, max_degree=1, include_constants=True)
        target_vars = {t.target_var for t in terms}
        assert target_vars == {0, 1, 2}


class TestBuildGeneratorGeneric:
    """build_generator для обобщённых термов."""

    def test_single_term(self):
        term = BasisTerm(alpha=(2, 1), target_var=0, n=2)
        gen = build_generator([(term, 3.0)])
        assert gen == {0: {(2, 1): 3.0}}

    def test_two_terms_different_vars(self):
        t1 = BasisTerm(alpha=(1, 0), target_var=0, n=2)
        t2 = BasisTerm(alpha=(0, 1), target_var=1, n=2)
        gen = build_generator([(t1, 1.0), (t2, -2.0)])
        assert gen[0][(1, 0)] == pytest.approx(1.0)
        assert gen[1][(0, 1)] == pytest.approx(-2.0)

    def test_two_terms_same_alpha_same_var_sum(self):
        t1 = BasisTerm(alpha=(1, 0), target_var=0, n=2)
        gen = build_generator([(t1, 1.0), (t1, 2.0)])
        assert gen[0][(1, 0)] == pytest.approx(3.0)

    def test_cancellation_removes_entry(self):
        t1 = BasisTerm(alpha=(1, 0), target_var=0, n=2)
        gen = build_generator([(t1, 1.0), (t1, -1.0)])
        assert gen == {} or 0 not in gen or (1, 0) not in gen.get(0, {})

    def test_n3_term(self):
        t = BasisTerm(alpha=(1, 0, 0), target_var=2, n=3)
        gen = build_generator([(t, 5.0)])
        assert gen == {2: {(1, 0, 0): 5.0}}


class TestIterCandidatesGeneric:
    """iter_candidates для обобщённого API."""

    def test_single_term_single_coeff(self):
        terms = [BasisTerm(alpha=(1, 0), target_var=0, n=2)]
        candidates = list(iter_candidates(terms, coeffs=[1.0], min_terms=1, max_terms=1))
        assert len(candidates) == 1
        assert candidates[0][0][0] == terms[0]
        assert candidates[0][0][1] == 1.0

    def test_two_terms_two_coeffs(self):
        terms = build_basis_terms(n=2, max_degree=1, include_constants=True)
        candidates = list(iter_candidates(terms[:3], coeffs=[1.0, -1.0], min_terms=1, max_terms=1))
        assert len(candidates) == 6  # 3 terms * 2 coeffs


class TestCandidateToSpecGeneric:
    """candidate_to_spec для обобщённого API."""

    def test_simple(self):
        t = BasisTerm(alpha=(2, 1), target_var=0, n=2)
        s = candidate_to_spec([(t, 1.0)])
        assert "x^2" in s
        assert "dx" in s

    def test_coefficient_shown(self):
        t = BasisTerm(alpha=(1, 0), target_var=0, n=2)
        s = candidate_to_spec([(t, 3.0)])
        assert "3" in s

    def test_custom_variables(self):
        t = BasisTerm(alpha=(1, 0), target_var=0, n=2)
        s = candidate_to_spec([(t, 1.0)], variables=["a", "b"])
        assert "da" in s


class TestEstimateAndParseCoefficientSharedFromOldModule:
    """Утилиты реэкспортируются и работают корректно."""

    def test_estimate_matches_small_example(self):
        assert estimate_candidate_count(6, 2, 1, 2) == 72

    def test_parse_coefficients_deduplicates(self):
        coeffs = parse_coefficients("1,-1,0,1")
        assert coeffs == [1.0, -1.0]

    def test_format_coefficient_integer(self):
        assert format_coefficient(3.0) == "3"

    def test_format_coefficient_float(self):
        s = format_coefficient(1.5)
        assert "1.5" in s
