"""
Тесты для lib/core/spec.py: DerivationSpec и примитивные кирпичики.
"""

import pytest

from lib.core.spec import (
    DerivationSpec,
    ad,
    combine_specs,
    lie_bracket,
    monomial_spec,
    partial_spec,
    random_spec,
    spec_degree,
    to_generator,
)


class TestDerivationSpec:
    """Тесты для базовой структуры DerivationSpec."""

    def test_valid_creation(self):
        spec = DerivationSpec(n=2, terms={0: {(1, 0): 1.0}})
        assert spec.n == 2
        assert spec.terms == {0: {(1, 0): 1.0}}

    def test_empty_terms(self):
        spec = DerivationSpec(n=2, terms={})
        assert spec.terms == {}

    def test_zero_coefficient_removed(self):
        spec = DerivationSpec(n=2, terms={0: {(1, 0): 0.0, (0, 1): 1.0}})
        assert (1, 0) not in spec.terms.get(0, {})
        assert spec.terms[0] == {(0, 1): 1.0}

    def test_zero_int_coefficient_removed(self):
        spec = DerivationSpec(n=2, terms={0: {(1, 0): 0, (0, 1): 2}})
        assert (1, 0) not in spec.terms.get(0, {})
        assert spec.terms[0] == {(0, 1): 2}

    def test_axis_cleared_if_all_zero(self):
        spec = DerivationSpec(n=2, terms={0: {(1, 0): 0.0}})
        assert 0 not in spec.terms

    def test_invalid_target_var(self):
        with pytest.raises(ValueError, match="target_var"):
            DerivationSpec(n=2, terms={2: {(0, 0): 1.0}})

    def test_negative_target_var(self):
        with pytest.raises(ValueError, match="target_var"):
            DerivationSpec(n=2, terms={-1: {(0, 0): 1.0}})

    def test_wrong_alpha_length(self):
        with pytest.raises(ValueError, match="len\\(alpha\\)"):
            DerivationSpec(n=2, terms={0: {(1, 0, 0): 1.0}})

    def test_invalid_n(self):
        with pytest.raises(ValueError, match="положительным"):
            DerivationSpec(n=0, terms={})

    def test_equality(self):
        s1 = DerivationSpec(n=2, terms={0: {(1, 0): 1.0}})
        s2 = DerivationSpec(n=2, terms={0: {(1, 0): 1.0}})
        assert s1 == s2

    def test_inequality_different_n(self):
        s1 = DerivationSpec(n=2, terms={0: {(1, 0): 1.0}})
        s2 = DerivationSpec(n=3, terms={0: {(1, 0, 0): 1.0}})
        assert s1 != s2

    def test_inequality_different_coeffs(self):
        s1 = DerivationSpec(n=2, terms={0: {(1, 0): 1.0}})
        s2 = DerivationSpec(n=2, terms={0: {(1, 0): 2.0}})
        assert s1 != s2

    def test_frozen(self):
        spec = DerivationSpec(n=2, terms={0: {(1, 0): 1.0}})
        with pytest.raises((AttributeError, TypeError)):
            spec.n = 3  # type: ignore


class TestMonomialSpec:
    """Тесты для monomial_spec."""

    def test_basic(self):
        spec = monomial_spec(2, axis=0, alpha=(2, 1), coeff=1.0)
        assert spec.n == 2
        assert spec.terms == {0: {(2, 1): 1.0}}

    def test_with_coeff(self):
        spec = monomial_spec(2, axis=1, alpha=(0, 3), coeff=-2.0)
        assert spec.terms == {1: {(0, 3): -2.0}}

    def test_int_coeff(self):
        spec = monomial_spec(2, axis=0, alpha=(1, 0), coeff=3)
        assert spec.terms == {0: {(1, 0): 3}}

    def test_wrong_alpha_length(self):
        with pytest.raises(ValueError, match="len\\(alpha\\)"):
            monomial_spec(2, axis=0, alpha=(1, 0, 0))

    def test_n3(self):
        spec = monomial_spec(3, axis=2, alpha=(1, 1, 1), coeff=2.0)
        assert spec.terms == {2: {(1, 1, 1): 2.0}}


class TestPartialSpec:
    """Тесты для partial_spec."""

    def test_partial_x(self):
        spec = partial_spec(2, axis=0)
        assert spec.terms == {0: {(0, 0): 1}}

    def test_partial_y(self):
        spec = partial_spec(2, axis=1)
        assert spec.terms == {1: {(0, 0): 1}}

    def test_equivalent_to_monomial_spec(self):
        ps = partial_spec(2, axis=0)
        ms = monomial_spec(2, axis=0, alpha=(0, 0), coeff=1)
        assert ps == ms

    def test_n3(self):
        spec = partial_spec(3, axis=2)
        assert spec.terms == {2: {(0, 0, 0): 1}}


class TestCombineSpecs:
    """Тесты для combine_specs."""

    def test_two_specs(self):
        u = partial_spec(2, axis=0)
        v = partial_spec(2, axis=1)
        combined = combine_specs(u, v)
        assert combined.terms == {0: {(0, 0): 1}, 1: {(0, 0): 1}}

    def test_cancellation(self):
        s1 = monomial_spec(2, axis=0, alpha=(1, 0), coeff=1.0)
        s2 = monomial_spec(2, axis=0, alpha=(1, 0), coeff=-1.0)
        combined = combine_specs(s1, s2)
        assert 0 not in combined.terms  # занулился

    def test_accumulation(self):
        s1 = monomial_spec(2, axis=0, alpha=(1, 0), coeff=2.0)
        s2 = monomial_spec(2, axis=0, alpha=(1, 0), coeff=3.0)
        combined = combine_specs(s1, s2)
        assert combined.terms[0][(1, 0)] == pytest.approx(5.0)

    def test_single_spec(self):
        u = partial_spec(2, axis=0)
        combined = combine_specs(u)
        assert combined == u

    def test_mismatched_n(self):
        s1 = partial_spec(2, axis=0)
        s2 = partial_spec(3, axis=0)
        with pytest.raises(ValueError, match="одинаковое n"):
            combine_specs(s1, s2)

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="хотя бы одного"):
            combine_specs()


class TestSpecDegree:
    """Тесты для spec_degree."""

    def test_partial_derivative(self):
        # ∂/∂x имеет степень sum((0,0))-1 = -1
        spec = partial_spec(2, axis=0)
        assert spec_degree(spec) == -1

    def test_monomial_degree_1(self):
        # x*∂/∂x: sum((1,0))-1 = 0 → степень 0
        spec = monomial_spec(2, axis=0, alpha=(1, 0))
        assert spec_degree(spec) == 0

    def test_monomial_degree_2(self):
        # x^2*y*∂/∂x: sum((2,1))-1 = 2 → степень 2
        spec = monomial_spec(2, axis=0, alpha=(2, 1))
        assert spec_degree(spec) == 2

    def test_max_over_components(self):
        u = monomial_spec(2, axis=0, alpha=(3, 0))  # степень 2
        v = monomial_spec(2, axis=1, alpha=(1, 1))  # степень 1
        combined = combine_specs(u, v)
        assert spec_degree(combined) == 2

    def test_empty_spec(self):
        spec = DerivationSpec(n=2)
        assert spec_degree(spec) == -1


class TestRandomSpec:
    """Тесты для random_spec."""

    def test_returns_valid_spec(self):
        spec = random_spec(n=2, max_degree=3, seed=42)
        assert spec.n == 2
        for axis, mono_dict in spec.terms.items():
            assert 0 <= axis < 2
            for alpha in mono_dict:
                assert len(alpha) == 2

    def test_degree_bounded(self):
        max_d = 3
        spec = random_spec(n=2, max_degree=max_d, seed=0)
        for mono_dict in spec.terms.values():
            for alpha in mono_dict:
                assert sum(alpha) - 1 <= max_d

    def test_reproducible_with_seed(self):
        s1 = random_spec(n=2, max_degree=3, seed=7)
        s2 = random_spec(n=2, max_degree=3, seed=7)
        assert s1 == s2

    def test_different_seeds(self):
        s1 = random_spec(n=2, max_degree=3, seed=1)
        s2 = random_spec(n=2, max_degree=3, seed=2)
        # С высокой вероятностью должны различаться
        assert s1 != s2

    def test_not_empty(self):
        spec = random_spec(n=2, max_degree=2, seed=100)
        assert len(spec.terms) > 0

    def test_n3(self):
        spec = random_spec(n=3, max_degree=2, seed=42)
        assert spec.n == 3
        for alpha_dict in spec.terms.values():
            for alpha in alpha_dict:
                assert len(alpha) == 3


class TestToGenerator:
    """Тесты для to_generator."""

    def test_basic(self):
        spec = monomial_spec(2, axis=0, alpha=(2, 1), coeff=3)
        gen = to_generator(spec)
        assert gen == {0: {(2, 1): 3.0}}

    def test_coefficients_are_float(self):
        spec = monomial_spec(2, axis=0, alpha=(1, 0), coeff=2)
        gen = to_generator(spec)
        assert isinstance(gen[0][(1, 0)], float)

    def test_multiple_components(self):
        u = partial_spec(2, axis=0)
        v = monomial_spec(2, axis=1, alpha=(1, 0), coeff=1.0)
        combined = combine_specs(u, v)
        gen = to_generator(combined)
        assert 0 in gen
        assert 1 in gen
        assert gen[0] == {(0, 0): 1.0}
        assert gen[1] == {(1, 0): 1.0}

    def test_empty_spec(self):
        spec = DerivationSpec(n=2)
        gen = to_generator(spec)
        assert gen == {}

    def test_compatibility_with_numeric_layer(self):
        """to_generator должен возвращать структуру совместимую с Generator type."""
        spec = monomial_spec(2, axis=0, alpha=(1, 1), coeff=2.5)
        gen = to_generator(spec)
        # Generator = Dict[int, Dict[Tuple[int,...], float]]
        for axis, mono_dict in gen.items():
            assert isinstance(axis, int)
            for alpha, coeff in mono_dict.items():
                assert isinstance(alpha, tuple)
                assert isinstance(coeff, float)


class TestLieBracket:
    """Тесты для lie_bracket и ad."""

    def test_partial_bracket_x_term(self):
        # [∂/∂x, x*∂/∂y] = ∂/∂y
        D1 = partial_spec(2, 0)
        D2 = monomial_spec(2, 1, (1, 0))
        result = lie_bracket(D1, D2)
        expected = partial_spec(2, 1)
        assert result == expected

    def test_anticommutativity(self):
        # [D1, D2] = -[D2, D1]
        D1 = monomial_spec(2, 0, (1, 0))  # x*∂/∂x
        D2 = monomial_spec(2, 1, (0, 1))  # y*∂/∂y
        b12 = lie_bracket(D1, D2)
        b21 = lie_bracket(D2, D1)
        # b12 + b21 == 0 → combine_specs should give zero spec
        combined = combine_specs(b12, b21)
        assert combined.terms == {}

    def test_self_bracket_is_zero(self):
        # [D, D] = 0
        D = monomial_spec(2, 0, (2, 1))
        result = lie_bracket(D, D)
        assert result.terms == {}

    def test_bracket_of_partials_is_zero(self):
        # [∂/∂x, ∂/∂y] = 0 (частные производные коммутируют)
        D1 = partial_spec(2, 0)
        D2 = partial_spec(2, 1)
        result = lie_bracket(D1, D2)
        assert result.terms == {}

    def test_euler_formula(self):
        # [x*∂/∂x + y*∂/∂y, x^2*∂/∂x] = x^2*∂/∂x (Euler: [E, D_k] = k*D_k для D_k однородного степени k)
        E = combine_specs(monomial_spec(2, 0, (1, 0)), monomial_spec(2, 1, (0, 1)))
        D = monomial_spec(2, 0, (2, 0))  # x^2*∂/∂x, степень 1
        result = lie_bracket(E, D)
        assert result == D

    def test_d_max_truncation(self):
        D1 = partial_spec(2, 0)
        D2 = monomial_spec(2, 1, (3, 0))  # x^3*∂/∂y, результат [D1, D2] = 3x^2*∂/∂y (степень 1)
        full = lie_bracket(D1, D2)
        trunc = lie_bracket(D1, D2, D_max=0)
        # Полный результат — степень 1, при D_max=0 должен быть обрезан
        assert full.terms != {}
        assert trunc.terms == {}

    def test_n_mismatch_raises(self):
        D1 = partial_spec(2, 0)
        D2 = partial_spec(3, 0)
        with pytest.raises(ValueError, match="n"):
            lie_bracket(D1, D2)

    def test_n3_basic(self):
        # [∂/∂x, x*∂/∂z] = ∂/∂z
        D1 = partial_spec(3, 0)
        D2 = monomial_spec(3, 2, (1, 0, 0))
        result = lie_bracket(D1, D2)
        assert result == partial_spec(3, 2)

    def test_ad_equals_lie_bracket(self):
        D1 = partial_spec(2, 0)
        D2 = monomial_spec(2, 1, (1, 0))
        assert ad(D1)(D2) == lie_bracket(D1, D2)

    def test_ad_with_d_max(self):
        D1 = partial_spec(2, 0)
        D2 = monomial_spec(2, 1, (3, 0))
        assert ad(D1)(D2, D_max=0) == lie_bracket(D1, D2, D_max=0)

    def test_jacobi_identity(self):
        # [D1,[D2,D3]] + [D2,[D3,D1]] + [D3,[D1,D2]] = 0
        D1 = partial_spec(2, 0)
        D2 = monomial_spec(2, 1, (1, 0))
        D3 = monomial_spec(2, 0, (0, 1))
        t1 = lie_bracket(D1, lie_bracket(D2, D3))
        t2 = lie_bracket(D2, lie_bracket(D3, D1))
        t3 = lie_bracket(D3, lie_bracket(D1, D2))
        result = combine_specs(t1, t2, t3) if any(
            s.terms for s in [t1, t2, t3]
        ) else DerivationSpec(n=2)
        assert result.terms == {}
