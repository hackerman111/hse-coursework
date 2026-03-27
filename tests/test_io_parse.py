"""
Тесты для lib/io/parse.py: разбор и форматирование генераторов.
"""

import pytest

from lib.core.spec import DerivationSpec, monomial_spec, partial_spec, combine_specs
from lib.io.parse import (
    parse_generator_spec,
    spec_from_string,
    spec_to_string,
    describe_generator,
    generator_term_count,
    generator_max_degree,
)


class TestParseGeneratorSpec:
    """Тесты для parse_generator_spec."""

    def test_simple_partial(self):
        gen = parse_generator_spec("dx", 2)
        assert gen == {0: {(0, 0): 1.0}}

    def test_partial_dy(self):
        gen = parse_generator_spec("dy", 2)
        assert gen == {1: {(0, 0): 1.0}}

    def test_monomial(self):
        gen = parse_generator_spec("x*dy", 2)
        assert gen == {1: {(1, 0): 1.0}}

    def test_monomial_with_power(self):
        gen = parse_generator_spec("x^2*dy", 2)
        assert gen == {1: {(2, 0): 1.0}}

    def test_sum(self):
        gen = parse_generator_spec("dx + x*dy", 2)
        assert gen == {0: {(0, 0): 1.0}, 1: {(1, 0): 1.0}}

    def test_difference(self):
        gen = parse_generator_spec("dx - dy", 2)
        assert gen[0][(0, 0)] == pytest.approx(1.0)
        assert gen[1][(0, 0)] == pytest.approx(-1.0)

    def test_coefficient(self):
        gen = parse_generator_spec("3*x*dx", 2)
        assert gen[0][(1, 0)] == pytest.approx(3.0)

    def test_cancellation(self):
        # При полном занулении коэффициентов — ValueError (оригинальное поведение)
        with pytest.raises(ValueError, match="занулился"):
            parse_generator_spec("dx - dx", 2)

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="пустая"):
            parse_generator_spec("", 2)

    def test_no_derivative_raises(self):
        with pytest.raises(ValueError, match="не указана производная"):
            parse_generator_spec("x*y", 2)

    def test_complex_expression(self):
        gen = parse_generator_spec("x^5*y^5*dx", 2)
        assert gen == {0: {(5, 5): 1.0}}


class TestSpecFromString:
    """Тесты для spec_from_string."""

    def test_basic(self):
        spec = spec_from_string("x^5*y^5*dx", 2)
        assert isinstance(spec, DerivationSpec)
        assert spec.n == 2
        assert spec.terms == {0: {(5, 5): 1.0}}

    def test_partial(self):
        spec = spec_from_string("dx", 2)
        assert spec == partial_spec(2, axis=0)

    def test_sum(self):
        spec = spec_from_string("dx + x*dy", 2)
        expected = combine_specs(partial_spec(2, 0), monomial_spec(2, 1, (1, 0)))
        assert spec == expected


class TestSpecToString:
    """Тесты для spec_to_string."""

    def test_partial_dx(self):
        spec = partial_spec(2, axis=0)
        s = spec_to_string(spec)
        assert "dx" in s

    def test_partial_dy(self):
        spec = partial_spec(2, axis=1)
        s = spec_to_string(spec)
        assert "dy" in s

    def test_monomial(self):
        spec = monomial_spec(2, axis=0, alpha=(2, 1))
        s = spec_to_string(spec)
        # Должен содержать x^2, y, dx
        assert "x^2" in s
        assert "y" in s
        assert "dx" in s

    def test_zero_spec(self):
        spec = DerivationSpec(n=2)
        assert spec_to_string(spec) == "0"

    def test_custom_variables(self):
        spec = partial_spec(2, axis=0)
        s = spec_to_string(spec, variables=["a", "b"])
        assert "da" in s

    def test_roundtrip_partial(self):
        """spec_from_string(spec_to_string(spec)) == spec."""
        spec = partial_spec(2, axis=0)
        s = spec_to_string(spec)
        recovered = spec_from_string(s, 2)
        assert recovered == spec

    def test_roundtrip_monomial(self):
        spec = monomial_spec(2, axis=1, alpha=(3, 2), coeff=1)
        s = spec_to_string(spec)
        recovered = spec_from_string(s, 2)
        assert recovered == spec

    def test_roundtrip_sum(self):
        spec = combine_specs(partial_spec(2, 0), monomial_spec(2, 1, (1, 0)))
        s = spec_to_string(spec)
        recovered = spec_from_string(s, 2)
        assert recovered == spec

    def test_roundtrip_negative_coeff(self):
        spec = monomial_spec(2, axis=0, alpha=(1, 1), coeff=-2)
        s = spec_to_string(spec)
        recovered = spec_from_string(s, 2)
        assert recovered == spec


class TestDescribeGenerator:
    """Тесты для describe_generator."""

    def test_partial(self):
        gen = {0: {(0, 0): 1.0}}
        s = describe_generator(gen, 2)
        assert "∂" in s or "d" in s.lower()

    def test_max_terms(self):
        gen = {0: {(1, 0): 1.0, (0, 1): 1.0, (2, 0): 1.0}}
        s = describe_generator(gen, 2, max_terms=2)
        assert "..." in s


class TestGeneratorStats:
    """Тесты для generator_term_count и generator_max_degree."""

    def test_term_count(self):
        gen = {0: {(1, 0): 1.0, (0, 1): 1.0}, 1: {(2, 0): 1.0}}
        assert generator_term_count(gen) == 3

    def test_max_degree(self):
        gen = {0: {(3, 0): 1.0}, 1: {(1, 1): 1.0}}
        # sum((3,0))-1 = 2, sum((1,1))-1 = 1 → max = 2
        assert generator_max_degree(gen) == 2

    def test_partial_degree(self):
        gen = {0: {(0, 0): 1.0}}
        # sum((0,0))-1 = -1
        assert generator_max_degree(gen) == -1


class TestBackwardsCompatibility:
    """Проверка обратной совместимости через lib.numeric.specs."""

    def test_parse_still_works_from_numeric(self):
        from lib.numeric.specs import parse_generator_spec as old_parse
        gen = old_parse("dx + x*dy", 2)
        assert gen == {0: {(0, 0): 1.0}, 1: {(1, 0): 1.0}}

    def test_describe_still_works_from_numeric(self):
        from lib.numeric.specs import describe_generator as old_describe
        gen = {0: {(0, 0): 1.0}}
        s = old_describe(gen, 2)
        assert isinstance(s, str)
