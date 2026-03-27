"""
Тесты для lib/generators/beldiev.py и lib/generators/andrist.py.

Проверяет:
- структуру DerivationSpec для каждого рецепта
- кросс-совместимость с numeric layer (make_beldiev_generators)
"""

import pytest

from lib.core.spec import partial_spec, spec_degree, to_generator
from lib.generators.beldiev import beldiev_specs
from lib.generators.andrist import andrist_specs
from lib.numeric.recipes import make_beldiev_generators


class TestBeldievSpecs:
    """Тесты для beldiev_specs."""

    def test_returns_two_specs(self):
        specs = beldiev_specs(2)
        assert len(specs) == 2

    def test_n2_u_is_partial_y(self):
        u, v = beldiev_specs(2)
        # U = ∂/∂z_1 (в W_2 это ∂/∂y)
        assert u == partial_spec(2, axis=1)

    def test_n2_v_structure(self):
        u, v = beldiev_specs(2)
        # V должен иметь компоненты для ∂/∂x и ∂/∂y
        assert 0 in v.terms  # компонента ∂/∂x
        assert 1 in v.terms  # компонента ∂/∂y

    def test_n2_v_dx_component(self):
        u, v = beldiev_specs(2)
        # V_0 = z_1^8 * ∂/∂z_0: power = 4*(0+2) = 8, alpha = (0, 8)
        assert v.terms[0] == {(0, 8): 1}

    def test_n2_v_dy_component(self):
        u, v = beldiev_specs(2)
        # V_1 = z_0^4 * z_1^4 * ∂/∂z_1
        assert v.terms[1] == {(4, 4): 1}

    def test_n3_structure(self):
        specs = beldiev_specs(3)
        u, v = specs
        assert u == partial_spec(3, axis=2)
        # V должен иметь 3 компоненты
        assert set(v.terms.keys()) == {0, 1, 2}

    def test_n3_v_components(self):
        u, v = beldiev_specs(3)
        # k=0: power=4*(0+2)=8, alpha=(0,8,8)
        assert v.terms[0] == {(0, 8, 8): 1}
        # k=1: power=4*(1+2)=12, alpha=(0,0,12)
        assert v.terms[1] == {(0, 0, 12): 1}
        # k=2: alpha=(4,4,4)
        assert v.terms[2] == {(4, 4, 4): 1}

    def test_requires_n_ge_2(self):
        with pytest.raises(ValueError, match="n >= 2"):
            beldiev_specs(1)

    def test_degree_of_v(self):
        u, v = beldiev_specs(2)
        # V_0 имеет степень sum((0,8))-1 = 7, V_1: sum((4,4))-1 = 7 → max = 7
        assert spec_degree(v) == 7

    def test_cross_check_n2_with_numeric(self):
        """Кросс-проверка: beldiev_specs совпадает с make_beldiev_generators."""
        specs = beldiev_specs(2)
        from_specs = [to_generator(s) for s in specs]
        from_numeric = make_beldiev_generators(2)

        assert len(from_specs) == len(from_numeric)
        for spec_gen, num_gen in zip(from_specs, from_numeric):
            assert set(spec_gen.keys()) == set(num_gen.keys())
            for axis in spec_gen:
                assert set(spec_gen[axis].keys()) == set(num_gen[axis].keys())
                for alpha in spec_gen[axis]:
                    assert abs(spec_gen[axis][alpha] - num_gen[axis][alpha]) < 1e-12

    def test_cross_check_n3_with_numeric(self):
        """Кросс-проверка для n=3."""
        specs = beldiev_specs(3)
        from_specs = [to_generator(s) for s in specs]
        from_numeric = make_beldiev_generators(3)

        for spec_gen, num_gen in zip(from_specs, from_numeric):
            assert set(spec_gen.keys()) == set(num_gen.keys())
            for axis in spec_gen:
                assert set(spec_gen[axis].keys()) == set(num_gen[axis].keys())

    def test_make_beldiev_generators_unchanged(self):
        """make_beldiev_generators должен работать как раньше."""
        result = make_beldiev_generators(2)
        assert len(result) == 2
        # U = ∂/∂y
        assert result[0] == {1: {(0, 0): 1.0}}
        # V_0 = z_1^8
        assert result[1][0] == {(0, 8): 1.0}
        # V_1 = z_0^4 z_1^4
        assert result[1][1] == {(4, 4): 1.0}


class TestAndristSpecs:
    """Тесты для andrist_specs."""

    def test_returns_three_specs(self):
        specs = andrist_specs(2)
        assert len(specs) == 3

    def test_n2_u_is_partial_y(self):
        u, v, w = andrist_specs(2)
        assert u == partial_spec(2, axis=1)

    def test_n2_v_structure(self):
        u, v, w = andrist_specs(2)
        # V для n=2: ∂_1 компонента = константа 1, ∂_0 компонента = z_1^3
        assert 0 in v.terms  # компонента ∂/∂x: z_y^3 ∂_x
        assert 1 in v.terms  # компонента ∂/∂y: константа

    def test_n2_v_components(self):
        u, v, w = andrist_specs(2)
        # k=0: alpha[1]=3, j in range(2,2)=empty → alpha=(0,3)
        assert v.terms[0] == {(0, 3): 1}
        # Константная компонента: alpha=(0,0)
        assert v.terms[1] == {(0, 0): 1}

    def test_n2_w_structure(self):
        u, v, w = andrist_specs(2)
        # W = z_0^2 * z_1 * ∂/∂z_1
        assert w.terms == {1: {(2, 1): 1}}

    def test_n3_structure(self):
        specs = andrist_specs(3)
        assert len(specs) == 3
        u, v, w = specs
        assert u == partial_spec(3, axis=2)

    def test_n3_v_components(self):
        u, v, w = andrist_specs(3)
        # Константная: target=2, alpha=(0,0,0)
        assert v.terms[2] == {(0, 0, 0): 1}
        # k=0: alpha[1]=3, alpha[2]=1 → (0,3,1)
        assert v.terms[0] == {(0, 3, 1): 1}
        # k=1: alpha[2]=3 → (0,0,3)
        assert v.terms[1] == {(0, 0, 3): 1}

    def test_n3_w_component(self):
        u, v, w = andrist_specs(3)
        # W = z_0^2 * z_1^2 * z_2 * ∂/∂z_2
        assert w.terms == {2: {(2, 2, 1): 1}}

    def test_requires_n_ge_2(self):
        with pytest.raises(ValueError, match="n >= 2"):
            andrist_specs(1)

    def test_n4_structure(self):
        specs = andrist_specs(4)
        assert len(specs) == 3
        u, v, w = specs
        assert u == partial_spec(4, axis=3)
        assert set(v.terms.keys()) == {0, 1, 2, 3}
