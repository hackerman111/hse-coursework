"""
Тесты для lib/backends/: numeric адаптер.

Sage-backend тесты пропускаются, если SageMath не установлен.
"""

import pytest
import numpy as np

from lib.core.spec import monomial_spec, partial_spec, combine_specs
from lib.backends.numeric import to_numeric_components
from lib.generators.beldiev import beldiev_specs
from lib.numeric.recipes import make_beldiev_generators
from lib.numeric.components import decompose_generator
from lib.core.spec import to_generator


class TestNumericBackend:
    """Тесты для to_numeric_components."""

    def test_partial_spec(self):
        spec = partial_spec(2, axis=0)
        components = to_numeric_components(spec, D_max=3)
        # ∂/∂x — степень -1 (константа)
        assert -1 in components
        assert isinstance(components[-1], np.ndarray)

    def test_monomial_spec(self):
        spec = monomial_spec(2, axis=0, alpha=(1, 0))
        components = to_numeric_components(spec, D_max=3)
        # x*∂/∂x — степень 0
        assert 0 in components

    def test_matches_decompose_generator(self):
        """to_numeric_components должен давать те же результаты что decompose_generator."""
        spec = beldiev_specs(2)[1]  # V-генератор Бельдиева
        gen = to_generator(spec)

        from_backend = to_numeric_components(spec, D_max=7)
        from_legacy = decompose_generator(gen, n=2, D_max=7)

        assert set(from_backend.keys()) == set(from_legacy.keys())
        for degree in from_backend:
            np.testing.assert_allclose(from_backend[degree], from_legacy[degree])

    def test_beldiev_v_generator_n2(self):
        """V-генератор Бельдиева для n=2 имеет компоненты степеней 7."""
        v_spec = beldiev_specs(2)[1]
        components = to_numeric_components(v_spec, D_max=8)
        # V_0 = z_1^8*∂_0: степень sum((0,8))-1 = 7
        # V_1 = z_0^4*z_1^4*∂_1: степень sum((4,4))-1 = 7
        assert 7 in components

    def test_degree_truncation(self):
        """Компоненты выше D_max не попадают в результат."""
        spec = monomial_spec(2, axis=0, alpha=(5, 5))  # степень 9
        components = to_numeric_components(spec, D_max=3)
        assert 9 not in components
        assert len(components) == 0  # всё обрезано

    def test_cross_check_beldiev_n3(self):
        """Кросс-проверка для n=3."""
        specs = beldiev_specs(3)
        for i, spec in enumerate(specs):
            gen = to_generator(spec)
            from_backend = to_numeric_components(spec, D_max=12)
            from_legacy = decompose_generator(gen, n=3, D_max=12)
            assert set(from_backend.keys()) == set(from_legacy.keys())


try:
    import sage  # type: ignore  # noqa: F401
    HAS_SAGE = True
except ImportError:
    HAS_SAGE = False


@pytest.mark.skipif(not HAS_SAGE, reason="SageMath не установлен")
class TestSageBackend:
    """Тесты для Sage backend (только при наличии SageMath)."""

    def _make_ring(self, n: int):
        from sage.all import PolynomialRing, QQ  # type: ignore
        varnames = ["x", "y", "z", "u", "v", "w"][:n]
        return PolynomialRing(QQ, varnames)

    def test_to_sage_partial(self):
        from lib.backends.sage import to_sage
        ring = self._make_ring(2)
        spec = partial_spec(2, axis=0)
        derivation = to_sage(spec, ring)
        # ∂/∂x(x) = 1, ∂/∂x(y) = 0
        gens = ring.gens()
        assert derivation(gens[0]) == 1
        assert derivation(gens[1]) == 0

    def test_to_sage_monomial(self):
        from lib.backends.sage import to_sage
        ring = self._make_ring(2)
        spec = monomial_spec(2, axis=1, alpha=(1, 0))  # x * ∂/∂y
        derivation = to_sage(spec, ring)
        gens = ring.gens()
        x, y = gens
        assert derivation(x) == 0
        assert derivation(y) == x

    def test_roundtrip(self):
        from lib.backends.sage import to_sage, from_sage
        ring = self._make_ring(2)
        spec = combine_specs(
            partial_spec(2, axis=0),
            monomial_spec(2, axis=1, alpha=(1, 0)),
        )
        derivation = to_sage(spec, ring)
        recovered = from_sage(derivation)
        assert recovered == spec
