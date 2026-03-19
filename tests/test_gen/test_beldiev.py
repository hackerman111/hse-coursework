import unittest

from sage.all import PolynomialRing

from lib.generator import check, get_Beldiev

from .fixtures import GeneratorTestMixin


class TestBeldievGenerators(GeneratorTestMixin, unittest.TestCase):
    def test_beldiev_theorem_n2(self):
        ring = PolynomialRing(self.base_field, "z1, z2")
        generators = get_Beldiev(ring)

        self.assertEqual(len(generators), 2)
        self.assertTrue(check(generators, max_iter=200))

    def test_beldiev_theorem_n3(self):
        ring = PolynomialRing(self.base_field, "z1, z2, z3")
        generators = get_Beldiev(ring)

        self.assertEqual(len(generators), 2)
        self.assertTrue(check(generators, max_iter=500))
