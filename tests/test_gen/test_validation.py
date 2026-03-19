import unittest

from sage.all import PolynomialRing

from lib.generator import get_Andristy, get_Beldiev

from .fixtures import GeneratorTestMixin


class TestGeneratorValidation(GeneratorTestMixin, unittest.TestCase):
    def test_argument_validation(self):
        ring = PolynomialRing(self.base_field, "z1")

        with self.assertRaises(ValueError):
            get_Beldiev(ring)

        with self.assertRaises(ValueError):
            get_Andristy(ring)
