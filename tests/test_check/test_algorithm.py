import unittest

from lib.generator import check

from .fixtures import DerivationFactoryMixin


class TestCheckAlgorithm(DerivationFactoryMixin, unittest.TestCase):
    """
    Тестирование самой функции проверки (check).
    """

    def test_check_with_full_basis(self):
        dx = self.make_derivation({self.x: 1})
        dy = self.make_derivation({self.y: 1})

        self.assertTrue(check([dx, dy]))

    def test_check_with_incomplete_basis(self):
        dx = self.make_derivation({self.x: 1})

        self.assertFalse(check([dx]))

    def test_check_with_commutator_generation(self):
        dx = self.make_derivation({self.x: 1})
        x_dy = self.make_derivation({self.y: self.x})

        self.assertTrue(check([dx, x_dy]))

    def test_check_with_subalgebra(self):
        x_dx = self.make_derivation({self.x: self.x})
        y_dy = self.make_derivation({self.y: self.y})

        self.assertFalse(check([x_dx, y_dy]))
