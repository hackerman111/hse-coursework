"""
Общие фикстуры для тестов алгоритма проверки.
"""

from sage.all import PolynomialRing, QQ

from lib.derivation import Derivation


class DerivationFactoryMixin:
    def setUp(self):
        self.R = PolynomialRing(QQ, "x, y")
        self.x, self.y = self.R.gens()

    def make_derivation(self, mapping):
        full_map = {self.x: 0, self.y: 0}
        full_map.update(mapping)
        return Derivation(self.R, full_map)
