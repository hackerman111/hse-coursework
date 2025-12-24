
import unittest
from sage.all import QQ
from lib.derivation import Derivation
from tests.common import MyAlgebra

class TestDerivationBasic(unittest.TestCase):
    def setUp(self):
        self.field = QQ
        self.names = 'x,y'
        self.alg = MyAlgebra(self.field, self.names)
        self.x, self.y = self.alg.Get_gens()

    def test_derivation_call(self):
        # D(x) = y, D(y) = x
        mapping = {self.x: self.y, self.y: self.x}
        D = Derivation(self.alg, mapping)
        
        # Test on gens
        self.assertEqual(D(self.x), self.y)
        self.assertEqual(D(self.y), self.x)
        
        # Test on polynomial: D(x^2 + xy) = 2x*y + (x*x + y*y) = 2xy + x^2 + y^2
        poly = self.x**2 + self.x*self.y
        expected = 2*self.x*self.y + self.x**2 + self.y**2
        self.assertEqual(D(poly), expected)

    def test_derivation_properties(self):
        # D(x) = y, D(y) = 2x
        D = Derivation(self.alg, {self.x: self.y, self.y: 2*self.x})
        f = self.x**2 + self.y
        g = self.x * self.y
        
        # 1. Linearity: D(f + g) = D(f) + D(g)
        self.assertEqual(D(f + g), D(f) + D(g))
        
        # 2. Linearity: D(c*f) = c*D(f)
        c = 5
        self.assertEqual(D(c * f), c * D(f))
        
        # 3. Leibniz Rule: D(f*g) = f*D(g) + g*D(f)
        self.assertEqual(D(f * g), f * D(g) + g * D(f))
        
        # 4. Constants: D(1) = 0
        self.assertEqual(D(1), 0)
        self.assertEqual(D(100), 0)

    def test_edge_cases(self):
        # 1. Zero derivation
        D_zero = Derivation(self.alg, {}) # Empty mapping implies 0
        self.assertEqual(D_zero(self.x), 0)
        self.assertEqual(D_zero(self.x**10 + self.y**5), 0)
        
        # 3. Invalid input types
        D_any = Derivation(self.alg, {self.x: 1, self.y: 1})
        with self.assertRaises(TypeError):
            D_any("not a polynomial")
