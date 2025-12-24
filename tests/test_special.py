
import unittest
from sage.all import QQ, Matrix
from lib.derivation import Derivation, LinearDerivation, WeitzenbockDerivation
from tests.common import MyAlgebra

class TestSpecialDerivations(unittest.TestCase):
    def setUp(self):
        self.field = QQ
        self.names = 'x,y,z'
        self.alg = MyAlgebra(self.field, self.names)
        self.x, self.y, self.z = self.alg.Get_gens()

    def test_triangle_derivation(self):
        # Triangular: D(x) in K[y,z], D(y) in K[z], D(z) in K
        # D(x) = y + z^2
        # D(y) = 3*z + 1
        # D(z) = 5
        D_tri = Derivation(self.alg, {self.x: self.y + self.z**2, self.y: 3*self.z + 1, self.z: 5})
        self.assertTrue(D_tri.is_triangular)
        
        # Not triangular: D(y) depends on x
        D_bad = Derivation(self.alg, {self.x: self.z, self.y: self.x, self.z: 0})
        self.assertFalse(D_bad.is_triangular)
        
        # Not triangular: D(z) depends on z (must be constant)
        D_bad2 = Derivation(self.alg, {self.x: self.y, self.y: self.z, self.z: self.z})
        self.assertFalse(D_bad2.is_triangular)

    def test_jacobian_derivation(self):
        # n=3 (x,y,z). n-1=2 polynomials.
        # Example from user: f1 = x^2 + y^2, f2 = z
        # Expected D = 2y*dx - 2x*dy + 0*dz
        f1 = self.x**2 + self.y**2
        f2 = self.z
        
        D = Derivation.from_jacobian(self.alg, [f1, f2])
        
        self.assertEqual(D(self.x), 2*self.y)
        self.assertEqual(D(self.y), -2*self.x)
        self.assertEqual(D(self.z), 0)

    def test_linear_derivation(self):
        # Need 2 vars for this test usually, but we have 3.
        # Let's reduce algebra context or just use 3x3 matrix for simplicity, 
        # OR just use top-left 2x2 block logic extended to 3.
        
        # Actually proper LinearDerivation requires matrix size == n.
        # So I'll use 3x3 matrix.
        # Matrix A = [[0, 1, 0], [1, 0, 0], [0, 0, 0]]
        # D(x) = y
        # D(y) = x
        # D(z) = 0
        matrix = Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
        D = LinearDerivation(self.alg, matrix)
        
        self.assertEqual(D(self.x), self.y)
        self.assertEqual(D(self.y), self.x)
        self.assertEqual(D(self.z), 0)
        
        # Check linearity
        # D(2x + y) = 2y + x
        poly = 2*self.x + self.y
        self.assertEqual(D(poly), 2*self.y + self.x)

    def test_weitzenbock_derivation(self):
        # Nilpotent matrix (ordering x, y, z)
        # [[0, 1, 0], 
        #  [0, 0, 1],
        #  [0, 0, 0]]
        # D(x) = y
        # D(y) = z
        # D(z) = 0
        matrix = Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
        D = WeitzenbockDerivation(self.alg, matrix)
        
        self.assertEqual(D(self.x), self.y)
        self.assertEqual(D(self.y), self.z)
        self.assertEqual(D(self.z), 0)
        
        # Should be triangular as well
        self.assertTrue(D.is_triangular)
        
        # Failing case: non-nilpotent matrix
        non_nil = Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]])
        with self.assertRaises(ValueError):
            WeitzenbockDerivation(self.alg, non_nil)
