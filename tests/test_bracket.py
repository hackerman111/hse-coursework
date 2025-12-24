
import unittest
from sage.all import QQ
from lib.derivation import Derivation
from tests.common import MyAlgebra

class TestBracket(unittest.TestCase):
    def setUp(self):
        self.field = QQ
        self.names = 'x,y'
        self.alg = MyAlgebra(self.field, self.names)
        self.x, self.y = self.alg.Get_gens()

    def test_bracket_basic(self):
        # D1: dx=y, dy=0
        D1 = Derivation(self.alg, {self.x: self.y, self.y: 0})
        # D2: dx=0, dy=x
        D2 = Derivation(self.alg, {self.x: 0, self.y: self.x})
        
        # [D1, D2](x) = D1(D2(x)) - D2(D1(x)) = D1(0) - D2(y) = 0 - x = -x
        bracket_D1_D2 = D1.bracket(D2)
        self.assertEqual(bracket_D1_D2(self.x), -self.x)
        
        # [D1, D2](y) = D1(D2(y)) - D2(D1(y)) = D1(x) - D2(0) = y - 0 = y
        self.assertEqual(bracket_D1_D2(self.y), self.y)

    def test_jacobi_identity(self):
        # Define 3 simple derivations
        # A: dx=1, dy=0
        # B: dx=y, dy=0
        # C: dx=0, dy=x
        A = Derivation(self.alg, {self.x: 1, self.y: 0})
        B = Derivation(self.alg, {self.x: self.y, self.y: 0})
        C = Derivation(self.alg, {self.x: 0, self.y: self.x})
        
        term1 = A.bracket(B.bracket(C))
        term2 = B.bracket(C.bracket(A))
        term3 = C.bracket(A.bracket(B))
        
        # Jacobi: [A, [B, C]] + [B, [C, A]] + [C, [A, B]] = 0
        sum_der = lambda el: term1(el) + term2(el) + term3(el)
        
        self.assertEqual(sum_der(self.x), 0)
        self.assertEqual(sum_der(self.y), 0)

    def test_bracket_properties(self):
        D1 = Derivation(self.alg, {self.x: self.y, self.y: 0})
        D2 = Derivation(self.alg, {self.x: 0, self.y: self.x})
        
        # 1. Anti-symmetry: [A, B] = -[B, A]
        bracket_12 = D1.bracket(D2)
        bracket_21 = D2.bracket(D1)
        
        self.assertEqual(bracket_12(self.x), -bracket_21(self.x))
        self.assertEqual(bracket_12(self.y), -bracket_21(self.y))
        
        # 2. Self-bracket: [A, A] = 0
        bracket_11 = D1.bracket(D1)
        self.assertEqual(bracket_11(self.x), 0)
        self.assertEqual(bracket_11(self.y), 0)
    
    def test_edge_cases_bracket(self):
        D_zero = Derivation(self.alg, {}) 
        D_any = Derivation(self.alg, {self.x: 1, self.y: 1})
        # 2. Bracket with zero is zero
        res = D_zero.bracket(D_any)
        self.assertEqual(res(self.x), 0)
