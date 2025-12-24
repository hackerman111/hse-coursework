
from sage.all import PolynomialRing
from lib.algebra import Algebra

class MyAlgebra(Algebra):
    def __init__(self, field, names):
        self.field = field
        self.gens = names
        self.n = len(names)
        self.algebra = PolynomialRing(field, names)

    def Get_gens(self):
        return self.algebra.gens()

    def Get_base_ring(self):
        return self.algebra.base_ring()
