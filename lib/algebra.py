"""
Base Algebra Class
"""


from abc import ABC, abstractmethod
from sage.all import PolynomialRing

class Algebra(ABC):
    """
    Abstract base class for algebras suitable for differentiation.
    """

    @abstractmethod
    def __init__(self, field, names):
        self.field = field
        self.gens = names
        self.n = len(names)
        self.algebra = PolynomialRing(field, names)

    @abstractmethod
    def Get_gens(self):
        """
        Возвращает генераторы алгебры
        """
        pass
    
    @abstractmethod
    def Get_base_ring(self):
        """
        Возвращает базовое кольцо алгебры
        """
        pass
