import unittest
from sage.all import PolynomialRing, QQ
from lib.generator import check
from lib.derivation import Derivation

class TestCheckAlgorithm(unittest.TestCase):
    """
    Тестирование самой функции проверки (check).
    Мы должны быть уверены, что она корректно определяет, порождают ли элементы алгебру.
    """
    
    def setUp(self):
        self.R = PolynomialRing(QQ, 'x, y')
        self.x, self.y = self.R.gens()

    def _make_deriv(self, mapping):
        """Вспомогательный метод для создания дифференцирования"""
        # Заполняем нулями отсутствующие ключи для удобства
        full_map = {self.x: 0, self.y: 0}
        full_map.update(mapping)
        return Derivation(self.R, full_map)

    def test_check_with_full_basis(self):
        """Тест 1: Подаем готовый базис {d/dx, d/dy}. Должен вернуть True."""
        dx = self._make_deriv({self.x: 1})
        dy = self._make_deriv({self.y: 1})
        
        result = check([dx, dy])
        self.assertTrue(result, "Функция check должна возвращать True для полного базиса")

    def test_check_with_incomplete_basis(self):
        """Тест 2: Подаем только {d/dx}. Должен вернуть False."""
        dx = self._make_deriv({self.x: 1})
        
        result = check([dx])
        self.assertFalse(result, "Функция check должна возвращать False, если базис не полон")

    def test_check_with_commutator_generation(self):
        """
        Тест 3: Подаем {d/dx, x*d/dy}.
        Коммутатор: [d/dx, x*d/dy] = 1 * d/dy.
        Таким образом, мы получаем d/dy и имеем полный набор. Должен вернуть True.
        """
        dx = self._make_deriv({self.x: 1})
        x_dy = self._make_deriv({self.y: self.x}) # x * d/dy
        
        result = check([dx, x_dy])
        self.assertTrue(result, "Функция check должна находить производные через коммутаторы")

    def test_check_with_subalgebra(self):
        """
        Тест 4: Подаем {x*d/dx, y*d/dy} (Эйлеровы поля).
        [x*dx, y*dy] = 0.
        Чистые производные d/dx и d/dy никогда не будут получены. Должен вернуть False.
        """
        x_dx = self._make_deriv({self.x: self.x})
        y_dy = self._make_deriv({self.y: self.y})
        
        result = check([x_dx, y_dy])
        self.assertFalse(result, "Функция check не должна возвращать True для подалгебры, не содержащей константных производных")