import unittest
from sage.all import PolynomialRing, QQ
from lib.generator import get_Beldiev, get_Andristy, check

class TestLieGenerators(unittest.TestCase):
    """
    Тестирование генераторов алгебры Ли полиномиальных векторных полей.
    Проверка основана на алгоритме замыкания относительно скобки Ли (функция check).
    """

    def setUp(self):
        # Используем поле рациональных чисел для точности и скорости
        self.base_field = QQ

    def test_beldiev_theorem_n2(self):
        """
        Проверка Теоремы 3.1 (Бельдиев) для размерности n=2.
        Согласно статье, алгебра порождается двумя элементами U и V.
        """
        print("\n--- Тест Бельдиева (n=2) ---")
        # Создаем кольцо полиномов K[z1, z2]
        R = PolynomialRing(self.base_field, 'z1, z2')
        
        # Получаем генераторы U, V
        generators = get_Beldiev(R)
        
        # Проверяем количество генераторов (должно быть 2)
        self.assertEqual(len(generators), 2, "Теорема Бельдиева утверждает наличие ровно 2 генераторов")
        
        # Проверяем, что они порождают базис частных производных
        # Это гарантирует, что Lie(U, V) = W_2
        is_generated = check(generators, max_iter=200)
        self.assertTrue(is_generated, "Генераторы Бельдиева не породили базис производных для n=2")

    def test_beldiev_theorem_n3(self):
        """
        Проверка Теоремы 3.1 (Бельдиев) для размерности n=3.
        Проверяем масштабируемость для n > 2.
        """
        print("\n--- Тест Бельдиева (n=3) ---")
        R = PolynomialRing(self.base_field, 'z1, z2, z3')
        generators = get_Beldiev(R)
        
        self.assertEqual(len(generators), 2)
        
        # Увеличиваем лимит итераций, так как для n=3 нужно больше коммутаторов
        is_generated = check(generators, max_iter=500)
        self.assertTrue(is_generated, "Генераторы Бельдиева не породили базис производных для n=3")

    def test_andrist_theorem_n2(self):
        """
        Проверка Теоремы 7 (Андрист) для размерности n=2.
        Согласно статье, алгебра порождается тремя полными векторными полями U, V, W.
        """
        print("\n--- Тест Андриста (n=2) ---")
        R = PolynomialRing(self.base_field, 'z1, z2')
        
        # Получаем генераторы U, V, W
        generators = get_Andristy(R)
        
        # Проверяем количество (должно быть 3)
        self.assertEqual(len(generators), 3, "Теорема Андриста утверждает наличие 3 генераторов")
        
        is_generated = check(generators, max_iter=200)
        self.assertTrue(is_generated, "Генераторы Андриста не породили базис производных для n=2")

    def test_andrist_theorem_n3(self):
        """
        Проверка Теоремы 7 (Андрист) для размерности n=3.
        """
        print("\n--- Тест Андриста (n=3) ---")
        R = PolynomialRing(self.base_field, 'z1, z2, z3')
        generators = get_Andristy(R)
        
        self.assertEqual(len(generators), 3)
        
        is_generated = check(generators, max_iter=500)
        self.assertTrue(is_generated, "Генераторы Андриста не породили базис производных для n=3")

    def test_argument_validation(self):
        """
        Проверка обработки некорректных входных данных.
        Статьи указывают ограничения на n.
        Андрист: n >= 2.
        """
        print("\n--- Тест валидации аргументов ---")
        R_small = PolynomialRing(self.base_field, 'z1') # n=1
        
        with self.assertRaises(ValueError):
            get_Beldiev(R_small)
            
        with self.assertRaises(ValueError):
            get_Andristy(R_small)

if __name__ == '__main__':
    unittest.main()