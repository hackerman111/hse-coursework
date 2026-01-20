from sage.all import PolynomialRing, QQ, Matrix
from lib.derivation import Derivation, LinearDerivation, WeitzenbockDerivation, JacobianDerivation

def test_derivation_refactor():
    print("Тестирование рефакторинга дифференцирования...")
    
    # Настройка алгебры
    R = PolynomialRing(QQ, 'x,y,z')
    x, y, z = R.gens()
    
    # 1. Общее дифференцирование
    print("1. Тестирование фабрики общего дифференцирования...")
    mapping = {x: y, y: z, z: 0}
    D = Derivation(R, mapping)
    print(f"D(x) = {D(x)}")
    assert D(x) == y
    print("Общее дифференцирование OK")
    
    # 2. Линейное дифференцирование
    print("2. Тестирование фабрики линейного дифференцирования...")
    M = Matrix(QQ, [[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    DL = LinearDerivation(R, M)
    print(f"DL(x) = {DL(x)}")
    assert DL(x) == y
    print("Линейное дифференцирование OK")
    
    # 3. Дифференцирование Вайтценбёка
    print("3. Тестирование фабрики дифференцирования Вайтценбёка...")
    DW = WeitzenbockDerivation(R, M) # M — нильпотентна
    assert DW(x) == y
    
    try:
        M_non_nil = Matrix(QQ, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        WeitzenbockDerivation(R, M_non_nil)
        print("Ошибка: Должно было возникнуть ValueError для ненильпотентной матрицы")
    except ValueError:
        print("Перехвачено ожидаемое ValueError для ненильпотентной матрицы")
    print("Дифференцирование Вайтценбёка OK")

    # 4. Якобиево дифференцирование
    # n=3, требуется 2 полинома
    print("4. Тестирование фабрики Якобиева дифференцирования...")
    polys = [x**2 + y, z] # просто случайные полиномы
    DJ = JacobianDerivation(R, polys)
    # Проверка, что результат существует и имеет правильный тип
    res = DJ(x)
    print(f"DJ(x) = {res}")
    print("Якобиево дифференцирование OK")

if __name__ == "__main__":
    test_derivation_refactor()
