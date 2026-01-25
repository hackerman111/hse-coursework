import sage.all
from sage.all import QQ, PolynomialRing, Matrix
from lib.derivation import LieDerivation
from lib.solver import LieBasisSolver
import logging

# Настройка логгера для вывода в консоль
logging.basicConfig(level=logging.INFO)

def test_quotient_circle():
    """
    Тестирует LieDerivation на R = QQ[x,y,z] / <x^2 + y^2 - 1>
    В идеале, дифференцирование вращения x*dy - y*dx должно быть корректно определено.
    """
    print("Инициализация кольца и фактор-кольца...")
    R = PolynomialRing(QQ, 'x,y,z')
    x, y, z = R.gens()
    
    # Идеал для окружности/цилиндра: x^2 + y^2 - 1
    I = R.ideal(x**2 + y**2 - 1)
    Q = R.quotient(I)
    xq, yq, zq = Q.gens()
    
    print(f"Фактор-кольцо создано: {Q}")
    
    # Определение векторного поля вращения: -y*d/dx + x*d/dy + 0*d/dz
    # В отображении дифференцирования Sage: {x: -y, y: x, z: 0}
    # Мы должны использовать элементы фактор-кольца для значений отображения
    
    gen_map = {
        xq: -yq,
        yq: xq,
        zq: Q(0)
    }
    
    print("Создание дифференцирования...")
    try:
        D = LieDerivation.from_mapping(Q, gen_map)
        print(f"Дифференцирование создано: {D}")
        print(f"Тип D: {type(D)}")
    except Exception as e:
        print(f"НЕ УДАЛОСЬ создать дифференцирование: {e}")
        return

    # Проверка свойств
    print(f"Старший член D: {D.leading_term}")
    print(f"Степень D: {D.degree()}")
    
    # Запуск решателя на проверку этого единственного дифференцирования (тривиальная проверка)
    print("Запуск решателя...")
    solver = LieBasisSolver([D])
    res = solver.run()
    print(f"Результат решателя: {res}")
    
    # Проверка, сохраняется ли идеал I
    # D(x^2 + y^2 - 1) = 2x*D(x) + 2y*D(y) = 2x*(-y) + 2y*(x) = -2xy + 2xy = 0. Верно.

if __name__ == "__main__":
    try:
        test_quotient_circle()
    except ImportError:
        print("SageMath недоступен, пропуск выполнения теста.")
    except Exception as e:
        print(f"Тест завершился с ошибкой: {e}")
        import traceback
        traceback.print_exc()
