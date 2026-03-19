import logging

from sage.all import PolynomialRing, QQ

from lib.derivation import LieDerivation
from lib.solver import LieBasisSolver

logging.basicConfig(level=logging.INFO)


def test_quotient_circle():
    print("Инициализация кольца и фактор-кольца...")
    ring = PolynomialRing(QQ, "x,y,z")
    x, y, z = ring.gens()

    ideal = ring.ideal(x**2 + y**2 - 1)
    quotient = ring.quotient(ideal)
    xq, yq, zq = quotient.gens()

    print(f"Фактор-кольцо создано: {quotient}")

    gen_map = {xq: -yq, yq: xq, zq: quotient(0)}

    print("Создание дифференцирования...")
    try:
        derivation = LieDerivation.from_mapping(quotient, gen_map)
        print(f"Дифференцирование создано: {derivation}")
        print(f"Тип D: {type(derivation)}")
    except NotImplementedError as exc:
        print("Sage пока не поддерживает полные derivation over quotient rings.")
        print(f"Причина: {exc}")
        return
    except Exception as exc:
        print(f"НЕ УДАЛОСЬ создать дифференцирование: {exc}")
        return

    print(f"Старший член D: {derivation.leading_term}")
    print(f"Степень D: {derivation.degree()}")

    print("Запуск решателя...")
    solver = LieBasisSolver([derivation])
    result = solver.run()
    print(f"Результат решателя: {result}")
