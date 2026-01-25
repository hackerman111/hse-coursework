from sage.all import QQ, PolynomialRing

R = PolynomialRing(QQ, 'x,y')
x, y = R.gens()

# Создание дифференцирования
# D(x) = y, D(y) = 2*x
try:
    D = R.derivation({x: y, y: 2*x})
    print(f"Type of D: {type(D)}")
    print(f"D(x) = {D(x)}")
    print(f"D(x^2 + y) = {D(x**2 + y)}")
    
    # Арифметика
    D2 = R.derivation({x: 1, y: 0})
    D_sum = D + D2
    print(f"(D+D2)(x) = {D_sum(x)}")
    
    # Скобка Ли
    try:
        # Дифференцирования Sage могут иметь метод bracket или lie_bracket
        if hasattr(D, 'bracket'):
            print("Has bracket method")
        if hasattr(D, 'commutator'):
            print("Has commutator method")
            
        # Попробуем вычислить скобку вручную, если необходимо, или проверим, работает ли операция
        # Обычно производная Ли или подобное.
        # Проверим доступные методы
        print(f"Dir of D: {dir(D)}")
        
    except Exception as e:
        print(f"Bracket check failed: {e}")

    # Старший член?
    try:
        if hasattr(D, 'leading_term'):
             print(f"Leading term: {D.leading_term()}")
        else:
             print("No leading_term method")
    except Exception as e:
        print(f"Leading term check failed: {e}")
        
except Exception as e:
    print(f"Creation failed: {e}")
