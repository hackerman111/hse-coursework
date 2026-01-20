"""
Класс дифференцирования
"""

from sage.all import Matrix
from functools import reduce

class _Derivation:
    """
    Класс дифференцирования алгебры (внутренняя реализация)

    Примеры:
        >>> from sage.all import QQ, PolynomialRing
        >>> R = PolynomialRing(QQ, 'x,y')
        >>> x, y = R.gens()
        >>> # D(x) = y, D(y) = x
        >>> D = _Derivation(R, {x: y, y: x})
        >>> D(x)
        y
        >>> D(x**2 + x*y)
        x^2 + 2*x*y + y^2
    """

    def __init__(self, algebra, gen_mapping=None):
        """
        Инициализация дифференцирования

        Аргументы:
            algebra: Алгебра (SageMath PolynomialRing), на которой определено дифференцирование
            gen_mapping (dict, optional): Словарь образов генераторов {g: D(g)}.
                                          Если None, производная считается нулевой.

        Примеры:
            >>> from sage.all import QQ, PolynomialRing
            >>> R = PolynomialRing(QQ, 'x,y')
            >>> D = _Derivation(R)
            >>> D(R.gen(0))
            0
        """
        self.algebra = algebra
        self.gen_mapping = gen_mapping if gen_mapping is not None else {}


    def __call__(self, element):
        """
        Применение дифференцирования к элементу алгебры.
        Использует линейность и правило Лейбница: D(f) = sum (df/dx_i) * D(x_i).

        Примеры:
            >>> from sage.all import QQ, PolynomialRing
            >>> x, y = PolynomialRing(QQ, 'x,y').gens()
            >>> R = x.parent()
            >>> D = _Derivation(R, {x: y, y: 2*x})
            >>> # Линейность: D(f + g) == D(f) + D(g)
            >>> D(x**2 + y) == D(x**2) + D(y)
            True
            >>> # Правило Лейбница: D(f*g) == f*D(g) + g*D(f)
            >>> D(x*y) == x*D(y) + y*D(x)
            True
            >>> # Константа: D(c) == 0
            >>> D(R(5))
            0
            >>> D("not a polynomial")
            Traceback (most recent call last):
                ...
            TypeError: Элемент not a polynomial не поддерживается.
        """
        if hasattr(element, 'is_constant') and element.is_constant():
             pass
        elif isinstance(element, (int, float)):
             return 0

        gens = self.algebra.gens()
        res = 0
        
        for g in gens:
            if g in self.gen_mapping:
                # D(x_i)
                val_on_gen = self.gen_mapping[g]
                # df/dx_i
                try:
                    partial = element.derivative(g)
                except AttributeError:
                    # Запасной вариант для типов, которые могут не иметь метода derivative, но являются константами
                    if isinstance(element, (str, list, dict, set)):
                         raise TypeError(f"Элемент {element} не поддерживается.")
                    partial = 0
                
                res += partial * val_on_gen
                
        return res


    def bracket(self, other):
        """
        Вычисление скобки Ли [self, other].
        Результат - новое дифференцирование D = D1*D2 - D2*D1.

        Примеры:
            >>> from sage.all import QQ, PolynomialRing
            >>> R = PolynomialRing(QQ, 'x,y')
            >>> x, y = R.gens()
            >>> D1 = _Derivation(R, {x: y, y: 0})
            >>> D2 = _Derivation(R, {x: 0, y: x})
            >>> bracket = D1.bracket(D2)
            >>> bracket(x)
            -x
            >>> bracket(y)
            y
            >>> # Антикоммутативность: [D1, D2] = -[D2, D1]
            >>> bracket_inv = D2.bracket(D1)
            >>> bracket(x) == -bracket_inv(x)
            True
            >>> # Тождество Якоби: [A, [B, C]] + [B, [C, A]] + [C, [A, B]] = 0
            >>> A = _Derivation(R, {x: 1, y: 0})
            >>> B = _Derivation(R, {x: y, y: 0})
            >>> C = _Derivation(R, {x: 0, y: x})
            >>> res = A.bracket(B.bracket(C))(x) + B.bracket(C.bracket(A))(x) + C.bracket(A.bracket(B))(x)
            >>> res
            0
        """
        # Сравнение алгебр может потребовать уточнения в зависимости от того, как сравниваются объекты Sage
        if self.algebra != other.algebra:
            pass # Или вызвать ошибку

        gens = self.algebra.gens()
        new_mapping = {}
        
        for g in gens:
            # [D1, D2](g) = D1(D2(g)) - D2(D1(g))
            val = self(other(g)) - other(self(g))
            new_mapping[g] = val
            
        return _Derivation(self.algebra, new_mapping)



    @property
    def is_triangular(self):
        """
        Проверяет, является ли дифференцирование треугольным.
        D треугольное, если D(x_i) зависит только от x_{i+1}, ..., x_n.
        Для x_n производная должна быть константой.

        Примеры:
            >>> from sage.all import QQ, PolynomialRing
            >>> R = PolynomialRing(QQ, 'x,y,z')
            >>> x, y, z = R.gens()
            >>> # D(x) = y + z**2, D(y) = 3*z + 1, D(z) = 5
            >>> D_tri = _Derivation(R, {x: y + z**2, y: 3*z + 1, z: 5})
            >>> D_tri.is_triangular
            True
            >>> # Не треугольное: D(y) зависит от x
            >>> D_bad = _Derivation(R, {x: z, y: x, z: 0})
            >>> D_bad.is_triangular
            False
            >>> # Не треугольное: D(z) зависит от z
            >>> D_bad2 = _Derivation(R, {x: y, y: z, z: z})
            >>> D_bad2.is_triangular
            False
        """
        gens = self.algebra.gens()
        n = len(gens)
        
        for i, g in enumerate(gens):
            allowed = set(gens[i+1:])
            
            val = self(g)
            
            if hasattr(val, 'is_constant') and val.is_constant():
                continue
            if isinstance(val, (int, float)):
                continue
                
            try:
                actual_vars = set(val.variables())
                if not actual_vars.issubset(allowed):
                    return False
            except AttributeError:
                pass
                
        return True


def Derivation(algebra, gen_mapping=None):
    """
    Фабричная функция для создания общего дифференцирования.
    
    Аргументы:
        algebra: SageMath PolynomialRing
        gen_mapping (dict, optional): отображение {генератор: образ}

    Примеры:
        >>> from sage.all import QQ, PolynomialRing
        >>> R = PolynomialRing(QQ, 'x,y')
        >>> x, y = R.gens()
        >>> D = Derivation(R, {x: y, y: x})
        >>> D(x)
        y
    """
    return _Derivation(algebra, gen_mapping)

def LinearDerivation(algebra, matrix):
    """
    Фабричная функция для создания линейного дифференцирования из матрицы.
    D(x_i) = sum(A[i,j] * x_j)
    
    Аргументы:
        algebra: Алгебра SageMath
        matrix: Квадратная матрица n x n

    Примеры:
        >>> from sage.all import QQ, PolynomialRing, Matrix
        >>> R = PolynomialRing(QQ, 'x,y,z')
        >>> x, y, z = R.gens()
        >>> mat = Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
        >>> D = LinearDerivation(R, mat)
        >>> D(x)
        y
        >>> D(y)
        x
        >>> D(z)
        0
    """
    gens = algebra.gens()
    n = len(gens)
    
    if matrix.nrows() != n or matrix.ncols() != n:
            raise ValueError(f"Матрица должна быть {n}x{n}, получено {matrix.nrows()}x{matrix.ncols()}")
            
    gen_mapping = {}
    # D(x_i) = sum_j A_ij * x_j
    for i, xi in enumerate(gens):
        val = 0
        for j, xj in enumerate(gens):
            val += matrix[i, j] * xj
        gen_mapping[xi] = val
    
    return _Derivation(algebra, gen_mapping)


def WeitzenbockDerivation(algebra, matrix):
    """
    Фабричная функция для создания дифференцирования Вайтценбёка (линейное с нильпотентной матрицей).

    Примеры:
        >>> from sage.all import QQ, PolynomialRing, Matrix
        >>> R = PolynomialRing(QQ, 'x,y,z')
        >>> x, y, z = R.gens()
        >>> mat = Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
        >>> D = WeitzenbockDerivation(R, mat)
        >>> D(x)
        y
        >>> D.is_triangular
        True
        >>> non_nil = Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]])
        >>> WeitzenbockDerivation(R, non_nil)
        Traceback (most recent call last):
            ...
        ValueError: Матрица должна быть нильпотентной для дифференцирования Вайтценбёка
    """
    if not matrix.is_nilpotent():
        raise ValueError("Матрица должна быть нильпотентной для дифференцирования Вайтценбёка")
        
    return LinearDerivation(algebra, matrix) # LinearDerivation возвращает _Derivation


def JacobianDerivation(algebra, polynomials):
    """
    Фабричная функция для создания Якобиева дифференцирования.
    D(x_j) = (-1)^(n+j+1) * det(J_j)
    
    Аргументы:
        algebra: Алгебра SageMath
        polynomials: список из n-1 полиномов

    Примеры:
        >>> from sage.all import QQ, PolynomialRing
        >>> R = PolynomialRing(QQ, 'x,y,z')
        >>> x, y, z = R.gens()
        >>> f1 = x**2 + y**2
        >>> f2 = z
        >>> D = JacobianDerivation(R, [f1, f2])
        >>> D(x)
        2*y
        >>> D(y)
        -2*x
        >>> D(z)
        0
    """
    gens = algebra.gens()
    n = len(gens)
    if len(polynomials) != n - 1:
        raise ValueError(f"Требуется n-1 полиномов ({n-1}), получено {len(polynomials)}.")
        
    #J[i][k] = d(poly_i)/d(gen_k)
    full_jac = []
    for p in polynomials:
        row = [p.derivative(g) for g in gens]
        full_jac.append(row)
        
    gen_mapping = {}
    
    for j, g in enumerate(gens):
        
        sub_matrix_data = []
        for row in full_jac:
            new_row = row[:j] + row[j+1:]
            sub_matrix_data.append(new_row)
            
        M_j = Matrix(sub_matrix_data)
        det = M_j.determinant()
        
        
        sign = (-1)**(n + j + 1)
        gen_mapping[g] = sign * det
        
    return _Derivation(algebra, gen_mapping)
