"""
Класс дифференцирования
"""

from sage.all import Matrix

class Derivation:
    """
    Класс дифференцирования алгебры
    """

    def __init__(self, algebra, gen_mapping=None):
        """
        Инициализация дифференцирования
        
        Args:
            algebra (Algebra): Алгебра на которой определено дифференцирование
            gen_mapping (dict, optional): Словарь образов генераторов {g: D(g)}.
                                          Если None, производная считается нулевой.
        """
        self.algebra = algebra
        self.gen_mapping = gen_mapping if gen_mapping is not None else {}


    def __call__(self, element):
        """
        Применение дифференцирования к элементу алгебры.
        Использует линейность и правило Лейбница: D(f) = sum (df/dx_i) * D(x_i).
        """
        # Handle scalars (int, float, sage types) which have derivative 0
        if hasattr(element, 'is_constant') and element.is_constant():
             # It's likely a polynomial constant
             pass
        elif isinstance(element, (int, float)):
             return 0

        gens = self.algebra.Get_gens()
        res = 0
        
        for g in gens:
            if g in self.gen_mapping:
                # D(x_i)
                val_on_gen = self.gen_mapping[g]
                # df/dx_i
                try:
                    partial = element.derivative(g)
                except AttributeError:
                    # If element is constant/scalar but missed by checks above, derivative is 0
                    # We should ensure it's not a garbage type (like str)
                    if isinstance(element, (str, list, dict, set)):
                         raise TypeError(f"Элемент {element} не поддерживается.")
                    partial = 0
                
                res += partial * val_on_gen
                
        return res


    def bracket(self, other):
        """
        Вычисление скобки Ли [self, other].
        Результат - новое дифференцирование D = D1*D2 - D2*D1.
        """
        if self.algebra != other.algebra:
            
            pass

        gens = self.algebra.Get_gens()
        new_mapping = {}
        
        for g in gens:
            # [D1, D2](g) = D1(D2(g)) - D2(D1(g))
            val = self(other(g)) - other(self(g))
            new_mapping[g] = val
            
        return Derivation(self.algebra, new_mapping)



    @property
    def is_triangular(self):
        """
        Проверяет, является ли дифференцирование треугольным.
        D треугольное, если D(x_i) зависит только от x_{i+1}, ..., x_n.
        Для x_n производная должна быть константой.
        """
        gens = self.algebra.Get_gens()
        n = len(gens)
        
        for i, g in enumerate(gens):
            # i is index of x_i (0..n-1). 
            # allowed vars are gens[i+1:]
            allowed = set(gens[i+1:])
            
            # Check D(g) variables
            val = self(g)
            
            # Constants are always fine (depend on empty set)
            if hasattr(val, 'is_constant') and val.is_constant():
                continue
            if isinstance(val, (int, float)):
                continue
                
            # If polynomial, check variables
            try:
                # Sage polynomials have .variables() method
                actual_vars = set(val.variables())
                if not actual_vars.issubset(allowed):
                    return False
            except AttributeError:
                # If we assume it works on polynomials, this might be scalar or something else
                # If simple scalar, loop continues. If unexpected object, fail safe?
                pass
                
        return True

    @staticmethod
    def from_jacobian(algebra, polynomials):
        """
        Создает Якобиево дифференцирование на основе набора из n-1 полиномов.
        D(x_j) = (-1)^(n+j+1) * det(J_j), где J_j - матрица Якоби без j-го столбца.
        
        Args:
            algebra: Алгебра
            polynomials: список из n-1 полиномов
        """
        gens = algebra.Get_gens()
        n = len(gens)
        if len(polynomials) != n - 1:
            raise ValueError(f"Требуется n-1 полиномов ({n-1}), получено {len(polynomials)}.")
            
        # 1. Compute full Jacobian matrix (n-1) x n
        # rows: polys, cols: gens. J[i][k] = d(poly_i)/d(gen_k)
        full_jac = []
        for p in polynomials:
            row = [p.derivative(g) for g in gens]
            full_jac.append(row)
            
        # 2. Compute mapping for new derivation
        gen_mapping = {}
        
        for j, g in enumerate(gens):
            # Form submatrix by removing column j
            # Matrix constructor: Matrix(ring, nrows, ncols, entries) or Matrix(list_of_lists)
            # Need to specify ring? generic Matrix can infer.
            
            sub_matrix_data = []
            for row in full_jac:
                new_row = row[:j] + row[j+1:]
                sub_matrix_data.append(new_row)
                
            M_j = Matrix(sub_matrix_data)
            det = M_j.determinant()
            
            # Sign: (-1)^(n + (j+1)) using 1-based indexing for formula D(x_j)
            # D = det( ... matrix ... row of partials )
            # position of partial_xj is (n, j+1). Cofactor is (-1)^(n + j+1) * det(M_j)
            
            sign = (-1)**(n + j + 1)
            gen_mapping[g] = sign * det
            
        return Derivation(algebra, gen_mapping)


class LinearDerivation(Derivation):
    """
    Класс линейного дифференцирования, заданного матрицей A.
    D(x_i) = sum(A[i,j] * x_j)
    """

    def __init__(self, algebra, matrix):
        """
        Инициализация линейного дифференцирования
        
        Args:
            algebra (Algebra): Алгебра
            matrix (Matrix): Квадратная матрица n x n
        """
        self.matrix = matrix
        gens = algebra.Get_gens()
        n = len(gens)
        
        # Check matrix dimensions
        if matrix.nrows() != n or matrix.ncols() != n:
             raise ValueError(f"Матрица должна быть {n}x{n}, получено {matrix.nrows()}x{matrix.ncols()}")
             
        gen_mapping = {}
        # D(x_i) = sum_j A_ij * x_j
        for i, xi in enumerate(gens):
            val = 0
            for j, xj in enumerate(gens):
                val += matrix[i, j] * xj
            gen_mapping[xi] = val
        
        super().__init__(algebra, gen_mapping)


class WeitzenbockDerivation(LinearDerivation):
    """
    Класс дифференцирования Вайтценбёка (линейное с нильпотентной матрицей).
    """
    
    def __init__(self, algebra, matrix):
        """
        Инициализация дифференцирования Вайтценбёка.
        
        Args:
            algebra (Algebra): Алгебра
            matrix (Matrix): Нильпотентная матрица
        """
        if not matrix.is_nilpotent():
            raise ValueError("Матрица должна быть нильпотентной для дифференцирования Вайтценбёка")
            
        super().__init__(algebra, matrix)
