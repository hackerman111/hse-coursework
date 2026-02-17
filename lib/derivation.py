"""
Класс дифференцирования
"""

from sage.all import Matrix
from functools import reduce
from typing import Optional, Tuple, Any, Dict, List, Union

class LieDerivation:
    """
    Обертка над дифференцированием SageMath для удобства работы.
    """
    def __init__(self, sage_derivation: Any) -> None:
        self._d = sage_derivation
        self._algebra = sage_derivation.codomain()

    def __call__(self, arg: Any) -> Any:
        return self._d(arg)

    def __add__(self, other: Union['LieDerivation', Any]) -> 'LieDerivation':
        if isinstance(other, LieDerivation):
            return LieDerivationFactory.create(self._d + other._d)
        return LieDerivationFactory.create(self._d + other)

    def __sub__(self, other: Union['LieDerivation', Any]) -> 'LieDerivation':
        if isinstance(other, LieDerivation):
            return LieDerivationFactory.create(self._d - other._d)
        return LieDerivationFactory.create(self._d - other)
        
    def __mul__(self, other: Union['LieDerivation', Any]) -> 'LieDerivation':
        # Умножение на скаляр или полином (справа)
        # Дифференцирование Sage * скаляр работает.
        if isinstance(other, LieDerivation):
             raise TypeError("Нельзя умножать два дифференцирования напрямую. Используйте скобку Ли (bracket).")
        return LieDerivationFactory.create(self._d * other)

    def __rmul__(self, other: Any) -> 'LieDerivation':
        # Скаляр * D
        return LieDerivationFactory.create(other * self._d)

    def codomain(self) -> Any:
        return self._algebra
        
    def bracket(self, other: 'LieDerivation') -> 'LieDerivation':
        """
        [D1, D2] = D1(D2) - D2(D1)
        """
        if not isinstance(other, LieDerivation):
             raise ValueError("Можно вычислять скобку только с другим LieDerivation")
             
        # Оптимизация: использовать явную формулу, если Sage это плохо поддерживает,
        # но переопределение логики из derivation_bracket безопаснее для сохранения того же поведения.
        # Исходная простая логика: 
        gens = self._algebra.gens()
        images = []
        for g in gens:
            val = self(other(g)) - other(self(g))
            images.append(val)
        return LieDerivationFactory.create(self._algebra.derivation(images))

    @property
    def leading_term(self) -> Optional[Tuple[int, Any, Any]]:
        """
        Возвращает старший член дифференцирования (term).
        """
        gens = self._algebra.gens()
        best_term = None

        for i, g in enumerate(gens):
            poly = self(g)
            if poly == 0:
                continue
                
            try:
                lm = poly.lm()
                lc = poly.lc()
            except AttributeError:
                # Константа
                lm = self._algebra(1)
                lc = poly
                if poly == 0: continue

            if best_term is None:
                best_term = (i, lm, lc)
            else:
                curr_lm = best_term[1]
                if lm > curr_lm:
                    best_term = (i, lm, lc)
                elif lm == curr_lm:
                     # При равенстве мономов выбираем компоненту с большим индексом
                    if i > best_term[0]:
                        best_term = (i, lm, lc)
        return best_term

    def is_triangular(self) -> bool:
        gens = self._algebra.gens()
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

    def degree(self) -> int:
        """
        Возвращает максимальную степень представителей коэффициентов.
        """
        gens = self._algebra.gens()
        max_deg = -1
        for g in gens:
            poly = self(g)
            if poly != 0:
                try:
                    d = poly.degree()
                    if d > max_deg:
                        max_deg = d
                except AttributeError:
                    # Для констант (чисел) считаем степень 0
                    if max_deg < 0: max_deg = 0
        return max_deg
        
    @staticmethod
    def from_mapping(algebra, gen_mapping=None):
        """
        Создает дифференцирование по отображению.
        """
        if gen_mapping is None:
            gen_mapping = {}
        gens = algebra.gens()
        images = [gen_mapping.get(g, 0) for g in gens]
        return LieDerivationFactory.create(algebra.derivation(images))

    @staticmethod
    def from_linear(algebra, matrix):
        """
        Создает линейное дифференцирование D(x_i) = sum(A[i,j] * x_j).
        """
        gens = algebra.gens()
        n = len(gens)
        
        if matrix.nrows() != n or matrix.ncols() != n:
                raise ValueError(f"Матрица должна быть {n}x{n}, получено {matrix.nrows()}x{matrix.ncols()}")
                
        images = []
        for i, xi in enumerate(gens):
            val = 0
            for j, xj in enumerate(gens):
                val += matrix[i, j] * xj
            images.append(val)
        
        return LieDerivationFactory.create(algebra.derivation(images))

    @staticmethod
    def from_weitzenbock(algebra, matrix):
        """
        Создает дифференцирование Вайтценбёка (линейное с нильпотентной матрицей).
        """
        if not matrix.is_nilpotent():
            raise ValueError("Матрица должна быть нильпотентной для дифференцирования Вайтценбёка")
            
        return LieDerivation.from_linear(algebra, matrix)

    @staticmethod
    def from_jacobian(algebra, polynomials):
        """
        Создает Якобиево дифференцирование.
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
            
        images = []
        
        for j, g in enumerate(gens):
            sub_matrix_data = []
            for row in full_jac:
                new_row = row[:j] + row[j+1:]
                sub_matrix_data.append(new_row)
            
            # В случае необходимости можно использовать алгебру для создания матрицы, или generic Matrix. Аргумент 'matrix' подразумевает использование sage.all.Matrix.
            M_j = Matrix(sub_matrix_data)
            det = M_j.determinant()
            
            sign = (-1)**(n + j + 1)
            images.append(sign * det)
            
        return LieDerivationFactory.create(algebra.derivation(images))

    def __repr__(self):
        return f"LieD({self._d})"


class QuotientLieDerivation(LieDerivation):
    """
    Специализированная версия дифференцирования для работы в фактор-кольцах (R/I).
    Автоматически поднимает (lift) элементы до исходного кольца для вычисления
    старших членов и степеней.
    """
    
    def __init__(self, sage_derivation: Any) -> None:
        super().__init__(sage_derivation)
        # Проверка, действительно ли это фактор-кольцо
        if not hasattr(self._algebra, 'cover'):
             raise ValueError("QuotientLieDerivation требует фактор-алгебру (с методом cover())")
        self._base_ring = self._algebra.cover()

    def _lift_poly(self, element: Any) -> Any:
        """Поднимает элемент из фактор-кольца в исходное кольцо многочленов."""
        if hasattr(element, 'lift'):
            return element.lift()
        return element

    @property
    def leading_term(self) -> Optional[Tuple[int, Any, Any]]:
        """
        Возвращает старший член, основываясь на представителе в исходном кольце.
        Возвращает: (index, lm_representative, lc_in_quotient)
        """
        gens = self._algebra.gens()
        best_term = None

        for i, g in enumerate(gens):
            val = self(g)
            if val == 0:
                continue
            
            # --- ЛОГИКА ФАКТОР-КОЛЬЦА ---
            # Поднимаем элемент, чтобы получить корректный мономиальный порядок
            lifted_val = self._lift_poly(val)
            
            try:
                lm = lifted_val.lm()
                # Коэффициент берем из lift, но возвращаем его проекцию в фактор-кольцо,
                # так как арифметика должна оставаться в R/I
                lc_lifted = lifted_val.lc()
                lc = self._algebra(lc_lifted) 
            except AttributeError:
                # Константа
                lm = self._base_ring(1) # Единица базового кольца для сравнения
                lc = val

            if best_term is None:
                best_term = (i, lm, lc)
            else:
                curr_lm = best_term[1]
                # Сравниваем мономы в базовом кольце (они уникальны)
                if lm > curr_lm:
                    best_term = (i, lm, lc)
                elif lm == curr_lm:
                    if i > best_term[0]:
                        best_term = (i, lm, lc)
                        
        return best_term

    def degree(self) -> int:
        """
        Возвращает максимальную степень представителей коэффициентов.
        Необходима для корректной сортировки в очереди Solver.
        """
        gens = self._algebra.gens()
        max_deg = -1
        for g in gens:
            val = self(g)
            if val != 0:
                lifted = self._lift_poly(val)
                try:
                    d = lifted.degree()
                    if d > max_deg:
                        max_deg = d
                except AttributeError:
                    if max_deg < 0: max_deg = 0
        return max_deg
        
    def __repr__(self):
        return f"QuotientLieD({self._d})"


class LieDerivationFactory:
    @staticmethod
    def create(sage_derivation: Any) -> LieDerivation:
        algebra = sage_derivation.codomain()
        # Проверяем наличие метода cover (признак фактор-кольца в Sage)
        if hasattr(algebra, 'cover'):
            return QuotientLieDerivation(sage_derivation)
        else:
            return LieDerivation(sage_derivation)





