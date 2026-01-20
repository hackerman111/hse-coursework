from .derivation import Derivation
from functools import reduce
import operator
from sage.all import QQ

def get_Beldiev(algebra):
    """
    Возвращает генераторы (U, V) для алгебры Ли векторных полей на K^n согласно статье Бельдиева.
    """
    gens = algebra.gens()
    n = len(gens)
    if n < 2:
        raise ValueError("Теорема Бельдиева требует n >= 2")
    
    zn = gens[n-1]
    
    # U = d/dzn
    u_map = {g: 0 for g in gens}
    u_map[zn] = 1
    U = Derivation(algebra, u_map)
    
    # V
    v_map = {g: 0 for g in gens}
    
    # Слагаемые для d/dzk (k < n-1)
    # k - индекс переменной, по которой дифференцируем (от 0 до n-2)
    # Соответствует z_1...z_{n-1}
    for k in range(n-1):
        tail_vars = gens[k+1:] 
        power = 4 * (k + 2)
        
        if tail_vars:
            prod = reduce(operator.mul, tail_vars, 1)
        else:
            prod = 1
            
        coeff = prod**power
        v_map[gens[k]] = coeff

    # Слагаемое для d/dzn
    full_prod = reduce(operator.mul, gens, 1)
    v_map[gens[n-1]] = full_prod**4
    
    V = Derivation(algebra, v_map)
    
    return [U, V]

def get_Andristy(algebra):
    """
    Возвращает генераторы (U, V, W) согласно статье Андриста.
    """
    gens = algebra.gens()
    n = len(gens)
    if n < 2:
        raise ValueError("Теорема Андриста требует n >= 2")
        
    zn = gens[n-1]
    
    # U = d/dzn
    u_map = {g: 0 for g in gens}
    u_map[zn] = 1
    U = Derivation(algebra, u_map)
    
    # V
    # Формула: V = d/dzn + z_n^3 * d/dz_{n-1} + z_n * z_{n-1}^3 * d/dz_{n-2} + ...
    # Общий член для d/dz_k (где z_k это gens[k-1]):
    # Coeff = (z_n * ... * z_{k+2}) * (z_{k+1})^3
    v_map = {g: 0 for g in gens}
    v_map[zn] = 1
    
    # Цикл по целевым переменным от z_{n-1} до z_1 (индексы n-2 ... 0)
    for i in range(n - 2, -1, -1):
        # i - индекс переменной, по которой идет дифференцирование (target)
        # В формуле участвуют переменные "выше" i.
        # Ближайшая переменная выше (gens[i+1]) входит в кубе.
        next_var = gens[i+1]
        
        # Остальные переменные выше (от i+2 до конца) входят в 1 степени
        tail_vars = gens[i+2:]
        if tail_vars:
            tail_prod = reduce(operator.mul, tail_vars, 1)
        else:
            tail_prod = 1
            
        coeff = tail_prod * (next_var**3)
        v_map[gens[i]] = coeff
        
    V = Derivation(algebra, v_map)
    
    # W
    # Формула: W = (z_1^2 ... z_{n-1}^2) * z_n * d/dzn
    w_map = {g: 0 for g in gens}
    
    head_vars = gens[:n-1]
    if head_vars:
        prod_head = reduce(operator.mul, head_vars, 1)
    else:
        prod_head = 1
        
    w_coeff = (prod_head**2) * zn
    # Направлен вдоль z_n (gens[n-1])
    w_map[gens[n-1]] = w_coeff
    
    W = Derivation(algebra, w_map)
    
    return [U, V, W]


def _deriv_add(D1, D2):
    """Вспомогательная функция сложения дифференцирований"""
    alg = D1.algebra
    gens = alg.gens()
    new_map = {}
    for g in gens:
        new_map[g] = D1(g) + D2(g)
    return Derivation(alg, new_map)

def _deriv_scale(D, scalar):
    """Вспомогательная функция умножения дифференцирования на скаляр"""
    alg = D.algebra
    gens = alg.gens()
    new_map = {}
    for g in gens:
        new_map[g] = scalar * D(g)
    return Derivation(alg, new_map)

def _deriv_sub_scaled(D1, D2, scalar):
    """Вспомогательная функция D1 - scalar * D2"""
    alg = D1.algebra
    gens = alg.gens()
    new_map = {}
    for g in gens:
        new_map[g] = D1(g) - scalar * D2(g)
    return Derivation(alg, new_map)

def _get_leading_term(D):
    """
    Возвращает старший член дифференцирования для редукции.
    Формат: (индекс_компоненты, мономиальный_ключ, коэффициент)
    Приоритет: сначала выбираем компоненту с "наибольшим" мономом.
    Порядок переменных: стандартный Sage.
    """
    gens = D.algebra.gens()
    best_term = None

    for i, g in enumerate(gens):
        poly = D(g)
        if poly == 0:
            continue
            
        lm = poly.lm()
        lc = poly.lc()
        
        if best_term is None:
            best_term = (i, lm, lc)
        else:
            curr_lm = best_term[1]
            if lm > curr_lm:
                best_term = (i, lm, lc)
            elif lm == curr_lm:
                if i > best_term[0]:
                    best_term = (i, lm, lc)
                    
    return best_term

def check(generators, max_iter=100):
    """
    Проверяет, порождают ли переданные генераторы базис частных производных (d/dx_1, ..., d/dx_n).
    Использует алгоритм, аналогичный вычислению базиса Грёбнера (замыкание относительно скобки Ли).
    
    Аргументы:
        generators: список объектов Derivation
        max_iter: ограничение на количество итераций (защита от зацикливания в бесконечной алгебре)
        
    Возвращает:
        True, если найдены все частные производные (константы).
        False, если превышен лимит итераций.
    """
    if not generators:
        return False
        
    algebra = generators[0].algebra
    gens = algebra.gens()
    n = len(gens)
    
    # Базис найденных элементов в "ступенчатом" виде.
    # Ключ: (component_index, monomial) -> Derivation
    # Храним элементы, нормированные на старший коэффициент.
    basis = {} 
    
    queue = list(generators)
    
    # Список уже обработанных элементов для формирования пар
    processed = []

    # Целевые элементы: d/dz_i (индекс i, моном 1)
    targets_found = {i: False for i in range(n)}
    
    iter_count = 0
    
    while queue and iter_count < max_iter:
        iter_count += 1
        current_D = queue.pop(0)
        
        # 1. Редукция
        while True:
            lt = _get_leading_term(current_D)
            if lt is None:
                break
                
            comp_idx, lm, lc = lt
            key = (comp_idx, lm)
            
            if key in basis:
                basis_D = basis[key]
                # D = D - (lc / 1) * basis_D (т.к. базисные нормированы)
                current_D = _deriv_sub_scaled(current_D, basis_D, lc)
            else:
                # Не можем редуцировать этот старший член
                break
        
        lt = _get_leading_term(current_D)
        if lt is None:
            continue
            
        # 2. Нормализация и добавление в базис
        comp_idx, lm, lc = lt
        norm_D = _deriv_scale(current_D, 1/lc)
        basis[(comp_idx, lm)] = norm_D
        
        # Проверка на целевые элементы (чистые частные производные)
        if lm == 1:
            targets_found[comp_idx] = True
            if all(targets_found.values()):
                print(f"Успех! Найдены все частные производные за {iter_count} шагов.")
                return True

        # 3. Генерация новых пар (коммутаторов)
        for old_D in processed:
            new_bracket = norm_D.bracket(old_D)
            queue.append(new_bracket)
            
        processed.append(norm_D)
        
    print(f"Остановка проверки: найдено {sum(targets_found.values())}/{n} производных. Лимит итераций.")
    return False