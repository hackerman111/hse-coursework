from .derivation import Derivation
from functools import reduce
import operator
import heapq
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
    v_map = {g: 0 for g in gens}
    v_map[zn] = 1
    
    # Цикл по целевым переменным
    for i in range(n - 2, -1, -1):
        next_var = gens[i+1]
        tail_vars = gens[i+2:]
        if tail_vars:
            tail_prod = reduce(operator.mul, tail_vars, 1)
        else:
            tail_prod = 1
            
        coeff = tail_prod * (next_var**3)
        v_map[gens[i]] = coeff
        
    V = Derivation(algebra, v_map)
    
    # W
    w_map = {g: 0 for g in gens}
    head_vars = gens[:n-1]
    if head_vars:
        prod_head = reduce(operator.mul, head_vars, 1)
    else:
        prod_head = 1
        
    w_coeff = (prod_head**2) * zn
    w_map[gens[n-1]] = w_coeff
    
    W = Derivation(algebra, w_map)
    
    return [U, V, W]


def _deriv_add(D1, D2):
    alg = D1.algebra
    gens = alg.gens()
    new_map = {}
    for g in gens:
        new_map[g] = D1(g) + D2(g)
    return Derivation(alg, new_map)

def _deriv_scale(D, scalar):
    alg = D.algebra
    gens = alg.gens()
    new_map = {}
    for g in gens:
        new_map[g] = scalar * D(g)
    return Derivation(alg, new_map)

def _deriv_sub_scaled(D1, D2, scalar):
    alg = D1.algebra
    gens = alg.gens()
    new_map = {}
    for g in gens:
        new_map[g] = D1(g) - scalar * D2(g)
    return Derivation(alg, new_map)

def _get_leading_term(D):
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

def check(generators, max_iter=1000):
    """
    Оптимизированная проверка (аналог базиса Грёбнера).
    Использует очередь с приоритетом (heapq) для обработки элементов
    с меньшей степенью в первую очередь.
    """
    if not generators:
        return False
        
    algebra = generators[0].algebra
    gens = algebra.gens()
    n = len(gens)
    
    # Базис найденных элементов
    basis = {} 
    
    # Очередь с приоритетом: (degree, unique_id, derivation)
    # unique_id нужен для стабильности сортировки при равных степенях
    queue = []
    unique_counter = 0

    def get_degree(D):
        """Вычисляет максимальную степень среди компонент векторного поля"""
        max_deg = -1
        for g in gens:
            poly = D(g)
            if poly != 0:
                try:
                    d = poly.degree()
                    if d > max_deg:
                        max_deg = d
                except:
                    pass
        return max_deg

    def add_to_queue(D):
        nonlocal unique_counter
        deg = get_degree(D)
        # Чем меньше степень, тем выше приоритет (min-heap)
        heapq.heappush(queue, (deg, unique_counter, D))
        unique_counter += 1

    for g in generators:
        add_to_queue(g)
    
    processed = []
    targets_found = {i: False for i in range(n)}
    iter_count = 0
    
    while queue and iter_count < max_iter:
        iter_count += 1
        # Извлекаем элемент с наименьшей степенью
        _, _, current_D = heapq.heappop(queue)
        
        # 1. Редукция
        while True:
            lt = _get_leading_term(current_D)
            if lt is None:
                break
                
            comp_idx, lm, lc = lt
            key = (comp_idx, lm)
            
            if key in basis:
                basis_D = basis[key]
                current_D = _deriv_sub_scaled(current_D, basis_D, lc)
            else:
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
                # print(f"Успех! Найдены все частные производные за {iter_count} шагов.")
                return True

        # 3. Генерация новых пар (коммутаторов)
        # Коммутируем только с уже добавленными в базис
        for old_D in processed:
            new_bracket = norm_D.bracket(old_D)
            add_to_queue(new_bracket)
            
        processed.append(norm_D)
        
    # print(f"Остановка проверки: найдено {sum(targets_found.values())}/{n}. Лимит итераций.")
    return False