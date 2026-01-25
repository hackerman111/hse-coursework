from .derivation import LieDerivation
from .solver import LieBasisSolver
from functools import reduce
import operator
import logging
from typing import List, Any

# Настройка логгера
logger = logging.getLogger(__name__)

def get_Beldiev(algebra: Any) -> List[LieDerivation]:
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
    U = LieDerivation.from_mapping(algebra, u_map)
    
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
    
    V = LieDerivation.from_mapping(algebra, v_map)
    
    return [U, V]

def get_Andristy(algebra: Any) -> List[LieDerivation]:
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
    U = LieDerivation.from_mapping(algebra, u_map)
    
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
        
    V = LieDerivation.from_mapping(algebra, v_map)
    
    # W
    w_map = {g: 0 for g in gens}
    head_vars = gens[:n-1]
    if head_vars:
        prod_head = reduce(operator.mul, head_vars, 1)
    else:
        prod_head = 1
        
    w_coeff = (prod_head**2) * zn
    w_map[gens[n-1]] = w_coeff
    
    W = LieDerivation.from_mapping(algebra, w_map)
    
    return [U, V, W]



def check(generators: List[LieDerivation], max_iter: int = 1000) -> bool:
    """
    Оптимизированная проверка (аналог базиса Грёбнера).
    Использует LieBasisSolver.
    """
    solver = LieBasisSolver(generators, max_iter)
    return solver.run()