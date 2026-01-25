from dataclasses import dataclass, field
import heapq
import logging
from typing import List, Tuple, Optional, Any, Dict
from .derivation import LieDerivation

logger = logging.getLogger(__name__)

@dataclass(order=True)
class QueueItem:
    degree: int
    unique_counter: int
    derivation: Any = field(compare=False)

class LieBasisSolver:
    """
    Класс, реализующий алгоритм проверки полноты набора дифференцирований
    (аналог поиска базиса Грёбнера для алгебры Ли векторных полей).
    """
    def __init__(self, generators: List[LieDerivation], max_iter: int = 1000):
        if not generators:
            raise ValueError("Список генераторов не может быть пустым")
        
        self.generators = generators
        self.max_iter = max_iter
        self.algebra = generators[0].codomain()
        self.gens = self.algebra.gens()
        self.n = len(self.gens)
        
        # Состояние
        self.basis: Dict[Tuple[int, Any], LieDerivation] = {} # {(comp_idx, lm): Derivation}
        self.queue: List[QueueItem] = []
        self.processed: List[LieDerivation] = []
        self.targets_found: Dict[int, bool] = {i: False for i in range(self.n)}
        
        self.unique_counter = 0
        self.iter_count = 0
        
        # Инициализация очереди
        for g in generators:
            self._add_to_queue(g)

    def _add_to_queue(self, D: LieDerivation):
        deg = D.degree()
        # min-heap по степени, затем по порядку добавления
        item = QueueItem(deg, self.unique_counter, D)
        heapq.heappush(self.queue, item)
        self.unique_counter += 1

    def _reduce(self, D: LieDerivation) -> LieDerivation:
        """Редуцирует элемент относительно текущего базиса."""
        while True:
            lt = D.leading_term
            if lt is None:
                break
                
            comp_idx, lm, lc = lt
            key = (comp_idx, lm)
            
            if key in self.basis:
                basis_D = self.basis[key]
                # D = D - (lc/basis_lc) * basis_D. 
                # Базисные элементы нормированы (lc=1), поэтому D - lc * basis_D
                D = D - (basis_D * lc)
            else:
                break
        return D

    def _generate_commutators(self, new_D: LieDerivation):
        """Добавляет коммутаторы [new_D, old_D] для всех old_D в processed."""
        for old_D in self.processed:
            new_bracket = new_D.bracket(old_D)
            self._add_to_queue(new_bracket)

    def step(self) -> bool:
        """
        Выполняет один шаг алгоритма.
        Возвращает True, если нужно продолжать, False если закончили.
        """
        if not self.queue or self.iter_count >= self.max_iter:
            return False
            
        self.iter_count += 1
        item = heapq.heappop(self.queue)
        current_D = item.derivation
        
        # 1. Редукция
        current_D = self._reduce(current_D)
        
        lt = current_D.leading_term
        if lt is None:
            return True # Элемент редуцировался в 0, продолжаем
            
        # 2. Нормализация и добавление в базис
        comp_idx, lm, lc = lt
        # Нормируем, чтобы старший коэффициент был 1
        norm_D = current_D * (1/lc)
        self.basis[(comp_idx, lm)] = norm_D
        
        logger.debug(f"Итерация {self.iter_count}: добавлен элемент с LT = d/d{self.gens[comp_idx]} * {lm}")

        # Проверка целевых производных (d/dx_i)
        # Они соответствуют моному 1
        if lm == 1:
            if not self.targets_found[comp_idx]:
                self.targets_found[comp_idx] = True
                logger.info(f"Найдена производная d/d{self.gens[comp_idx]}")
                
            if all(self.targets_found.values()):
                logger.info(f"Успех! Найдены все {self.n} частных производных за {self.iter_count} шагов.")
                return False # Успех, останавливаемся

        # 3. Генерация пар
        self._generate_commutators(norm_D)
        self.processed.append(norm_D)
        
        return True

    def run(self) -> bool:
        """Запускает цикл вычислений."""
        logger.info(f"Начало проверки базиса. Генераторов: {len(self.generators)}, переменных: {self.n}")
        
        running = True
        while running:
            running = self.step()
            
        # По завершении проверяем результат
        found_count = sum(self.targets_found.values())
        success = all(self.targets_found.values())
        
        if not success:
            logger.warning(f"Остановка проверки: найдено {found_count}/{self.n}. Очередь пуста или лимит {self.max_iter}.")
            
        return success
