"""
Решатель для проверки полноты набора дифференцирований.
"""

from __future__ import annotations

import heapq
import logging
from typing import Any, Callable, Dict, List, Optional, Tuple

from ..derivation import LieDerivation
from .operations import generate_commutators, reduce_derivation
from .queue_item import QueueItem

logger = logging.getLogger(__name__)


class LieBasisSolver:
    """
    Аналог поиска базиса Грёбнера для алгебры Ли векторных полей.
    """

    def __init__(
        self,
        generators: List[LieDerivation],
        max_iter: int = 1000,
        stop_condition: Optional[Callable[["LieBasisSolver"], bool]] = None,
    ):
        if not generators:
            raise ValueError("Список генераторов не может быть пустым")

        self.generators = generators
        self.max_iter = max_iter
        self.stop_condition = stop_condition
        self.algebra = generators[0].codomain()
        self.gens = self.algebra.gens()
        self.n = len(self.gens)

        self.basis: Dict[Tuple[int, Any], LieDerivation] = {}
        self.queue: List[QueueItem] = []
        self.processed: List[LieDerivation] = []
        self.targets_found: Dict[int, bool] = {index: False for index in range(self.n)}

        self.unique_counter = 0
        self.iter_count = 0

        for generator in generators:
            self._add_to_queue(generator)

    def _add_to_queue(self, derivation: LieDerivation) -> None:
        item = QueueItem(derivation.degree(), self.unique_counter, derivation)
        heapq.heappush(self.queue, item)
        self.unique_counter += 1

    def _reduce(self, derivation: LieDerivation) -> LieDerivation:
        return reduce_derivation(derivation, self.basis)

    def _generate_commutators(self, derivation: LieDerivation) -> None:
        for bracket in generate_commutators(derivation, self.processed):
            self._add_to_queue(bracket)

    def step(self) -> bool:
        if not self.queue or self.iter_count >= self.max_iter:
            return False

        self.iter_count += 1
        current = heapq.heappop(self.queue).derivation
        current = self._reduce(current)

        leading_term = current.leading_term
        if leading_term is None:
            return True

        comp_idx, monomial, coeff = leading_term
        normalized = current * (1 / coeff)
        self.basis[(comp_idx, monomial)] = normalized

        logger.debug(
            "Итерация %s: добавлен элемент с LT = d/d%s * %s",
            self.iter_count,
            self.gens[comp_idx],
            monomial,
        )

        if self.stop_condition is not None and self.stop_condition(self):
            logger.info(
                "Остановка по stop_condition на итерации %s.",
                self.iter_count,
            )
            return False

        if monomial == 1:
            if not self.targets_found[comp_idx]:
                self.targets_found[comp_idx] = True
                logger.info("Найдена производная d/d%s", self.gens[comp_idx])

            if all(self.targets_found.values()):
                logger.info(
                    "Успех! Найдены все %s частных производных за %s шагов.",
                    self.n,
                    self.iter_count,
                )
                return False

        self._generate_commutators(normalized)
        self.processed.append(normalized)
        return True

    def run(self) -> bool:
        logger.info(
            "Начало проверки базиса. Генераторов: %s, переменных: %s",
            len(self.generators),
            self.n,
        )

        running = True
        while running:
            running = self.step()

        success = all(self.targets_found.values())
        if not success:
            found_count = sum(self.targets_found.values())
            logger.warning(
                "Остановка проверки: найдено %s/%s. Очередь пуста или лимит %s.",
                found_count,
                self.n,
                self.max_iter,
            )
        return success
