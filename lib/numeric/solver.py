"""
Основной решатель: проверка Lie(G_1, …, G_m) ⊇ W_n^{(≤d)}.

Алгоритм Buchberger-style: скобки вычисляются между полными
(неоднородными) элементами, а не между степенными пулами.
Оптимизация new-vs-existing: каждый новый элемент коммутируется
только с уже обработанными.  SVD-базис для устойчивого ранга.
"""

from __future__ import annotations

import time
from typing import Dict, List, Optional, Tuple, TYPE_CHECKING

import numpy as np

from .components import decompose_generator
from .indexing import dim_Wn_k
from .results import DegreeStatus, NumericResult
from .structure import StructureConstantCache
from .subspace import OrthonormalBasis
from .types import Generator
from ..utils import LibraryLogger


# ---------------------------------------------------------------------------
# Основной класс
# ---------------------------------------------------------------------------

class NumericLieGenerationChecker:
    """
    Проверяет: Lie(G_1, …, G_m) ⊇ W_n^{(≤d)}?

    Параметры
    ---------
    n : число переменных
    d : проверяем порождённость до этой степени включительно
    generators : список генераторов
        Каждый генератор — dict {i: {alpha: coeff}},
        где i (0-based) — номер переменной ∂_{z_i},
        alpha — кортеж-мультииндекс.
    D_max : максимальная степень для усечения (если None, берётся max(d, max deg генераторов))
    verbose : выводить промежуточные результаты
    log_every_s : раз в сколько секунд печатать полную сводку; None/0 отключает
    logger : переиспользуемый логгер библиотеки; если не задан, создаётся автоматически
    """

    def __init__(
        self,
        n: int,
        d: int,
        generators: List[Generator],
        D_max: Optional[int] = None,
        verbose: bool = True,
        log_every_s: Optional[float] = None,
        logger: Optional[LibraryLogger] = None,
        structure_constants: Optional[StructureConstantCache] = None,
    ) -> None:
        self.n = n
        self.d = d
        self.generators = generators
        self.logger = logger if logger is not None else LibraryLogger(enabled=verbose)
        self.verbose = verbose
        self.log_every_s = log_every_s if log_every_s and log_every_s > 0 else None

        # Определяем D_max
        max_gen_deg = -1
        for gen in generators:
            for i, mono_dict in gen.items():
                for alpha in mono_dict:
                    k = sum(alpha) - 1
                    if k > max_gen_deg:
                        max_gen_deg = k
        self.D_max = D_max if D_max is not None else max(d, max_gen_deg)

        # Инициализация подпространств B[k] для k = -1 … D_max
        self.bases: Dict[int, OrthonormalBasis] = {}
        for k in range(-1, self.D_max + 1):
            N_k = dim_Wn_k(n, k)
            self.bases[k] = OrthonormalBasis(N_k)

        # Разложить генераторы в полные элементы и засеять подпространства
        self._elements: List[Dict[int, np.ndarray]] = []
        for gen in generators:
            elem = decompose_generator(gen, n, self.D_max)
            self._elements.append(elem)
            for k, vec in elem.items():
                if k in self.bases:
                    self.bases[k].try_add(vec)

        self.structure_constants = (
            structure_constants
            if structure_constants is not None
            else StructureConstantCache(n=n, D_max=self.D_max)
        )
        if self.structure_constants.n != n:
            raise ValueError(
                "structure_constants cache was built for a different number of variables"
            )
        if self.structure_constants.D_max < self.D_max:
            raise ValueError(
                "structure_constants cache does not cover the requested truncation degree"
            )

    def run(self) -> NumericResult:
        """Запустить итеративное замыкание и вернуть результат."""
        from .bracket import bracket_full_elements

        t0 = time.perf_counter()
        iteration = 0
        next_detail_log_at = (
            t0 + self.log_every_s
            if self.verbose and self.logger.enabled and self.log_every_s is not None
            else None
        )

        if self.verbose and self.logger.enabled:
            self._print_status(iteration)

        # Алгоритм «новые-против-существующих» (Buchberger-style):
        # processed — элементы, уже прошедшие через все скобки;
        # pending   — элементы, которые ещё нужно прокоммутировать.
        processed: List[Dict[int, np.ndarray]] = []
        pending: List[Dict[int, np.ndarray]] = list(self._elements)

        while pending:
            iteration += 1
            next_pending: List[Dict[int, np.ndarray]] = []

            # Скобки: каждый pending × (все processed + остальные pending)
            for i, elem_new in enumerate(pending):
                # С каждым processed
                for elem_old in processed:
                    next_detail_log_at = self._maybe_print_detailed_status(
                        iteration, t0, next_detail_log_at,
                    )
                    bracket = bracket_full_elements(
                        elem_new, elem_old, self.n, self.structure_constants,
                    )
                    self._try_add_bracket(bracket, next_pending)

                # С другими pending (только i < j, чтобы не дублировать)
                for j in range(i + 1, len(pending)):
                    next_detail_log_at = self._maybe_print_detailed_status(
                        iteration, t0, next_detail_log_at,
                    )
                    bracket = bracket_full_elements(
                        elem_new, pending[j], self.n, self.structure_constants,
                    )
                    self._try_add_bracket(bracket, next_pending)

            processed.extend(pending)
            pending = next_pending

            if self.verbose and self.logger.enabled:
                self._print_status(iteration)
            next_detail_log_at = self._maybe_print_detailed_status(
                iteration, t0, next_detail_log_at,
            )

            # Ранняя остановка: все целевые степени полны
            if all(
                self.bases[k].rank >= dim_Wn_k(self.n, k)
                for k in range(-1, self.d + 1)
            ):
                break

        elapsed = time.perf_counter() - t0

        # Собираем результат
        degrees: List[DegreeStatus] = []
        success = True
        for k in range(-1, self.d + 1):
            N_k = dim_Wn_k(self.n, k)
            rank = self.bases[k].rank
            ds = DegreeStatus(degree=k, rank=rank, target_dim=N_k)
            degrees.append(ds)
            if not ds.is_full:
                success = False

        return NumericResult(
            n=self.n,
            d=self.d,
            D_max=self.D_max,
            success=success,
            degrees=degrees,
            iterations=iteration,
            elapsed_s=elapsed,
        )

    def _try_add_bracket(
        self,
        bracket: Dict[int, np.ndarray],
        next_pending: List[Dict[int, np.ndarray]],
    ) -> None:
        """Попробовать добавить компоненты скобки в базисы.

        Если хотя бы одна компонента увеличила ранг соответствующего
        подпространства, весь элемент добавляется в next_pending
        для дальнейших итераций.
        """
        if not bracket:
            return
        added_any = False
        for k, vec in bracket.items():
            if k not in self.bases:
                continue
            norm = np.linalg.norm(vec)
            if norm < 1e-15:
                continue
            normalized = vec / norm
            if self.bases[k].try_add(normalized):
                added_any = True
        if added_any:
            next_pending.append(bracket)

    def _print_status(self, iteration: int) -> None:
        """Вывести краткий статус."""
        parts = []
        for k in range(-1, min(self.d, self.D_max) + 1):
            N_k = dim_Wn_k(self.n, k)
            r = self.bases[k].rank
            if r >= N_k:
                parts.append(f"[{k}:✓]")
            else:
                parts.append(f"[{k}:{r}/{N_k}]")
        status = " ".join(parts)
        self.logger.info(f"iter {iteration:3d}: {status}")

    def _maybe_print_detailed_status(
        self,
        iteration: int,
        started_at: float,
        next_log_at: Optional[float],
        active_degree: Optional[int] = None,
        active_pair: Optional[Tuple[int, int]] = None,
    ) -> Optional[float]:
        """Периодически печатать полную сводку по всем степеням."""
        return self.logger.maybe_emit_periodic(
            next_log_at=next_log_at,
            interval_s=self.log_every_s,
            message_factory=lambda: self._build_detailed_status(
                iteration=iteration,
                elapsed_s=time.perf_counter() - started_at,
                active_degree=active_degree,
                active_pair=active_pair,
            ),
        )

    def _build_detailed_status(  # noqa: E301
        self,
        iteration: int,
        elapsed_s: float,
        active_degree: Optional[int] = None,
        active_pair: Optional[Tuple[int, int]] = None,
    ) -> str:
        """Построить подробную сводку по состоянию замыкания."""
        target_full = 0
        closure_full = 0
        lines = [
            "=== Progress Summary ===",
            f"elapsed = {elapsed_s:.3f}s, iteration = {iteration}",
            f"target range = -1..{self.d}, closure range = -1..{self.D_max}",
        ]

        if active_degree is not None:
            lines.append(f"active degree = {active_degree}")
        if active_pair is not None:
            lines.append(f"active pair = ({active_pair[0]}, {active_pair[1]})")

        for k in range(-1, self.D_max + 1):
            N_k = dim_Wn_k(self.n, k)
            rank = self.bases[k].rank
            is_full = rank >= N_k
            if k <= self.d and is_full:
                target_full += 1
            if is_full:
                closure_full += 1
            tag = "FULL" if is_full else f"{rank}/{N_k}"
            scope = "target" if k <= self.d else "closure"
            lines.append(f"Degree {k:3d} [{scope:7s}]: {tag}")

        lines.insert(
            3,
            f"full layers: target {target_full}/{self.d + 2}, closure {closure_full}/{self.D_max + 2}",
        )
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Публичные примитивы алгоритма
# ---------------------------------------------------------------------------


def try_add_to_basis(basis: OrthonormalBasis, vectors: np.ndarray) -> int:
    """
    Нормализовать строки и добавить новые в ортонормальный базис.

    Это извлечённый строительный блок из NumericLieGenerationChecker.run():
    нормализация + фильтрация нулевых строк + добавление через SVD.

    Args:
        basis:   ортонормальный базис подпространства
        vectors: матрица строк для добавления

    Returns:
        Количество добавленных новых базисных векторов.
    """
    if vectors.shape[0] == 0:
        return 0
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    mask = norms.ravel() > 1e-15
    if not mask.any():
        return 0
    normalized = vectors[mask] / norms[mask]
    return basis.try_add_many(normalized)


def check_degree_fullness(
    bases: "Dict[int, OrthonormalBasis]",
    n: int,
    d: int,
) -> "List[DegreeStatus]":
    """
    Сформировать список статусов для степеней -1..d.

    Извлечённый строительный блок из конца NumericLieGenerationChecker.run().

    Args:
        bases: словарь {degree: OrthonormalBasis}
        n:     размерность пространства
        d:     максимальная проверяемая степень

    Returns:
        Список DegreeStatus по каждой степени.
    """
    result: List[DegreeStatus] = []
    for k in range(-1, d + 1):
        N_k = dim_Wn_k(n, k)
        rank = bases[k].rank if k in bases else 0
        result.append(DegreeStatus(degree=k, rank=rank, target_dim=N_k))
    return result
