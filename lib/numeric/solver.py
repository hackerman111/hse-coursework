"""
Основной решатель: проверка Lie(G_1, …, G_m) ⊇ W_n^{(≤d)}.

Реализует алгоритм из Numeric.md §3.4 с оптимизациями:
  - нисходящий обход степеней для каскадов вида ad_U^k(V);
  - нормализация строк со скобками перед пополнением базиса;
  - SVD-базис вместо RREF для устойчивого численного ранга;
  - послойное раннее прекращение (§5.4);
  - кеширование структурных констант.
"""

from __future__ import annotations

import time
from typing import Dict, List, Optional, Tuple

import numpy as np

from .bracket import batch_bracket
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

        # Разложить генераторы и засеять подпространства
        for gen in generators:
            components = decompose_generator(gen, n, self.D_max)
            for k, vec in components.items():
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
        t0 = time.perf_counter()
        iteration = 0
        next_detail_log_at = (
            t0 + self.log_every_s
            if self.verbose and self.logger.enabled and self.log_every_s is not None
            else None
        )

        if self.verbose and self.logger.enabled:
            self._print_status(iteration)

        while True:
            iteration += 1
            changed = False

            # Обрабатываем степени в НИСХОДЯЩЕМ порядке: каскад ad_U
            # идёт от высоких степеней к низким, и при нисходящем порядке
            # новые элементы доступны для вычисления скобок на ТЕКУЩЕЙ итерации.
            for r in range(self.D_max, -2, -1):
                N_r = dim_Wn_k(self.n, r)
                if self.bases[r].is_full(N_r):
                    continue  # Степень r уже полна — скобки не добавят нового

                new_rows_list: List[np.ndarray] = []

                for p in range(max(-1, r - self.D_max), min(self.D_max, r + 1) + 1):
                    q = r - p
                    next_detail_log_at = self._maybe_print_detailed_status(
                        iteration,
                        t0,
                        next_detail_log_at,
                        active_degree=r,
                        active_pair=(p, q),
                    )
                    if q < -1 or q > self.D_max or q < p:
                        continue

                    Bp = self.bases.get(p)
                    Bq = self.bases.get(q)
                    if Bp is None or Bq is None or Bp.rank == 0 or Bq.rank == 0:
                        continue

                    sc = self.structure_constants.get(p, q)
                    if sc is None:
                        continue

                    # Вычисляем все скобки [B_p, B_q]
                    all_p = Bp.all_rows
                    all_q = Bq.all_rows
                    if all_p.shape[0] > 0 and all_q.shape[0] > 0:
                        brackets = batch_bracket(all_p, p, all_q, q, self.n, sc)
                        if brackets.shape[0] > 0:
                            new_rows_list.append(brackets)

                if new_rows_list:
                    all_new = np.vstack(new_rows_list)
                    # Нормализуем строки для числовой устойчивости SVD
                    norms = np.linalg.norm(all_new, axis=1, keepdims=True)
                    mask = norms.ravel() > 1e-15
                    if mask.any():
                        all_new = all_new[mask] / norms[mask]
                        added = self.bases[r].try_add_many(all_new)
                        if added > 0:
                            changed = True

            if self.verbose and self.logger.enabled:
                self._print_status(iteration)
            next_detail_log_at = self._maybe_print_detailed_status(
                iteration,
                t0,
                next_detail_log_at,
            )

            if not changed:
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

    def _build_detailed_status(
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
