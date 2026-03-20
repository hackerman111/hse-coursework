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
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .bracket import batch_bracket, precompute_structure_constants
from .indexing import dim_Wn_k, dim_Wn_leq_d, multiindices, multiindex_lookup
from .subspace import OrthonormalBasis
from ..utils import LibraryLogger


# ---------------------------------------------------------------------------
# Результат проверки
# ---------------------------------------------------------------------------

@dataclass
class DegreeStatus:
    """Статус для одной степени k."""
    degree: int
    rank: int
    target_dim: int

    @property
    def is_full(self) -> bool:
        return self.rank >= self.target_dim

    def __repr__(self) -> str:
        tag = "FULL" if self.is_full else f"{self.rank}/{self.target_dim}"
        return f"Degree {self.degree:3d}: {tag}"


@dataclass
class NumericResult:
    """Результат работы NumericLieGenerationChecker."""
    n: int
    d: int
    D_max: int
    success: bool
    degrees: List[DegreeStatus] = field(default_factory=list)
    iterations: int = 0
    elapsed_s: float = 0.0

    def summary(self) -> str:
        lines = [
            f"=== Numeric Lie Generation Check ===",
            f"n = {self.n}, d = {self.d}, D_max = {self.D_max}",
            f"Итераций: {self.iterations}, Время: {self.elapsed_s:.3f}s",
            f"Результат: {'PASS — все степени ≤ d полны' if self.success else 'FAIL — есть неполные степени'}",
            "",
        ]
        for ds in self.degrees:
            lines.append(str(ds))
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Представление генератора
# ---------------------------------------------------------------------------

def _decompose_generator(
    gen_dict: Dict[int, Dict[Tuple[int, ...], float]],
    n: int,
    D_max: int,
) -> Dict[int, np.ndarray]:
    """
    Разложить генератор на однородные компоненты.

    gen_dict : {i: {alpha: coeff}} — отображение ∂_{z_i} → сумма coeff · z^α
               i — 0-based номер переменной
               alpha — кортеж-мультииндекс, |α| = (степень мономиального поля) + 1

    Возвращает {k: vector_in_Wn_k} для каждой ненулевой компоненты степени k.
    """
    components: Dict[int, np.ndarray] = {}

    for i, mono_dict in gen_dict.items():
        for alpha, coeff in mono_dict.items():
            s = sum(alpha)
            k = s - 1  # degree
            if k < -1 or k > D_max:
                continue
            N_k = dim_Wn_k(n, k)
            if k not in components:
                components[k] = np.zeros(N_k, dtype=np.float64)
            lookup = multiindex_lookup(n, s)
            m = len(lookup)
            idx = i * m + lookup[alpha]
            components[k][idx] += coeff

    return components


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
        generators: List[Dict[int, Dict[Tuple[int, ...], float]]],
        D_max: Optional[int] = None,
        verbose: bool = True,
        log_every_s: Optional[float] = None,
        logger: Optional[LibraryLogger] = None,
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
            components = _decompose_generator(gen, n, self.D_max)
            for k, vec in components.items():
                if k in self.bases:
                    self.bases[k].try_add(vec)

        # Предвычислить структурные константы для всех нужных пар
        self._sc_cache: Dict[Tuple[int, int], Any] = {}
        for p in range(-1, self.D_max + 1):
            for q in range(p, self.D_max + 1):
                r = p + q
                if r < -1 or r > self.D_max:
                    continue
                self._sc_cache[(p, q)] = precompute_structure_constants(n, p, q)

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

                    sc = self._sc_cache.get((p, q))
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


# ---------------------------------------------------------------------------
# Вспомогательные конструкторы генераторов
# ---------------------------------------------------------------------------

def make_beldiev_generators(n: int) -> List[Dict[int, Dict[Tuple[int, ...], float]]]:
    """
    Генераторы Бельдиева для W_n:
      U = ∂_{z_n}
      V = z_n^{4n} ∂_{z_{n-1}} + (z_n z_{n-1})^{4(n-1)} ∂_{z_{n-2}} + …
          + (z_n … z_2)^8 ∂_{z_1} + (z_n … z_1)^4 ∂_{z_n}
    """
    # U = ∂_{z_{n-1}}  (0-based: переменная n-1)
    U: Dict[int, Dict[Tuple[int, ...], float]] = {
        (n - 1): {tuple(0 if j != n - 1 else 0 for j in range(n)): 1.0}
    }
    # α = (0, …, 0) с |α| = 0, то есть z^0 ∂_{z_{n-1}} = ∂_{z_{n-1}}
    # Но |α| должен быть 0, и 0-based, так что alpha = (0,)*n
    U = {(n - 1): {(0,) * n: 1.0}}

    V: Dict[int, Dict[Tuple[int, ...], float]] = {}

    # Компоненты V:
    # Для index = 0, …, n-2 (0-based):
    #   V отображает ∂_{z_index} с коэффициентом (z_{index+1} · … · z_{n-1})^{4(n - index - 1) / ...}
    # По формуле Бельдиева (0-based):
    #   V = z_{n-1}^{4n} ∂_{z_{n-2}}
    #     + (z_{n-1} z_{n-2})^{4(n-1)} ∂_{z_{n-3}}
    #     + …
    #     + (z_{n-1} z_{n-2} … z_1)^8 ∂_{z_0}
    #     + (z_{n-1} z_{n-2} … z_0)^4 ∂_{z_{n-1}}

    # Общий паттерн (1-based в статье, переводим в 0-based):
    # Слагаемое для ∂_{z_{n-2-j}} (j = 0, …, n-2):
    #   коэффициент = (z_{n-1} · z_{n-2} · … · z_{n-1-j})^{4(n-j)}
    # Слагаемое для ∂_{z_{n-1}}:
    #   коэффициент = (z_{n-1} · z_{n-2} · … · z_0)^4

    for j in range(n - 1):
        # ∂_{z_{n-2-j}}: коэфф = (z_{n-1} · … · z_{n-1-j})^{power}
        # tail_vars: индексы n-1, n-2, …, n-1-j  (j+1 переменных)
        target_var = n - 2 - j  # куда идёт производная
        power = 4 * (n - j)
        # Мультииндекс α: α_k = power если k ∈ {n-1, n-2, …, n-1-j}, иначе 0
        alpha = [0] * n
        for k in range(n - 1, n - 2 - j, -1):
            alpha[k] = power
        alpha_tuple = tuple(alpha)

        if target_var not in V:
            V[target_var] = {}
        V[target_var][alpha_tuple] = 1.0

    # Последнее слагаемое: (z_{n-1} · … · z_0)^4 ∂_{z_{n-1}}
    alpha = tuple([4] * n)
    target_var = n - 1
    if target_var not in V:
        V[target_var] = {}
    V[target_var][alpha] = 1.0

    return [U, V]


def make_monomial_generator(
    n: int,
    target_var: int,
    alpha: Tuple[int, ...],
    coeff: float = 1.0,
) -> Dict[int, Dict[Tuple[int, ...], float]]:
    """
    Создать генератор вида coeff · z^α ∂_{z_{target_var}}.
    """
    return {target_var: {alpha: coeff}}


def make_random_generator(
    n: int,
    max_degree: int,
    rng: np.random.Generator | None = None,
    sparsity: float = 0.5,
) -> Dict[int, Dict[Tuple[int, ...], float]]:
    """
    Случайный генератор в W_n с мономами до степени max_degree.

    sparsity : доля ненулевых коэффициентов (0 < sparsity ≤ 1).
    """
    if rng is None:
        rng = np.random.default_rng()

    gen: Dict[int, Dict[Tuple[int, ...], float]] = {}
    for i in range(n):
        mono_dict: Dict[Tuple[int, ...], float] = {}
        for s in range(0, max_degree + 2):  # |α| = s → deg = s-1
            for alpha in multiindices(n, s):
                if rng.random() < sparsity:
                    mono_dict[alpha] = float(rng.standard_normal())
        if mono_dict:
            gen[i] = mono_dict

    return gen
