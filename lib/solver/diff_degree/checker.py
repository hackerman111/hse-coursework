"""
Алгоритм GENCHECK в точной символьной форме.
"""

from __future__ import annotations

import time
from typing import Iterable

from ...derivation import LieDerivation
from ...utils import LibraryLogger
from .dimensions import dim_Wn_coeff_k, dim_Wn_coeff_leq_d
from .result import DiffDegreeLayerStatus, DiffDegreeResult
from .tracker import InsertResult, SymbolicSubspaceTracker, TrackedDerivation


class DiffDegreeLieGenerationChecker:
    """
    Символьная проверка ``Lie(G) ⊇ W_n^{<=d}`` в градуировке по степени
    коэффициентов из ``Diff_deg_D.md``.
    """

    def __init__(
        self,
        generators: Iterable[LieDerivation],
        target_degree: int,
        *,
        work_degree: int | None = None,
        max_iter: int = 5000,
        verbose: bool = True,
        log_every_s: float | None = 5.0,
        logger: LibraryLogger | None = None,
    ) -> None:
        self.generators = list(generators)
        if not self.generators:
            raise ValueError("Список генераторов не может быть пустым")

        self.algebra = self.generators[0].codomain()
        for generator in self.generators[1:]:
            if generator.codomain() != self.algebra:
                raise ValueError("Все генераторы должны лежать в одной алгебре")

        self.variables = self.algebra.gens()
        self.n = len(self.variables)
        self.target_degree = target_degree
        self.max_iter = max_iter
        self.verbose = verbose
        self.log_every_s = None if log_every_s is None or log_every_s <= 0 else log_every_s
        self.logger = logger if logger is not None else LibraryLogger(enabled=verbose)

        max_generator_degree = max(generator.degree() for generator in self.generators)
        base_work_degree = max(target_degree, max_generator_degree)
        self.work_degree = base_work_degree if work_degree is None else max(
            work_degree, base_work_degree
        )

        self.target_dim = dim_Wn_coeff_leq_d(self.n, self.target_degree)
        self.state = SymbolicSubspaceTracker(self.algebra, self.work_degree)
        self.partials = {
            axis: self._make_partial_derivative(axis)
            for axis in range(self.n)
        }

        self.iterations = 0
        self._sufficient_condition_hit = False
        self._next_log_at: float | None = None
        self._start_time = 0.0

    def run(self) -> DiffDegreeResult:
        self._start_time = time.perf_counter()
        if self.log_every_s is not None:
            self._next_log_at = self._start_time + self.log_every_s

        self._log_phase("F0", "initialization")
        for index, generator in enumerate(self.generators):
            self._insert(generator, source=f"G{index + 1}")
            if self._is_success():
                return self._build_result("F0")

        success = self._phase_one()
        if success:
            return self._build_result("F1")

        success = self._phase_two()
        if success:
            return self._build_result("F2")

        success = self._phase_three()
        if success:
            return self._build_result("F3")

        success = self._phase_four()
        if success:
            return self._build_result("F4")

        return self._build_result("F5")

    def _phase_one(self) -> bool:
        self._log_phase("F1", "peeling using discovered partial derivatives")
        _, success = self._peel_fixpoint()
        if success:
            return True

        for left_index, left in enumerate(self.generators):
            for right_index in range(left_index + 1, len(self.generators)):
                if self._limit_reached():
                    return self._is_success()

                right = self.generators[right_index]
                self.iterations += 1
                outcome = self._insert(
                    left.bracket(right),
                    source=f"[G{left_index + 1}, G{right_index + 1}]",
                )
                if outcome.tracked is not None:
                    _, success = self._peel_fixpoint([outcome.tracked])
                    if success:
                        return True
                elif self._is_success():
                    return True

                self._maybe_emit_periodic("F1", active=f"[G{left_index + 1}, G{right_index + 1}]")

        return self._is_success()

    def _phase_two(self) -> bool:
        if not self.state.partial_axes:
            return False

        self._log_phase("F2", "degree-1 extraction and closure")
        for element in tuple(self.state.elements):
            if element.degree < 2:
                continue
            for axis in sorted(self.state.partial_axes):
                if self._limit_reached():
                    return self._is_success()

                self.iterations += 1
                derived = self.partials[axis].bracket(element.derivation)
                if derived.degree() != 1:
                    continue

                outcome = self._insert(derived, source=f"deg1 peel #{element.id} by z{axis}")
                if outcome.partial_added:
                    _, success = self._peel_fixpoint()
                    if success:
                        return True
                if self._is_success():
                    return True

        changed = True
        while changed and not self._limit_reached():
            changed = False
            degree_one = tuple(self.state.elements_of_degree(1))
            for left_index, left in enumerate(degree_one):
                for right in degree_one[left_index + 1 :]:
                    self.iterations += 1
                    bracket = left.derivation.bracket(right.derivation)
                    if bracket.degree() != 1:
                        continue
                    outcome = self._insert(
                        bracket,
                        source=f"deg1 bracket #{left.id} #{right.id}",
                    )
                    if self._mark_change(outcome):
                        changed = True
                    if self._is_success():
                        return True

            for pivot in range(self.n):
                pivot_monomial = self._unit(pivot)
                for src in range(self.n):
                    if not self.state.has_monomial(self._unit(src), pivot):
                        continue
                    for dst in range(self.n):
                        if src == dst:
                            continue
                        if not self.state.has_monomial(pivot_monomial, dst):
                            continue
                        outcome = self._insert_monomial(
                            dst,
                            self._unit(src),
                            source=f"cross term z{src}*d{dst} via z{pivot}",
                        )
                        if self._mark_change(outcome):
                            changed = True
                        if self._is_success():
                            return True

        return self._is_success()

    def _phase_three(self) -> bool:
        if self.state.degree_rank(1) < self.n * self.n:
            return False

        self._log_phase("F3", "symbolic saturation by explicit identities")
        for degree in range(2, self.target_degree + 1):
            changed = True
            while changed and not self._limit_reached():
                changed = False

                if self._mark_change(self._symbolic_power_raise(degree)):
                    changed = True
                if self._is_success():
                    return True

                if self._mark_change(self._symbolic_direction_change(degree)):
                    changed = True
                if self._is_success():
                    return True

                if self._mark_change(self._symbolic_variable_mix(degree)):
                    changed = True
                if self._is_success():
                    return True

                layer_changed, success = self._gl_action(degree)
                changed |= layer_changed
                if success:
                    return True

                if self.state.degree_rank(degree) >= dim_Wn_coeff_k(self.n, degree):
                    break

        return self._is_success()

    def _phase_four(self) -> bool:
        self._log_phase("F4", "fallback bracket BFS")
        queue = sorted(self.state.elements, key=lambda item: (item.degree, item.id))
        visited: set[tuple[int, int]] = set()
        stagnant_rounds = 0
        last_rank = self.state.rank_upto(self.target_degree)

        cursor = 0
        while cursor < len(queue) and not self._limit_reached():
            left = queue[cursor]
            cursor += 1
            growth = False

            for right in tuple(self.state.elements):
                pair = (min(left.id, right.id), max(left.id, right.id))
                if pair in visited:
                    continue
                visited.add(pair)

                if left.degree + right.degree - 1 > self.work_degree:
                    continue

                self.iterations += 1
                outcome = self._insert(
                    left.derivation.bracket(right.derivation),
                    source=f"fallback bracket #{left.id} #{right.id}",
                )
                if outcome.tracked is not None:
                    queue.append(outcome.tracked)
                    queue.sort(key=lambda item: (item.degree, item.id))
                    growth = True
                    peel_changed, success = self._peel_fixpoint([outcome.tracked])
                    growth |= peel_changed
                    if success:
                        return True
                elif outcome.partial_added:
                    growth = True

                if self._is_success():
                    return True
                self._maybe_emit_periodic("F4", active=f"[{left.id}, {right.id}]")

            current_rank = self.state.rank_upto(self.target_degree)
            if not growth and current_rank == last_rank:
                stagnant_rounds += 1
            else:
                stagnant_rounds = 0
                last_rank = current_rank

            if stagnant_rounds >= 2:
                break

        return self._is_success()

    def _peel_fixpoint(
        self,
        seeds: Iterable[TrackedDerivation] | None = None,
    ) -> tuple[bool, bool]:
        queue = list(seeds if seeds is not None else self.state.elements)
        queued_ids = {item.id for item in queue}
        processed: set[tuple[int, int]] = set()
        changed = False
        cursor = 0

        while cursor < len(queue) and not self._limit_reached():
            element = queue[cursor]
            cursor += 1

            for axis in sorted(self.state.partial_axes):
                key = (element.id, axis)
                if key in processed:
                    continue
                processed.add(key)

                current = element.derivation
                while current.degree() > 0 and not self._limit_reached():
                    self.iterations += 1
                    current = self.partials[axis].bracket(current)
                    outcome = self._insert(
                        current,
                        source=f"peel #{element.id} by z{axis}",
                    )
                    if self._mark_change(outcome):
                        changed = True
                    if outcome.tracked is not None and outcome.tracked.id not in queued_ids:
                        queue.append(outcome.tracked)
                        queued_ids.add(outcome.tracked.id)
                    if outcome.partial_added:
                        for existing in self.state.elements:
                            if existing.id not in queued_ids:
                                queue.append(existing)
                                queued_ids.add(existing.id)
                    if self._is_success():
                        return changed, True
                    self._maybe_emit_periodic("F1", active=f"peel #{element.id} by z{axis}")

        return changed, False

    def _symbolic_power_raise(self, limit_degree: int) -> InsertResult:
        result = InsertResult(None, False, False)
        for axis in range(self.n):
            base = [0] * self.n
            base[axis] = 2
            if not self.state.has_monomial(tuple(base), axis):
                continue

            for power in range(3, limit_degree):
                current = [0] * self.n
                current[axis] = power
                if not self.state.has_monomial(tuple(current), axis):
                    continue

                current[axis] += 1
                result = self._merge_outcomes(
                    result,
                    self._insert_monomial(
                        axis,
                        tuple(current),
                        source=f"power raise z{axis}^{power + 1}*d{axis}",
                    ),
                )
        return result

    def _symbolic_direction_change(self, limit_degree: int) -> InsertResult:
        result = InsertResult(None, False, False)
        for axis in range(self.n):
            source_key = self._unit(axis)
            for target in range(self.n):
                if axis == target:
                    continue
                if not self.state.has_monomial(source_key, target):
                    continue

                for power in range(0, limit_degree + 1):
                    exponents = [0] * self.n
                    exponents[axis] = power
                    if not self.state.has_monomial(tuple(exponents), axis):
                        continue
                    result = self._merge_outcomes(
                        result,
                        self._insert_monomial(
                            target,
                            tuple(exponents),
                            source=f"direction change z{axis}^{power}*d{target}",
                        ),
                    )
        return result

    def _symbolic_variable_mix(self, limit_degree: int) -> InsertResult:
        result = InsertResult(None, False, False)
        for exponents, component in tuple(self.state.monomials):
            current_degree = sum(exponents)
            if current_degree >= limit_degree:
                continue

            lifted = list(exponents)
            lifted[component] += 1
            if not self.state.has_monomial(tuple(lifted), component):
                continue

            for axis in range(self.n):
                if axis == component:
                    continue

                for power in range(1, limit_degree - current_degree + 1):
                    pure = [0] * self.n
                    pure[axis] = power
                    if not self.state.has_monomial(tuple(pure), component):
                        continue

                    mixed = list(exponents)
                    mixed[axis] += power
                    result = self._merge_outcomes(
                        result,
                        self._insert_monomial(
                            component,
                            tuple(mixed),
                            source=f"mix z{axis}^{power} into d{component}",
                        ),
                    )
        return result

    def _gl_action(self, degree: int) -> tuple[bool, bool]:
        changed = False
        while not self._limit_reached():
            local_changed = False
            layer = tuple(self.state.elements_of_degree(degree))
            degree_one = tuple(self.state.elements_of_degree(1))
            if not layer or not degree_one:
                return changed, False

            for linear in degree_one:
                for element in layer:
                    self.iterations += 1
                    bracket = linear.derivation.bracket(element.derivation)
                    if bracket.degree() != degree:
                        continue
                    outcome = self._insert(
                        bracket,
                        source=f"gl action #{linear.id} on #{element.id}",
                    )
                    if self._mark_change(outcome):
                        changed = True
                        local_changed = True
                    if self._is_success():
                        return changed, True
            if not local_changed:
                return changed, False

        return changed, self._is_success()

    def _insert(self, derivation: LieDerivation, *, source: str) -> InsertResult:
        return self.state.insert(derivation, source=source)

    def _insert_monomial(
        self,
        component: int,
        exponents: tuple[int, ...],
        *,
        source: str,
    ) -> InsertResult:
        self.iterations += 1
        derivation = self.state.make_monomial_derivation(component, exponents)
        return self._insert(derivation, source=source)

    def _build_result(self, phase: str) -> DiffDegreeResult:
        layers = tuple(
            DiffDegreeLayerStatus(
                degree=degree,
                rank=self.state.degree_rank(degree),
                dimension=dim_Wn_coeff_k(self.n, degree),
            )
            for degree in range(0, self.target_degree + 1)
        )
        return DiffDegreeResult(
            success=self._is_success(),
            phase=phase,
            sufficient_condition=self._sufficient_condition(),
            iterations=self.iterations,
            elapsed_s=time.perf_counter() - self._start_time,
            target_rank=self.state.rank_upto(self.target_degree),
            target_dim=self.target_dim,
            total_rank=self.state.total_rank,
            work_degree=self.work_degree,
            partial_axes=tuple(sorted(self.state.partial_axes)),
            monomial_count=len(self.state.monomials),
            degree_one_rank=self.state.degree_rank(1),
            layers=layers,
        )

    def _make_partial_derivative(self, axis: int) -> LieDerivation:
        mapping = {gen: 0 for gen in self.variables}
        mapping[self.variables[axis]] = 1
        return LieDerivation.from_mapping(self.algebra, mapping)

    def _sufficient_condition(self) -> bool:
        if self._sufficient_condition_hit:
            return True

        partials_ok = all(self.state.has_monomial(self._zero(), axis) for axis in range(self.n))
        linear_ok = all(
            self.state.has_monomial(self._unit(axis), component)
            for axis in range(self.n)
            for component in range(self.n)
        )
        quadratic_ok = all(
            self.state.has_monomial(self._unit(axis, power=2), axis)
            for axis in range(self.n)
        )
        self._sufficient_condition_hit = partials_ok and linear_ok and quadratic_ok
        return self._sufficient_condition_hit

    def _is_success(self) -> bool:
        return self._sufficient_condition() or self.state.rank_upto(self.target_degree) >= self.target_dim

    def _mark_change(self, outcome: InsertResult) -> bool:
        return (
            outcome.tracked is not None
            or outcome.monomial_added
            or outcome.partial_added
        )

    def _merge_outcomes(self, left: InsertResult, right: InsertResult) -> InsertResult:
        tracked = right.tracked if right.tracked is not None else left.tracked
        return InsertResult(
            tracked=tracked,
            monomial_added=left.monomial_added or right.monomial_added,
            partial_added=left.partial_added or right.partial_added,
        )

    def _limit_reached(self) -> bool:
        return self.iterations >= self.max_iter

    def _log_phase(self, phase: str, message: str) -> None:
        if not (self.verbose and self.logger.enabled):
            return
        self.logger.rule(phase)
        self.logger.info(message)

    def _maybe_emit_periodic(self, phase: str, *, active: str) -> None:
        if not (self.verbose and self.logger.enabled):
            return

        self._next_log_at = self.logger.maybe_emit_periodic(
            next_log_at=self._next_log_at,
            interval_s=self.log_every_s,
            message_factory=lambda: self._build_status(phase, active=active),
        )

    def _build_status(self, phase: str, *, active: str) -> str:
        lines = [
            "=== Diff-Degree Status ===",
            f"phase = {phase}, iter = {self.iterations}, active = {active}",
            (
                f"rank <= {self.target_degree}: "
                f"{self.state.rank_upto(self.target_degree)}/{self.target_dim}"
            ),
            f"partials = {tuple(sorted(self.state.partial_axes))}",
            f"degree-1 rank = {self.state.degree_rank(1)}",
            f"tracked monomials = {len(self.state.monomials)}",
        ]
        for degree in range(0, self.target_degree + 1):
            lines.append(
                f"deg {degree:2d}: {self.state.degree_rank(degree)}/{dim_Wn_coeff_k(self.n, degree)}"
            )
        return "\n".join(lines)

    def _zero(self) -> tuple[int, ...]:
        return (0,) * self.n

    def _unit(self, axis: int, *, power: int = 1) -> tuple[int, ...]:
        data = [0] * self.n
        data[axis] = power
        return tuple(data)
