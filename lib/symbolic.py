"""
High-level public workflows for symbolic Lie-generation experiments.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from lib.core.results import SymbolicResult
from lib.core.spec import DerivationSpec
from lib.solver.basis import LieBasisSolver
from lib.solver.criterion import (
    enumerate_targets_up_to_degree2,
    format_target,
    make_degree2_stop_condition,
)

__all__ = [
    "SymbolicCriterionResult",
    "check_symbolic_degree2_criterion",
    "check_symbolic_generation",
    "enumerate_targets_up_to_degree2",
    "format_target",
    "make_degree2_stop_condition",
    "make_symbolic_algebra",
    "spec_to_symbolic_derivation",
]


@dataclass(frozen=True)
class SymbolicCriterionResult:
    """Structured result for the degree-2 symbolic stopping criterion."""

    success: bool
    iterations: int = 0
    basis_size: int = 0
    found_targets: tuple[str, ...] = field(default_factory=tuple)
    missing_targets: tuple[str, ...] = field(default_factory=tuple)

    @property
    def total_targets(self) -> int:
        return len(self.found_targets) + len(self.missing_targets)

    def summary(self) -> str:
        lines = [
            "=== Symbolic Degree-2 Criterion Check ===",
            f"Итераций: {self.iterations}, Базис: {self.basis_size} элементов",
            f"Найдено целей: {len(self.found_targets)}/{self.total_targets}",
            (
                "Результат: PASS — все дифференцирования степени <= 2 найдены"
                if self.success
                else "Результат: FAIL — не все дифференцирования степени <= 2 найдены"
            ),
        ]
        if self.missing_targets:
            lines.append("Не найдены:")
            lines.extend(f"  {target}" for target in self.missing_targets)
        return "\n".join(lines)


def make_symbolic_algebra(n: int) -> Any:
    """
    Построить полиномиальную алгебру Sage для symbolic solver-а.

    На вход принимает число переменных ``n``.
    На выходе возвращает ``PolynomialRing(QQ, z0,...,z{n-1})`` с порядком ``degrevlex``.

    >>> ring = make_symbolic_algebra(2)
    >>> tuple(str(gen) for gen in ring.gens())
    ('z0', 'z1')
    """
    _ensure_sage_runtime_dir()
    from sage.all import PolynomialRing, QQ

    var_names = ",".join(f"z{i}" for i in range(n))
    return PolynomialRing(QQ, var_names, order="degrevlex")


def spec_to_symbolic_derivation(spec: DerivationSpec, algebra: Any) -> Any:
    """
    Конвертировать ``DerivationSpec`` в ``LieDerivation`` для symbolic solver-а.

    На вход принимает ``DerivationSpec`` и полиномиальную алгебру Sage ``algebra``.
    На выходе возвращает объект ``LieDerivation``, готовый для ``LieBasisSolver``.

    >>> ring = make_symbolic_algebra(2)
    >>> derivation = spec_to_symbolic_derivation(DerivationSpec(n=2, terms={0: {(0, 0): 1.0}}), ring)
    >>> type(derivation).__name__
    'LieDerivation'
    """
    from lib.backends.sage import to_sage

    return to_sage(spec, algebra)


def check_symbolic_generation(
    generators: Iterable[Any],
    *,
    n: int | None = None,
    algebra: Any | None = None,
    max_iter: int = 1000,
) -> SymbolicResult:
    """
    Запустить symbolic solver на наборе ``DerivationSpec`` или ``LieDerivation``.

    На вход принимает итерируемый набор symbolic-генераторов, а также опциональные ``n``, ``algebra`` и ``max_iter``.
    На выходе возвращает ``SymbolicResult`` с числом итераций, размером базиса и покрытием базовых производных.

    >>> specs = [DerivationSpec(n=2, terms={0: {(0, 0): 1.0}}), DerivationSpec(n=2, terms={1: {(0, 0): 1.0}})]
    >>> result = check_symbolic_generation(specs, max_iter=10)
    >>> isinstance(result.success, bool)
    True
    """
    symbolic_generators, _ = _resolve_symbolic_generators(
        generators,
        n=n,
        algebra=algebra,
    )
    solver = LieBasisSolver(symbolic_generators, max_iter=max_iter)
    success = solver.run()

    return SymbolicResult(
        success=success,
        iterations=solver.iter_count,
        basis_size=len(solver.basis),
        targets_found=dict(solver.targets_found),
    )


def check_symbolic_degree2_criterion(
    generators: Iterable[Any],
    *,
    n: int | None = None,
    algebra: Any | None = None,
    max_iter: int = 1000,
) -> SymbolicCriterionResult:
    """
    Проверить symbolic-критерий покрытия всех дифференцирований степени ``<= 2``.

    На вход принимает итерируемый набор ``DerivationSpec`` или ``LieDerivation`` и параметры symbolic solver-а.
    На выходе возвращает ``SymbolicCriterionResult`` со списками найденных и отсутствующих целевых направлений.

    >>> specs = [DerivationSpec(n=2, terms={0: {(0, 0): 1.0}}), DerivationSpec(n=2, terms={1: {(0, 0): 1.0}})]
    >>> result = check_symbolic_degree2_criterion(specs, max_iter=10)
    >>> isinstance(result.success, bool)
    True
    """
    symbolic_generators, resolved_algebra = _resolve_symbolic_generators(
        generators,
        n=n,
        algebra=algebra,
    )
    stop_condition, targets_status = make_degree2_stop_condition(resolved_algebra)
    solver = LieBasisSolver(
        symbolic_generators,
        max_iter=max_iter,
        stop_condition=stop_condition,
    )
    solver.run()

    for target_key in targets_status:
        if not targets_status[target_key] and target_key in solver.basis:
            targets_status[target_key] = True

    found_targets = tuple(
        format_target(axis, monomial, resolved_algebra)
        for (axis, monomial), is_found in targets_status.items()
        if is_found
    )
    missing_targets = tuple(
        format_target(axis, monomial, resolved_algebra)
        for (axis, monomial), is_found in targets_status.items()
        if not is_found
    )

    return SymbolicCriterionResult(
        success=not missing_targets,
        iterations=solver.iter_count,
        basis_size=len(solver.basis),
        found_targets=found_targets,
        missing_targets=missing_targets,
    )


def _resolve_symbolic_generators(
    generators: Iterable[Any],
    *,
    n: int | None,
    algebra: Any | None,
) -> tuple[list[Any], Any]:
    generator_list = list(generators)
    if not generator_list:
        raise ValueError("Нужен хотя бы один symbolic generator")

    first = generator_list[0]
    if isinstance(first, DerivationSpec):
        if not all(isinstance(item, DerivationSpec) for item in generator_list):
            raise ValueError(
                "Все symbolic generators должны быть либо DerivationSpec, "
                "либо LieDerivation"
            )
        effective_n = first.n if n is None else n
        if effective_n != first.n:
            raise ValueError(
                f"n={effective_n} не совпадает с размерностью spec.n={first.n}"
            )
        resolved_algebra = algebra if algebra is not None else make_symbolic_algebra(effective_n)
        if hasattr(resolved_algebra, "gens") and len(resolved_algebra.gens()) != effective_n:
            raise ValueError(
                "Размерность symbolic algebra не совпадает с размерностью генераторов"
            )
        return [
            spec_to_symbolic_derivation(spec, resolved_algebra)
            for spec in generator_list
        ], resolved_algebra

    if _is_lie_derivation(first):
        if not all(_is_lie_derivation(item) for item in generator_list):
            raise ValueError(
                "Все symbolic generators должны быть либо DerivationSpec, "
                "либо LieDerivation"
            )
        resolved_algebra = algebra if algebra is not None else first.codomain()
        return generator_list, resolved_algebra

    raise ValueError(
        "Symbolic API принимает list[DerivationSpec] или list[LieDerivation]"
    )


def _is_lie_derivation(obj: Any) -> bool:
    return type(obj).__name__ in ("LieDerivation", "QuotientLieDerivation")


def _ensure_sage_runtime_dir() -> None:
    current = os.environ.get("DOT_SAGE")
    if current:
        current_path = Path(current)
        if current_path.exists() and os.access(current_path, os.W_OK):
            return

    fallback = Path("/tmp") / f"codex-sage-{os.getuid()}"
    fallback.mkdir(parents=True, exist_ok=True)
    os.environ["DOT_SAGE"] = str(fallback)
