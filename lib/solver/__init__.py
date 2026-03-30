"""
Публичный API solver: алгоритмы замыкания по скобке Ли.

Основная точка входа — check_generation(), которая принимает генераторы
в любом формате и диспетчеризует на подходящий solver.
"""

from __future__ import annotations

from importlib import import_module
from typing import Union, TYPE_CHECKING

from .queue_item import QueueItem

if TYPE_CHECKING:
    from lib.core.spec import DerivationSpec
    from lib.core.results import NumericResult, SymbolicResult

__all__ = [
    "LieBasisSolver",
    "QueueItem",
    "check_generation",
]


def check_generation(
    generators,
    *,
    mode: str = "auto",
    n: int | None = None,
    d: int | None = None,
    D_max: int | None = None,
    max_iter: int = 1000,
    verbose: bool = False,
    logger=None,
    structure_constants=None,
):
    """
    Единая точка входа для проверки порождения алгебры Ли W_n.

    Принимает генераторы в любом из форматов:
    - list[DerivationSpec]   — каноническое представление
    - list[LieDerivation]    — объекты SageMath
    - list[Generator]        — dict[int, dict[tuple, float]] (numeric format)

    Args:
        generators:          список генераторов
        mode:                "auto", "numeric" или "symbolic"
        n:                   размерность (нужна для numeric mode с DerivationSpec)
        d:                   максимальная степень (только для numeric mode)
        D_max:               обрезание структурных констант (только для numeric)
        max_iter:            максимальное число итераций (только для symbolic)
        verbose:             вывод прогресса
        logger:              LibraryLogger (только для numeric)
        structure_constants: предвычисленный StructureConstantCache (только для numeric)

    На вход принимает список генераторов, режим и параметры solver.
    На выходе возвращает `NumericResult` или `SymbolicResult`.

    >>> from lib.generators.beldiev import beldiev_specs
    >>> result = check_generation(beldiev_specs(2), mode="numeric", d=2)
    >>> result.success, result.is_full_up_to(2)
    (True, True)
    """
    from lib.core.spec import DerivationSpec
    from lib.core.results import SymbolicResult

    if not generators:
        raise ValueError("Нужен хотя бы один генератор")

    # Определяем тип входных данных
    first = generators[0]
    is_spec = isinstance(first, DerivationSpec)
    is_lie = _is_lie_derivation(first)
    is_generator_dict = isinstance(first, dict)

    # Определяем режим
    effective_mode = mode
    if mode == "auto":
        if d is not None:
            effective_mode = "numeric"
        elif is_lie:
            effective_mode = "symbolic"
        elif is_spec and d is None:
            effective_mode = "symbolic"
        else:
            effective_mode = "numeric"

    if effective_mode == "numeric":
        return _run_numeric(
            generators=generators,
            is_spec=is_spec,
            is_lie=is_lie,
            n=n,
            d=d,
            D_max=D_max,
            verbose=verbose,
            logger=logger,
            structure_constants=structure_constants,
        )
    elif effective_mode == "symbolic":
        return _run_symbolic(
            generators=generators,
            is_spec=is_spec,
            is_lie=is_lie,
            n=n,
            max_iter=max_iter,
        )
    else:
        raise ValueError(f"Неизвестный режим: {mode!r}. Допустимые: 'auto', 'numeric', 'symbolic'")


def _run_numeric(
    generators,
    *,
    is_spec: bool,
    is_lie: bool,
    n: int | None,
    d: int | None,
    D_max: int | None,
    verbose: bool,
    logger,
    structure_constants,
):
    """Запуск numeric solver."""
    from lib.backends.numeric import to_generator
    from lib.numeric.workflows import check_numeric_generation

    # Приводим вход к словарному формату Generator
    if is_spec:
        specs = generators
        effective_n = specs[0].n if n is None else n
        numeric_generators = [to_generator(s) for s in specs]
    elif is_lie:
        raise ValueError(
            "Для numeric mode с LieDerivation нужно сначала преобразовать "
            "в DerivationSpec через lib.backends.sage.from_sage()"
        )
    else:
        # Формат уже совпадает с dict[int, dict[tuple, float]]
        numeric_generators = list(generators)
        if n is None:
            raise ValueError("Для numeric mode с Generator dict необходимо указать n=")
        effective_n = n

    effective_d = d if d is not None else 2

    return check_numeric_generation(
        n=effective_n,
        d=effective_d,
        generators=numeric_generators,
        D_max=D_max,
        verbose=verbose,
        logger=logger,
        structure_constants=structure_constants,
    )


def _run_symbolic(generators, *, is_spec: bool, is_lie: bool, n: int | None, max_iter: int):
    """Запуск symbolic solver."""
    from lib.core.results import SymbolicResult
    from .basis import LieBasisSolver

    if is_spec:
        raise ValueError(
            "Для symbolic mode с DerivationSpec нужно сначала преобразовать "
            "в LieDerivation через lib.backends.sage.to_sage()"
        )

    if not is_lie:
        raise ValueError(
            "Symbolic mode требует list[LieDerivation]. "
            "Для dict-генераторов используйте mode='numeric'."
        )

    solver = LieBasisSolver(generators, max_iter=max_iter)
    success = solver.run()

    return SymbolicResult(
        success=success,
        iterations=solver.iter_count,
        basis_size=len(solver.basis),
        targets_found=dict(solver.targets_found),
    )


def _is_lie_derivation(obj) -> bool:
    """Проверка, является ли объект LieDerivation (без импорта Sage)."""
    return type(obj).__name__ in ("LieDerivation", "QuotientLieDerivation")


def __getattr__(name: str):
    if name == "LieBasisSolver":
        value = getattr(import_module("lib.solver.basis"), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
