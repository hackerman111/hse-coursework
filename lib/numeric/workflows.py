"""
High-level public workflows for numeric Lie-generation experiments.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .criteria import (
    FULL_LIE_GENERATION_DEGREE,
    LieGenerationCriterionResult,
    evaluate_full_lie_generation,
)
from .recipes import make_beldiev_generators
from .results import NumericResult
from .solver import NumericLieGenerationChecker
from .specs import (
    describe_generator,
    describe_random_generator,
    generator_max_degree,
    generator_term_count,
    make_monomial_D,
    make_random_generator,
    parse_generator_spec,
)
from .structure import StructureConstantCache
from .types import Generator
from ..utils import LibraryLogger


@dataclass
class HypothesisTrialResult:
    """One trial of the two-generator numeric hypothesis workflow."""

    index: int
    second_generator: Generator
    numeric_result: NumericResult


@dataclass
class HypothesisRunResult:
    """Aggregate result for several two-generator trials."""

    n: int
    d: int
    trials: list[HypothesisTrialResult] = field(default_factory=list)

    @property
    def successes(self) -> int:
        return sum(1 for trial in self.trials if trial.numeric_result.success)

    @property
    def total_trials(self) -> int:
        return len(self.trials)

    @property
    def all_passed(self) -> bool:
        return self.successes == self.total_trials


@dataclass
class FullLieGenerationCheck:
    """Raw numeric result together with the project-level criterion decision."""

    numeric_result: NumericResult
    criterion_result: LieGenerationCriterionResult

    @property
    def success(self) -> bool:
        return self.criterion_result.satisfied


def resolve_hypothesis_generator(
    n: int,
    spec: str | None,
    a: int,
    b: int,
) -> Generator:
    """
    Разрешить второй генератор гипотезы из явной спецификации или legacy-параметров.

    На вход принимает размерность ``n``, строку ``spec`` и показатели ``a, b``.
    На выходе возвращает generator, который затем можно передать в numeric solver.

    >>> resolve_hypothesis_generator(2, "dy", 1, 2)
    {1: {(0, 0): 1.0}}
    """
    if spec:
        return parse_generator_spec(spec, n)
    return make_monomial_D(n, a, b)


def check_numeric_generation(
    n: int,
    d: int,
    generators: list[Generator],
    *,
    D_max: int | None = None,
    verbose: bool = True,
    log_every_s: float | None = None,
    logger: LibraryLogger | None = None,
    structure_constants: StructureConstantCache | None = None,
) -> NumericResult:
    """
    Запустить основной numeric solver через стабильную публичную обёртку.

    На вход принимает размерность ``n``, глубину проверки ``d`` и список generator-ов с дополнительными параметрами.
    На выходе возвращает ``NumericResult`` со статусом полноты по степеням.

    >>> result = check_numeric_generation(2, 0, [{0: {(0, 0): 1.0}, 1: {(0, 0): 1.0}}], verbose=False)
    >>> result.d, result.n
    (0, 2)
    """
    checker = NumericLieGenerationChecker(
        n,
        d,
        generators,
        D_max=D_max,
        verbose=verbose,
        log_every_s=log_every_s,
        logger=logger,
        structure_constants=structure_constants,
    )
    return checker.run()


def check_full_lie_generation(
    n: int,
    generators: list[Generator],
    *,
    required_degree: int = FULL_LIE_GENERATION_DEGREE,
    D_max: int | None = None,
    verbose: bool = False,
    log_every_s: float | None = None,
    logger: LibraryLogger | None = None,
    structure_constants: StructureConstantCache | None = None,
) -> FullLieGenerationCheck:
    """
    Применить проектный критерий полной порождённости к набору генераторов.

    На вход принимает размерность ``n``, список generator-ов и параметры усечения.
    На выходе возвращает ``FullLieGenerationCheck`` с raw numeric-результатом и итоговым критерием.

    >>> check = check_full_lie_generation(2, [{0: {(0, 0): 1.0}}], verbose=False)
    >>> hasattr(check, "criterion_result")
    True
    """
    numeric_result = check_numeric_generation(
        n=n,
        d=required_degree,
        generators=generators,
        D_max=D_max,
        verbose=verbose,
        log_every_s=log_every_s,
        logger=logger,
        structure_constants=structure_constants,
    )
    criterion_result = evaluate_full_lie_generation(
        numeric_result,
        required_degree=required_degree,
    )
    return FullLieGenerationCheck(
        numeric_result=numeric_result,
        criterion_result=criterion_result,
    )


def run_beldiev(
    n: int,
    d: int,
    *,
    D_max: int | None = None,
    log_every_s: float | None = None,
    logger: LibraryLogger | None = None,
) -> NumericResult:
    """
    Проверить семейство Бельдиева на порождение до степени ``d``.

    На вход принимает размерность ``n``, степень усечения ``d`` и параметры логирования.
    На выходе возвращает ``NumericResult`` для численного эксперимента.

    >>> result = run_beldiev(2, 0)
    >>> result.n, result.d
    (2, 0)
    """
    logger = logger or LibraryLogger()
    logger.banner(
        "Проверка генераторов Бельдиева",
        [f"n = {n}, d = {d}"],
    )

    generators = make_beldiev_generators(n)
    effective_D_max = D_max if D_max is not None else max(d, 4 * n - 1)
    result = check_numeric_generation(
        n=n,
        d=d,
        generators=generators,
        D_max=effective_D_max,
        verbose=True,
        log_every_s=log_every_s,
        logger=logger,
    )
    logger.line()
    logger.line(result.summary(), indent=0)
    return result


def run_hypothesis(
    n: int,
    d: int,
    first_gen: Generator,
    *,
    second_gen: Generator | None = None,
    num_trials: int = 5,
    E_max_degree: int | None = None,
    seed: int = 42,
    log_every_s: float | None = None,
    logger: LibraryLogger | None = None,
) -> HypothesisRunResult:
    """
    Проверить гипотезу двухгенераторной порождённости на серии численных запусков.

    На вход принимает размерность ``n``, степень ``d``, первый generator и параметры выбора второго generator-а.
    На выходе возвращает ``HypothesisRunResult`` со всеми проведёнными испытаниями.

    >>> result = run_hypothesis(2, 0, {0: {(0, 0): 1.0}}, second_gen={1: {(0, 0): 1.0}})
    >>> result.total_trials
    1
    """
    if E_max_degree is None:
        E_max_degree = d + 2
    trial_count = 1 if second_gen is not None else num_trials

    logger = logger or LibraryLogger()
    logger.banner(
        "Проверка гипотезы 1.5-порождаемости",
        [f"n = {n}, d = {d}, trials = {trial_count}"],
    )
    logger.kv("G1", describe_generator(first_gen, n))
    logger.info(
        f"G1 stats: terms={generator_term_count(first_gen)}, "
        f"max_degree={generator_max_degree(first_gen)}"
    )

    if second_gen is None:
        logger.info(
            f"G2 mode = random (trials={num_trials}, max_degree={E_max_degree}, seed={seed})"
        )
    else:
        logger.info("G2 mode = fixed")
        logger.kv("G2", describe_generator(second_gen, n))
        logger.info(
            f"G2 stats: terms={generator_term_count(second_gen)}, "
            f"max_degree={generator_max_degree(second_gen)}"
        )
    logger.line()

    rng = np.random.default_rng(seed)
    run_result = HypothesisRunResult(n=n, d=d)

    if second_gen is not None and num_trials != 1:
        logger.info(f"note: fixed G2 задан, поэтому --trials={num_trials} игнорируется")
        logger.line()

    for trial_index in range(trial_count):
        logger.rule(label=f"Trial {trial_index + 1}/{trial_count}")
        current_second = second_gen
        if current_second is None:
            current_second = make_random_generator(
                n,
                E_max_degree,
                rng=rng,
                sparsity=0.3,
            )
            logger.info(f"G2 random: {describe_random_generator(current_second, n)}")

        current_max_degree = max(
            generator_max_degree(first_gen),
            generator_max_degree(current_second),
            E_max_degree if second_gen is None else -1,
        )
        numeric_result = check_numeric_generation(
            n=n,
            d=d,
            generators=[first_gen, current_second],
            D_max=max(d, current_max_degree + 1),
            verbose=True,
            log_every_s=log_every_s,
            logger=logger,
        )
        run_result.trials.append(
            HypothesisTrialResult(
                index=trial_index,
                second_generator=current_second,
                numeric_result=numeric_result,
            )
        )

        if numeric_result.success:
            logger.info(
                f"PASS ({numeric_result.elapsed_s:.3f}s, {numeric_result.iterations} iter)"
            )
        else:
            logger.info(
                f"FAIL ({numeric_result.elapsed_s:.3f}s, {numeric_result.iterations} iter)"
            )
            for status in numeric_result.missing_degrees():
                logger.info(str(status), indent=4)
        logger.line()

    logger.info(f"Итого: {run_result.successes}/{run_result.total_trials} успешных")
    if run_result.all_passed:
        logger.info("Гипотеза ПОДТВЕРЖДЕНА (численно, для данных E)")
    elif run_result.successes > 0:
        logger.info("Частичное подтверждение")
    else:
        logger.info("Все попытки ПРОВАЛЕНЫ - кандидат на контрпример?")
        logger.info(
            "(Попробуйте увеличить d, E_max_degree или число попыток)",
            indent=4,
        )

    return run_result
