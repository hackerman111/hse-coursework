#!/usr/bin/env python3
"""
diff_deg_check.py — символьная проверка порождения W_n по Diff_deg_D.md.

Режимы:
  1. beldiev  — проверка известных генераторов Бельдиева
  2. hypo     — проверка пары G1, G2; если G2 не задан, берётся случайный
  3. dims     — сверка размерностей в градуировке degree(coefficients)
"""

from __future__ import annotations

import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sage.all import PolynomialRing, QQ

from lib.generator import get_Andristy, get_Beldiev
from lib.solver import DiffDegreeLieGenerationChecker, dim_Wn_coeff_leq_d
from lib.utils import (
    LibraryLogger,
    describe_generator_mapping,
    derivation_to_mapping,
    generator_coeff_degree,
    generator_term_count,
    mapping_to_derivation,
    parse_generator_spec,
)


class RichHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    """Help formatter with preserved examples and visible defaults."""


def make_polynomial_ring(n: int):
    names = [f"z{i}" for i in range(n)]
    return PolynomialRing(QQ, names)


def make_monomial_derivation(algebra, a: int, b: int, target_var: int = 0):
    n = len(algebra.gens())
    exponents = [0] * n
    if n >= 1:
        exponents[0] = a
    if n >= 2:
        exponents[1] = b
    return mapping_to_derivation(algebra, {target_var: {tuple(exponents): 1}})


def make_random_derivation(algebra, max_degree: int, rng: random.Random):
    n = len(algebra.gens())
    terms = max(1, 2 * n)
    mapping: dict[int, dict[tuple[int, ...], int]] = {}

    for _ in range(terms):
        component = rng.randrange(n)
        degree = rng.randrange(max_degree + 1)
        exponents = [0] * n
        for _ in range(degree):
            exponents[rng.randrange(n)] += 1
        coefficient = rng.choice([-2, -1, 1, 2])
        mapping.setdefault(component, {})
        key = tuple(exponents)
        mapping[component][key] = mapping[component].get(key, 0) + coefficient
        if mapping[component][key] == 0:
            del mapping[component][key]
        if not mapping[component]:
            del mapping[component]

    if not mapping:
        mapping = {0: {(0,) * n: 1}}
    return mapping_to_derivation(algebra, mapping)


def run_beldiev(
    n: int,
    d: int,
    *,
    work_degree: int | None = None,
    log_every_s: float | None = None,
    logger: LibraryLogger | None = None,
):
    logger = logger or LibraryLogger()
    algebra = make_polynomial_ring(n)
    generators = get_Beldiev(algebra)

    logger.banner(
        "Diff-Degree Symbolic Check",
        [f"mode = beldiev, n = {n}, d = {d}"],
    )
    logger.kv(f"dim W_{n}[coeff<={d}]", dim_Wn_coeff_leq_d(n, d))
    logger.line()

    checker = DiffDegreeLieGenerationChecker(
        generators,
        d,
        work_degree=work_degree,
        logger=logger,
        verbose=True,
        log_every_s=log_every_s,
    )
    result = checker.run()

    logger.line()
    logger.line(result.summary(), indent=0)
    return result


def run_andrist(
    n: int,
    d: int,
    *,
    work_degree: int | None = None,
    log_every_s: float | None = None,
    logger: LibraryLogger | None = None,
):
    logger = logger or LibraryLogger()
    algebra = make_polynomial_ring(n)
    generators = get_Andristy(algebra)

    logger.banner(
        "Diff-Degree Symbolic Check",
        [f"mode = andrist, n = {n}, d = {d}"],
    )
    logger.kv(f"dim W_{n}[coeff<={d}]", dim_Wn_coeff_leq_d(n, d))
    logger.line()

    checker = DiffDegreeLieGenerationChecker(
        generators,
        d,
        work_degree=work_degree,
        logger=logger,
        verbose=True,
        log_every_s=log_every_s,
    )
    result = checker.run()

    logger.line()
    logger.line(result.summary(), indent=0)
    return result


def run_hypothesis(
    n: int,
    d: int,
    first_gen,
    *,
    second_gen=None,
    num_trials: int = 5,
    second_max_degree: int | None = None,
    seed: int = 42,
    work_degree: int | None = None,
    log_every_s: float | None = None,
    logger: LibraryLogger | None = None,
):
    logger = logger or LibraryLogger()
    rng = random.Random(seed)
    trial_count = 1 if second_gen is not None else num_trials

    logger.banner(
        "Diff-Degree Symbolic Check",
        [f"mode = hypo, n = {n}, d = {d}, trials = {trial_count}"],
    )
    logger.kv("G1", describe_generator_mapping(derivation_to_mapping(first_gen), n))
    if second_gen is None:
        logger.info(
            f"G2 mode = random (trials={num_trials}, max_degree={second_max_degree}, seed={seed})"
        )
    else:
        logger.info("G2 mode = fixed")
        logger.kv("G2", describe_generator_mapping(derivation_to_mapping(second_gen), n))
    logger.kv(f"dim W_{n}[coeff<={d}]", dim_Wn_coeff_leq_d(n, d))
    logger.line()

    successes = 0
    for trial in range(trial_count):
        logger.rule(f"Trial {trial + 1}/{trial_count}")
        current_second = second_gen
        if current_second is None:
            current_second = make_random_derivation(
                first_gen.codomain(),
                second_max_degree if second_max_degree is not None else d,
                rng,
            )
            mapping = derivation_to_mapping(current_second)
            logger.info(
                "G2 random: "
                f"terms={generator_term_count(mapping)}, "
                f"max_degree={generator_coeff_degree(mapping)}, "
                f"sample={describe_generator_mapping(mapping, n, max_terms=4)}"
            )

        checker = DiffDegreeLieGenerationChecker(
            [first_gen, current_second],
            d,
            work_degree=work_degree,
            logger=logger,
            verbose=True,
            log_every_s=log_every_s,
        )
        result = checker.run()
        if result.success:
            successes += 1
            logger.info(f"PASS ({result.elapsed_s:.3f}s, {result.iterations} iter)")
        else:
            logger.info(f"FAIL ({result.elapsed_s:.3f}s, {result.iterations} iter)")
            for layer in result.layers:
                if not layer.is_full:
                    logger.info(str(layer), indent=4)
        logger.line()

    logger.info(f"Итого: {successes}/{trial_count} успешных")
    if successes == trial_count:
        logger.info("Гипотеза ПОДТВЕРЖДЕНА (символьно)")
    elif successes > 0:
        logger.info("Частичное подтверждение")
    else:
        logger.info("Все попытки ПРОВАЛЕНЫ - кандидат на контрпример?")

    return successes


def check_dimensions(logger: LibraryLogger | None = None):
    logger = logger or LibraryLogger()
    logger.line("Проверка таблицы размерностей для degree(coefficients) <= d:", indent=0)
    expected = [
        (2, 5, 42),
        (2, 10, 132),
        (3, 5, 168),
        (3, 10, 858),
        (5, 5, 1260),
        (5, 10, 15015),
    ]
    all_ok = True
    for n, d, dim in expected:
        actual = dim_Wn_coeff_leq_d(n, d)
        ok = "✓" if actual == dim else "✗"
        if actual != dim:
            all_ok = False
        logger.info(f"{ok} dim W_{n}[coeff<={d}] = {actual} (ожидалось {dim})")
    logger.line()
    return all_ok


def _resolve_hypothesis_generator(algebra, spec: str | None, a: int, b: int):
    if spec:
        return mapping_to_derivation(algebra, parse_generator_spec(spec, len(algebra.gens())))
    return make_monomial_derivation(algebra, a, b)


def main():
    description = (
        "Символьная проверка порождения W_n по Diff_deg_D.md.\n\n"
        "Режимы:\n"
        "  beldiev : эталонная проверка генераторов Бельдиева.\n"
        "  andrist : эталонная проверка генераторов Андриста.\n"
        "  hypo    : проверка пары G1, G2. Если G2 не задан,\n"
        "            он генерируется случайно на каждом trial.\n"
        "  dims    : сверка размерностей в degree(coefficients).\n\n"
        "Синтаксис generator-spec:\n"
        "  сумма членов вида [coeff*]monomial*differential\n"
        "  примеры: x^5*y^5*dx, dx + x*dy, 2*z0^3*d1 - z1*dz0\n"
    )
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=RichHelpFormatter,
        epilog=__doc__,
    )

    common = parser.add_argument_group("Общие параметры")
    common.add_argument(
        "--mode",
        choices=["beldiev", "andrist", "hypo", "dims"],
        default="beldiev",
        help="Режим запуска",
    )
    common.add_argument("--n", type=int, default=2, help="Число переменных")
    common.add_argument(
        "--d",
        type=int,
        default=5,
        help="Проверять порождение до degree(coefficients) <= d",
    )
    common.add_argument("--Dmax", type=int, default=None, help="Рабочее усечение D")
    common.add_argument(
        "--log-every",
        type=float,
        default=5.0,
        help="Раз в сколько секунд печатать сводку; 0 отключает периодические сводки",
    )

    hypo_group = parser.add_argument_group("Параметры режима hypo")
    hypo_group.add_argument("--g1", type=str, default=None, help="Первый генератор")
    hypo_group.add_argument("--g2", type=str, default=None, help="Второй генератор")
    hypo_group.add_argument(
        "--Da",
        type=int,
        default=5,
        help="Legacy shorthand для G1: степень при z0, если --g1 не задан",
    )
    hypo_group.add_argument(
        "--Db",
        type=int,
        default=5,
        help="Legacy shorthand для G1: степень при z1, если --g1 не задан",
    )
    hypo_group.add_argument(
        "--trials",
        type=int,
        default=5,
        help="Число случайных G2; используется только без --g2",
    )
    hypo_group.add_argument(
        "--E-max-deg",
        type=int,
        default=None,
        help="Максимальная степень коэффициентов случайного G2",
    )
    hypo_group.add_argument("--seed", type=int, default=42, help="Seed для случайного G2")

    args = parser.parse_args()

    if args.mode == "dims":
        sys.exit(0 if check_dimensions() else 1)

    if args.mode == "beldiev":
        result = run_beldiev(
            args.n,
            args.d,
            work_degree=args.Dmax,
            log_every_s=args.log_every,
        )
        sys.exit(0 if result.success else 1)

    if args.mode == "andrist":
        result = run_andrist(
            args.n,
            args.d,
            work_degree=args.Dmax,
            log_every_s=args.log_every,
        )
        sys.exit(0 if result.success else 1)

    algebra = make_polynomial_ring(args.n)
    try:
        first_gen = _resolve_hypothesis_generator(algebra, args.g1, args.Da, args.Db)
        second_gen = (
            mapping_to_derivation(algebra, parse_generator_spec(args.g2, args.n))
            if args.g2
            else None
        )
    except ValueError as exc:
        parser.error(str(exc))

    successes = run_hypothesis(
        n=args.n,
        d=args.d,
        first_gen=first_gen,
        second_gen=second_gen,
        num_trials=args.trials,
        second_max_degree=args.E_max_deg,
        seed=args.seed,
        work_degree=args.Dmax,
        log_every_s=args.log_every,
    )
    sys.exit(0 if successes > 0 else 1)


if __name__ == "__main__":
    main()
