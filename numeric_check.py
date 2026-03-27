#!/usr/bin/env python3
"""
numeric_check.py — CLI для численной проверки порождения W_n.

Режимы:
  1. beldiev   — проверка генераторов Бельдиева;
  2. hypo      — проверка пары генераторов G1, G2 до степени d;
  3. criterion — проверка критерия "вся алгебра Ли порождена, если
                 порождены все дифференцирования до степени 2 включительно";
  4. dims      — сверка таблицы размерностей из Numeric.md.

Примеры:
  python numeric_check.py --mode beldiev --n 2 --d 7
  python numeric_check.py --mode hypo --n 2 --d 8 --g1 "x^5*y^5*dx" --g2 "dx + x*dy"
  python numeric_check.py --mode criterion --n 2 --g1 "x^5*y^5*dx" --g2 "dx + x*dy"
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib.numeric import (
    FULL_LIE_GENERATION_DEGREE,
    check_full_lie_generation,
    describe_generator,
    parse_generator_spec,
    resolve_hypothesis_generator,
    run_beldiev,
    run_hypothesis,
)
from lib.numeric.reference import check_dimension_examples
from lib.utils import LibraryLogger


class RichHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    """Help formatter with preserved examples and visible defaults."""


def _print_criterion_result(
    logger: LibraryLogger,
    *,
    g1,
    g2,
    n: int,
    criterion_degree: int,
    check_result,
) -> None:
    logger.banner(
        "Проверка критерия полной порождаемости",
        [f"n = {n}, criterion_degree = {criterion_degree}"],
    )
    logger.kv("G1", describe_generator(g1, n))
    logger.kv("G2", describe_generator(g2, n))
    logger.line()
    logger.line(check_result.numeric_result.summary(), indent=0)
    logger.line()
    logger.info(check_result.criterion_result.summary(), indent=0)
    logger.info(check_result.criterion_result.reason, indent=0)


def main() -> None:
    description = (
        "Численная проверка порождения W_n.\n\n"
        "Режимы:\n"
        "  beldiev   : эталонная проверка известных генераторов Бельдиева.\n"
        "  hypo      : проверка пары генераторов G1, G2 до степени d.\n"
        "  criterion : проверка критерия полной порождаемости по степеням <= 2.\n"
        "  dims      : сверка таблицы размерностей из Numeric.md.\n\n"
        "Синтаксис generator-spec:\n"
        "  сумма членов вида [coeff*]monomial*differential\n"
        "  примеры: x^5*y^5*dx, dx + x*dy, 2*z0^3*d1 - 0.5*z1*dz0\n"
        "  допустимые переменные: x,y,z,u,v,w и z0,z1,...\n"
        "  допустимые производные: dx,dy,..., d0,d1,..., dz0,dz1,..."
    )
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=RichHelpFormatter,
        epilog=__doc__,
    )
    common = parser.add_argument_group("Общие параметры")
    common.add_argument(
        "--mode",
        choices=["beldiev", "hypo", "criterion", "dims"],
        default="beldiev",
        help="Режим запуска",
    )
    common.add_argument("--n", type=int, default=2, help="Число переменных")
    common.add_argument(
        "--d",
        type=int,
        default=7,
        help="Проверять до степени d включительно",
    )
    common.add_argument(
        "--Dmax",
        type=int,
        default=None,
        help="Максимальное усечение D_max",
    )
    common.add_argument(
        "--log-every",
        type=float,
        default=5.0,
        help="Раз в сколько секунд печатать полную сводку по степеням; 0 отключает",
    )
    common.add_argument(
        "--criterion-degree",
        type=int,
        default=FULL_LIE_GENERATION_DEGREE,
        help="Степенной порог project-criterion для полной порождаемости",
    )

    hypo_group = parser.add_argument_group("Параметры режимов hypo / criterion")
    hypo_group.add_argument(
        "--g1",
        type=str,
        default=None,
        help="Первый генератор в формате generator-spec",
    )
    hypo_group.add_argument(
        "--g2",
        type=str,
        default=None,
        help="Второй генератор в формате generator-spec",
    )
    hypo_group.add_argument(
        "--Da",
        type=int,
        default=5,
        help="Legacy shorthand для G1: экспонента z0^a при отсутствии --g1",
    )
    hypo_group.add_argument(
        "--Db",
        type=int,
        default=5,
        help="Legacy shorthand для G1: экспонента z1^b при отсутствии --g1",
    )
    hypo_group.add_argument(
        "--trials",
        type=int,
        default=5,
        help="Число случайных G2; используется только когда --g2 не задан",
    )
    hypo_group.add_argument(
        "--E-max-deg",
        type=int,
        default=None,
        help="Максимальная степень случайного G2",
    )
    hypo_group.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Seed для случайного G2",
    )

    args = parser.parse_args()

    if args.mode == "dims":
        ok = check_dimension_examples()
        sys.exit(0 if ok else 1)

    if args.mode == "beldiev":
        result = run_beldiev(
            args.n,
            args.d,
            D_max=args.Dmax,
            log_every_s=args.log_every,
        )
        sys.exit(0 if result.success else 1)

    try:
        first_gen = resolve_hypothesis_generator(args.n, args.g1, args.Da, args.Db)
        second_gen = parse_generator_spec(args.g2, args.n) if args.g2 else None
    except ValueError as exc:
        parser.error(str(exc))

    if args.mode == "hypo":
        run_result = run_hypothesis(
            n=args.n,
            d=args.d,
            first_gen=first_gen,
            second_gen=second_gen,
            num_trials=args.trials,
            E_max_degree=args.E_max_deg,
            seed=args.seed,
            log_every_s=args.log_every,
        )
        sys.exit(0 if run_result.successes > 0 else 1)

    if second_gen is None:
        parser.error("--mode criterion требует заданный --g2")

    logger = LibraryLogger()
    check_result = check_full_lie_generation(
        n=args.n,
        generators=[first_gen, second_gen],
        required_degree=args.criterion_degree,
        D_max=args.Dmax,
        verbose=True,
        log_every_s=args.log_every,
        logger=logger,
    )
    _print_criterion_result(
        logger,
        g1=first_gen,
        g2=second_gen,
        n=args.n,
        criterion_degree=args.criterion_degree,
        check_result=check_result,
    )
    sys.exit(0 if check_result.success else 1)


if __name__ == "__main__":
    main()
