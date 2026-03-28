#!/usr/bin/env python3
"""
symbolic_check.py — CLI для символьной проверки порождения W_n (требует SageMath).

Режимы:
  1. hypo      — проверка пары генераторов G1, G2 через LieBasisSolver;
  2. criterion — проверка: порождены ли все дифференцирования степени <= 2
                 (dx_i, x_j*dx_i, x_j*x_k*dx_i для всех допустимых индексов).

Примеры:
  sage -python symbolic_check.py --mode hypo --n 2 --g1 "dx" --g2 "x*dy"
  sage -python symbolic_check.py --mode hypo --n 2 --g1 "x^5*y^5*dx" --seed 42
  sage -python symbolic_check.py --mode criterion --n 2 --g1 "dx" --g2 "x*dy"
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sage.all import PolynomialRing, QQ

from lib import spec_from_string, spec_to_string, random_spec
from lib.backends.sage import to_sage
from lib.core.results import SymbolicResult
from lib.solver.basis import LieBasisSolver
from lib.utils import LibraryLogger


class RichHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    """Help formatter with preserved examples and visible defaults."""


def _make_algebra(n: int):
    """Создать кольцо многочленов K[z0, ..., z_{n-1}]."""
    var_names = ",".join(f"z{i}" for i in range(n))
    return PolynomialRing(QQ, var_names, order="degrevlex")


def _spec_to_lie(spec, algebra):
    """DerivationSpec → LieDerivation."""
    return to_sage(spec, algebra)


def _print_result(
    logger: LibraryLogger,
    result: SymbolicResult,
    *,
    g1_str: str,
    g2_str: str,
    n: int,
    mode: str,
) -> None:
    """Вывод результатов символьной проверки."""
    logger.banner(
        f"Символьная проверка порождения W_{n}",
        [f"Режим: {mode}"],
    )
    logger.kv("G1", g1_str)
    logger.kv("G2", g2_str)
    logger.line()
    logger.line(result.summary(), indent=0)


def run_hypo(
    *,
    n: int,
    g1_spec,
    g2_spec,
    max_iter: int,
    logger: LibraryLogger,
) -> SymbolicResult:
    """Запуск режима hypo: полная символьная проверка пары генераторов."""
    algebra = _make_algebra(n)

    g1 = _spec_to_lie(g1_spec, algebra)
    g2 = _spec_to_lie(g2_spec, algebra)

    solver = LieBasisSolver([g1, g2], max_iter=max_iter)
    success = solver.run()

    return SymbolicResult(
        success=success,
        iterations=solver.iter_count,
        basis_size=len(solver.basis),
        targets_found=dict(solver.targets_found),
    )


def main() -> None:
    description = (
        "Символьная проверка порождения W_n (требует SageMath).\n\n"
        "Режимы:\n"
        "  hypo      : проверка пары генераторов G1, G2 через LieBasisSolver.\n"
        "  criterion : проверка порождения всех дифференцирований степени <= 2.\n\n"
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
        choices=["hypo", "criterion"],
        default="hypo",
        help="Режим запуска",
    )
    common.add_argument("--n", type=int, default=2, help="Число переменных")
    common.add_argument(
        "--max-iter",
        type=int,
        default=1000,
        help="Максимальное число итераций солвера",
    )

    gen_group = parser.add_argument_group("Параметры генераторов")
    gen_group.add_argument(
        "--g1",
        type=str,
        required=True,
        help="Первый генератор в формате generator-spec",
    )
    gen_group.add_argument(
        "--g2",
        type=str,
        default=None,
        help="Второй генератор в формате generator-spec (если не задан — случайный)",
    )
    gen_group.add_argument(
        "--trials",
        type=int,
        default=5,
        help="Число случайных G2 (когда --g2 не задан)",
    )
    gen_group.add_argument(
        "--E-max-deg",
        type=int,
        default=3,
        help="Максимальная степень случайного G2",
    )
    gen_group.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Seed для случайного G2",
    )

    args = parser.parse_args()
    logger = LibraryLogger()

    # Парсим G1
    try:
        g1_spec = spec_from_string(args.g1, args.n)
    except ValueError as exc:
        parser.error(f"Ошибка разбора --g1: {exc}")

    # G2: явный или случайный
    if args.g2 is not None:
        try:
            g2_spec = spec_from_string(args.g2, args.n)
        except ValueError as exc:
            parser.error(f"Ошибка разбора --g2: {exc}")

        g1_str = spec_to_string(g1_spec)
        g2_str = spec_to_string(g2_spec)

        if args.mode == "hypo":
            result = run_hypo(
                n=args.n,
                g1_spec=g1_spec,
                g2_spec=g2_spec,
                max_iter=args.max_iter,
                logger=logger,
            )
            _print_result(logger, result, g1_str=g1_str, g2_str=g2_str,
                          n=args.n, mode="hypo")
            sys.exit(0 if result.success else 1)
        else:
            # criterion mode — will be added in Task 3
            parser.error("Режим criterion ещё не реализован (будет в Task 3)")

    else:
        # Случайные G2
        g1_str = spec_to_string(g1_spec)
        successes = 0

        for trial in range(args.trials):
            trial_seed = args.seed + trial
            g2_spec = random_spec(
                args.n,
                max_degree=args.E_max_deg,
                seed=trial_seed,
            )
            g2_str = spec_to_string(g2_spec)

            logger.info(f"--- Попытка {trial + 1}/{args.trials} (seed={trial_seed}) ---")
            logger.kv("G2", g2_str)

            if args.mode == "hypo":
                result = run_hypo(
                    n=args.n,
                    g1_spec=g1_spec,
                    g2_spec=g2_spec,
                    max_iter=args.max_iter,
                    logger=logger,
                )
                _print_result(logger, result, g1_str=g1_str, g2_str=g2_str,
                              n=args.n, mode="hypo")
                if result.success:
                    successes += 1
            else:
                # criterion mode — will be added in Task 3
                parser.error("Режим criterion ещё не реализован (будет в Task 3)")

            logger.line()

        logger.info(f"Итого: {successes}/{args.trials} успешных попыток")
        sys.exit(0 if successes > 0 else 1)


if __name__ == "__main__":
    main()
