#!/usr/bin/env python3
"""
numeric_check.py — быстрая проверка гипотез о порождении W_n.

Режимы:
  1. beldiev  — проверка генераторов Бельдиева (известный результат)
  2. hypo     — проверка гипотезы 1.5-порождаемости:
                фиксируем D, выбираем случайный E, проверяем Lie(D, E) ⊇ W_n^{(≤d)}

Примеры:
  python numeric_check.py --mode beldiev --n 2 --d 7
  python numeric_check.py --mode hypo --n 2 --d 10 --D "x^5*y^5*dx"
  python numeric_check.py --mode hypo --n 2 --d 8 --Da 5 --Db 5
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

# Добавляем корень проекта в путь
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib.numeric.solver import (
    NumericLieGenerationChecker,
    make_beldiev_generators,
    make_random_generator,
)
from lib.numeric.indexing import dim_Wn_k, dim_Wn_leq_d


# ──────────────────────────────────────────────────────────────────────────────
# Вспомогательные функции для создания генераторов
# ──────────────────────────────────────────────────────────────────────────────

def make_monomial_D(n: int, a: int, b: int, target_var: int = 0):
    """
    Создаёт генератор D = x^a * y^b * ∂_{z_{target_var}} для n = 2.
    Для n > 2 — D = z_0^a * z_1^b * ∂_{z_{target_var}}.
    """
    alpha = [0] * n
    if n >= 1:
        alpha[0] = a
    if n >= 2:
        alpha[1] = b
    return {target_var: {tuple(alpha): 1.0}}


def make_general_monomial_D(n: int, exponents: list, target_var: int = 0):
    """
    Создаёт генератор D = z_0^{e_0} · z_1^{e_1} · … · z_{n-1}^{e_{n-1}} · ∂_{z_{target_var}}.
    """
    alpha = [0] * n
    for i, e in enumerate(exponents):
        if i < n:
            alpha[i] = e
    return {target_var: {tuple(alpha): 1.0}}


def make_random_E(n: int, max_degree: int, seed: int | None = None):
    """Случайный генератор E достаточно общего вида."""
    rng = np.random.default_rng(seed)
    return make_random_generator(n, max_degree, rng=rng, sparsity=0.3)


# ──────────────────────────────────────────────────────────────────────────────
# Режим 1: проверка генераторов Бельдиева
# ──────────────────────────────────────────────────────────────────────────────

def run_beldiev(n: int, d: int, D_max: int | None = None):
    """Проверить генераторы Бельдиева для W_n до степени d."""
    print(f"╔══════════════════════════════════════════╗")
    print(f"║  Проверка генераторов Бельдиева          ║")
    print(f"║  n = {n}, d = {d}                           ║")
    print(f"╚══════════════════════════════════════════╝")
    print()

    generators = make_beldiev_generators(n)

    # Показать размерности
    total_dim = dim_Wn_leq_d(n, d)
    print(f"  dim W_{n}^(≤{d}) = {total_dim}")
    print()

    if D_max is None:
        D_max = max(d, 4 * n - 1)

    checker = NumericLieGenerationChecker(n, d, generators, D_max=D_max, verbose=True)
    result = checker.run()

    print()
    print(result.summary())
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Режим 2: проверка гипотезы 1.5-порождаемости
# ──────────────────────────────────────────────────────────────────────────────

def run_hypothesis(
    n: int,
    d: int,
    D_gen,
    num_trials: int = 5,
    E_max_degree: int | None = None,
    seed: int = 42,
):
    """
    Для фиксированного D проверить, порождает ли Lie(D, E) = W_n^{(≤d)}
    для нескольких случайных E.
    """
    if E_max_degree is None:
        E_max_degree = d + 2

    print(f"╔══════════════════════════════════════════╗")
    print(f"║  Проверка гипотезы 1.5-порождаемости     ║")
    print(f"║  n = {n}, d = {d}, trials = {num_trials}          ║")
    print(f"╚══════════════════════════════════════════╝")
    print()

    # Описание D
    print(f"  D = {_describe_generator(D_gen, n)}")
    print(f"  dim W_{n}^(≤{d}) = {dim_Wn_leq_d(n, d)}")
    print()

    successes = 0
    rng = np.random.default_rng(seed)

    for trial in range(num_trials):
        print(f"  ──── Trial {trial + 1}/{num_trials} ────")
        E_gen = make_random_generator(n, E_max_degree, rng=rng, sparsity=0.3)

        checker = NumericLieGenerationChecker(
            n, d,
            [D_gen, E_gen],
            D_max=max(d, E_max_degree + 1),
            verbose=True,
        )
        result = checker.run()

        if result.success:
            successes += 1
            print(f"  ✓ PASS ({result.elapsed_s:.3f}s, {result.iterations} iter)")
        else:
            print(f"  ✗ FAIL ({result.elapsed_s:.3f}s, {result.iterations} iter)")
            # Показать неполные степени
            for ds in result.degrees:
                if not ds.is_full:
                    print(f"    {ds}")
        print()

    print(f"  Итого: {successes}/{num_trials} успешных")
    if successes == num_trials:
        print(f"  ⟹ Гипотеза ПОДТВЕРЖДЕНА (численно, для данных E)")
    elif successes > 0:
        print(f"  ⟹ Частичное подтверждение")
    else:
        print(f"  ⟹ Все попытки ПРОВАЛЕНЫ — кандидат на контрпример?")
        print(f"     (Попробуйте увеличить d, E_max_degree или число попыток)")

    return successes


def _describe_generator(gen, n):
    """Краткое описание генератора."""
    vars = [f"z{i}" for i in range(n)]
    parts = []
    for i, mono_dict in gen.items():
        for alpha, coeff in mono_dict.items():
            mono_str = ""
            for j, a in enumerate(alpha):
                if a == 1:
                    mono_str += vars[j]
                elif a > 1:
                    mono_str += f"{vars[j]}^{a}"
            if not mono_str:
                mono_str = "1"
            c_str = "" if coeff == 1.0 else f"{coeff}·"
            parts.append(f"{c_str}{mono_str}·∂_{vars[i]}")
    return " + ".join(parts) if parts else "0"


# ──────────────────────────────────────────────────────────────────────────────
# Быстрый тест размерностей
# ──────────────────────────────────────────────────────────────────────────────

def check_dimensions():
    """Проверить таблицу размерностей из Numeric.md §2.3."""
    print("Проверка таблицы размерностей из Numeric.md:")
    expected = [
        (2, 5, 56), (2, 10, 156),
        (3, 5, 168), (3, 10, 1092),
        (5, 5, 2310), (5, 10, 48048),
    ]
    all_ok = True
    for n, d, exp in expected:
        actual = dim_Wn_leq_d(n, d)
        ok = "✓" if actual == exp else "✗"
        if actual != exp:
            all_ok = False
        print(f"  {ok} dim W_{n}^(≤{d}) = {actual} (ожидалось {exp})")
    print()
    return all_ok


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Числ. проверка порождения W_n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--mode", choices=["beldiev", "hypo", "dims"], default="beldiev",
                        help="Режим: beldiev, hypo (гипотеза), dims (проверка размерностей)")
    parser.add_argument("--n", type=int, default=2, help="Число переменных")
    parser.add_argument("--d", type=int, default=7, help="Проверять до степени d включительно")
    parser.add_argument("--Dmax", type=int, default=None, help="Максимальное усечение D_max")

    # Параметры для hypo-режима
    parser.add_argument("--Da", type=int, default=5, help="Экспонента a для D = z0^a·z1^b·∂_0")
    parser.add_argument("--Db", type=int, default=5, help="Экспонента b для D = z0^a·z1^b·∂_0")
    parser.add_argument("--trials", type=int, default=5, help="Число случайных E")
    parser.add_argument("--E-max-deg", type=int, default=None, help="Макс. степень E")
    parser.add_argument("--seed", type=int, default=42, help="Seed для генератора случайных чисел")

    args = parser.parse_args()

    if args.mode == "dims":
        ok = check_dimensions()
        sys.exit(0 if ok else 1)

    elif args.mode == "beldiev":
        result = run_beldiev(args.n, args.d, D_max=args.Dmax)
        sys.exit(0 if result.success else 1)

    elif args.mode == "hypo":
        D_gen = make_monomial_D(args.n, args.Da, args.Db)
        successes = run_hypothesis(
            n=args.n,
            d=args.d,
            D_gen=D_gen,
            num_trials=args.trials,
            E_max_degree=args.E_max_deg,
            seed=args.seed,
        )
        sys.exit(0 if successes > 0 else 1)


if __name__ == "__main__":
    main()
