#!/usr/bin/env python3
"""
numeric_check.py — быстрая проверка гипотез о порождении W_n.

Режимы:
  1. beldiev  — проверка генераторов Бельдиева (известный результат)
  2. hypo     — проверка гипотезы 1.5-порождаемости:
                задаём G1 и либо фиксируем G2, либо выбираем случайный G2
                и проверяем Lie(G1, G2) ⊇ W_n^{(≤d)}

Примеры:
  python numeric_check.py --mode beldiev --n 2 --d 7
  python numeric_check.py --mode hypo --n 2 --d 10 --g1 "x^5*y^5*dx" --trials 5 --E-max-deg 12
  python numeric_check.py --mode hypo --n 2 --d 8 --g1 "x^5*y^5*dx" --g2 "dx + x*dy"
  python numeric_check.py --mode hypo --n 1 --d 5 --g1 "x^2*dx" --g2 "dx + x^3*dx"
"""

from __future__ import annotations

import argparse
import os
import re
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from lib.numeric.solver import (
    NumericLieGenerationChecker,
    make_beldiev_generators,
    make_random_generator,
)
from lib.numeric.indexing import dim_Wn_leq_d


VARIABLE_ALIASES = ("x", "y", "z", "u", "v", "w")


class RichHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    """Help formatter with preserved examples and visible defaults."""


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


def make_general_monomial_generator(n: int, exponents: list[int], target_var: int = 0):
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


def parse_generator_spec(spec: str, n: int):
    """
    Разобрать строку вида:
      x^5*y^5*dx
      dx + x*dy
      2*z0^3*d1 - 0.5*z1*dz0
    """
    text = spec.replace(" ", "")
    if not text:
        raise ValueError("пустая спецификация генератора")

    normalized = text if text[0] in "+-" else f"+{text}"
    terms = re.findall(r"[+-][^+-]+", normalized)
    if not terms:
        raise ValueError(f"не удалось разобрать генератор: {spec!r}")

    generator = {}
    for raw_term in terms:
        sign = -1.0 if raw_term[0] == "-" else 1.0
        body = raw_term[1:]
        factors = [factor for factor in body.split("*") if factor]
        if not factors:
            raise ValueError(f"пустой член в генераторе: {spec!r}")

        coeff = sign
        alpha = [0] * n
        target_var = None

        for factor in factors:
            if re.fullmatch(r"\d+(?:\.\d+)?", factor):
                coeff *= float(factor)
                continue

            deriv_index = _parse_derivative_factor(factor, n)
            if deriv_index is not None:
                if target_var is not None:
                    raise ValueError(
                        f"в одном члене больше одной производной: {raw_term!r}"
                    )
                target_var = deriv_index
                continue

            var_index, exponent = _parse_variable_factor(factor, n)
            if var_index is None:
                raise ValueError(f"неизвестный фактор {factor!r} в {raw_term!r}")
            alpha[var_index] += exponent

        if target_var is None:
            raise ValueError(f"в члене {raw_term!r} не указана производная")

        generator.setdefault(target_var, {})
        alpha_tuple = tuple(alpha)
        generator[target_var][alpha_tuple] = (
            generator[target_var].get(alpha_tuple, 0.0) + coeff
        )
        if abs(generator[target_var][alpha_tuple]) <= 1e-15:
            del generator[target_var][alpha_tuple]
        if not generator[target_var]:
            del generator[target_var]

    if not generator:
        raise ValueError(f"генератор {spec!r} занулился после сокращения коэффициентов")
    return generator


def _parse_derivative_factor(token: str, n: int):
    """Распознать d0, dz0, dx, dy, ..."""
    if token.startswith("d") and token[1:].isdigit():
        index = int(token[1:])
        return _validate_axis(index, n, token)

    if token.startswith("dz") and token[2:].isdigit():
        index = int(token[2:])
        return _validate_axis(index, n, token)

    alias_index = _alias_index(token[1:]) if token.startswith("d") else None
    if alias_index is not None:
        return _validate_axis(alias_index, n, token)

    return None


def _parse_variable_factor(token: str, n: int):
    """Распознать z0^3, x^5, y, ..."""
    exponent = 1
    base = token
    if "^" in token:
        base, exponent_str = token.split("^", 1)
        if not exponent_str.isdigit():
            raise ValueError(f"степень должна быть целой неотрицательной: {token!r}")
        exponent = int(exponent_str)

    if base.startswith("z") and base[1:].isdigit():
        index = int(base[1:])
        return _validate_axis(index, n, token), exponent

    alias_index = _alias_index(base)
    if alias_index is not None:
        return _validate_axis(alias_index, n, token), exponent

    if base == "1":
        return 0, 0
    return None, exponent


def _alias_index(name: str):
    try:
        return VARIABLE_ALIASES.index(name)
    except ValueError:
        return None


def _validate_axis(index: int, n: int, token: str) -> int:
    if index < 0 or index >= n:
        raise ValueError(f"индекс {index} вне диапазона 0..{n - 1} в {token!r}")
    return index


def _resolve_hypothesis_generator(n: int, spec: str | None, a: int, b: int):
    """Приоритет: явная спецификация, иначе legacy-пара --Da/--Db."""
    if spec:
        return parse_generator_spec(spec, n)
    return make_monomial_D(n, a, b)


def _generator_term_count(gen) -> int:
    return sum(len(mono_dict) for mono_dict in gen.values())


def _generator_max_degree(gen) -> int:
    max_degree = -1
    for mono_dict in gen.values():
        for alpha in mono_dict:
            max_degree = max(max_degree, sum(alpha) - 1)
    return max_degree


def _describe_random_generator(gen, n: int) -> str:
    return (
        f"random terms={_generator_term_count(gen)}, "
        f"max_degree={_generator_max_degree(gen)}, "
        f"targets={sorted(gen.keys())}, "
        f"example={_describe_generator(gen, n, max_terms=4)}"
    )


# ──────────────────────────────────────────────────────────────────────────────
# Режим 1: проверка генераторов Бельдиева
# ──────────────────────────────────────────────────────────────────────────────

def run_beldiev(
    n: int,
    d: int,
    D_max: int | None = None,
    log_every_s: float | None = None,
):
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

    checker = NumericLieGenerationChecker(
        n,
        d,
        generators,
        D_max=D_max,
        verbose=True,
        log_every_s=log_every_s,
    )
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
    first_gen,
    second_gen=None,
    num_trials: int = 5,
    E_max_degree: int | None = None,
    seed: int = 42,
    log_every_s: float | None = None,
):
    """
    Проверить, порождает ли Lie(G1, G2) = W_n^{(≤d)}.

    Если second_gen не задан, G2 берётся случайным на каждом trial.
    """
    if E_max_degree is None:
        E_max_degree = d + 2
    trial_count = 1 if second_gen is not None else num_trials

    print(f"╔══════════════════════════════════════════╗")
    print(f"║  Проверка гипотезы 1.5-порождаемости     ║")
    print(f"║  n = {n}, d = {d}, trials = {trial_count}          ║")
    print(f"╚══════════════════════════════════════════╝")
    print()

    print(f"  G1 = {_describe_generator(first_gen, n)}")
    print(
        f"  G1 stats: terms={_generator_term_count(first_gen)}, "
        f"max_degree={_generator_max_degree(first_gen)}"
    )
    if second_gen is None:
        print(
            f"  G2 mode = random (trials={num_trials}, max_degree={E_max_degree}, seed={seed})"
        )
    else:
        print(f"  G2 mode = fixed")
        print(f"  G2 = {_describe_generator(second_gen, n)}")
        print(
            f"  G2 stats: terms={_generator_term_count(second_gen)}, "
            f"max_degree={_generator_max_degree(second_gen)}"
        )
    print(f"  dim W_{n}^(≤{d}) = {dim_Wn_leq_d(n, d)}")
    print()

    successes = 0
    rng = np.random.default_rng(seed)
    if second_gen is not None and num_trials != 1:
        print(f"  note: fixed G2 задан, поэтому --trials={num_trials} игнорируется")
        print()

    for trial in range(trial_count):
        print(f"  ──── Trial {trial + 1}/{trial_count} ────")
        if second_gen is None:
            current_second = make_random_generator(
                n,
                E_max_degree,
                rng=rng,
                sparsity=0.3,
            )
            print(f"  G2 random: {_describe_random_generator(current_second, n)}")
        else:
            current_second = second_gen

        current_max_degree = max(
            _generator_max_degree(first_gen),
            _generator_max_degree(current_second),
            E_max_degree if second_gen is None else -1,
        )

        checker = NumericLieGenerationChecker(
            n,
            d,
            [first_gen, current_second],
            D_max=max(d, current_max_degree + 1),
            verbose=True,
            log_every_s=log_every_s,
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

    print(f"  Итого: {successes}/{trial_count} успешных")
    if successes == trial_count:
        print(f"  ⟹ Гипотеза ПОДТВЕРЖДЕНА (численно, для данных E)")
    elif successes > 0:
        print(f"  ⟹ Частичное подтверждение")
    else:
        print(f"  ⟹ Все попытки ПРОВАЛЕНЫ — кандидат на контрпример?")
        print(f"     (Попробуйте увеличить d, E_max_degree или число попыток)")

    return successes


def _describe_generator(gen, n, max_terms: int | None = None):
    """Краткое описание генератора."""
    vars = [f"z{i}" for i in range(n)]
    parts = []
    for i, mono_dict in sorted(gen.items()):
        for alpha, coeff in sorted(mono_dict.items()):
            factors = []
            for j, a in enumerate(alpha):
                if a == 1:
                    factors.append(vars[j])
                elif a > 1:
                    factors.append(f"{vars[j]}^{a}")
            mono_str = "*".join(factors) if factors else "1"
            c_str = "" if coeff == 1.0 else f"{coeff:g}·"
            parts.append(f"{c_str}{mono_str}·∂_{vars[i]}")
    if max_terms is not None and len(parts) > max_terms:
        return " + ".join(parts[:max_terms]) + f" + ... ({len(parts)} terms)"
    return " + ".join(parts) if parts else "0"


# ──────────────────────────────────────────────────────────────────────────────
# Быстрый тест размерностей
# ──────────────────────────────────────────────────────────────────────────────

def check_dimensions():
    """Проверить таблицу размерностей из Numeric.md §2.3."""
    print("Проверка таблицы размерностей из Numeric.md:")
    expected = [
        (2, 5, 56), (2, 10, 156),
        (3, 5, 252), (3, 10, 1092),
        (5, 5, 2310), (5, 10, 21840),
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
    description = (
        "Численная проверка порождения W_n.\n\n"
        "Режимы:\n"
        "  beldiev : эталонная проверка известных генераторов Бельдиева.\n"
        "  hypo    : проверка пары генераторов G1, G2. Если G2 не задан,\n"
        "            он генерируется случайно на каждом trial.\n"
        "  dims    : сверка таблицы размерностей из Numeric.md.\n\n"
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
        choices=["beldiev", "hypo", "dims"],
        default="beldiev",
        help="Режим запуска",
    )
    common.add_argument("--n", type=int, default=2, help="Число переменных")
    common.add_argument("--d", type=int, default=7, help="Проверять до степени d включительно")
    common.add_argument("--Dmax", type=int, default=None, help="Максимальное усечение D_max")
    common.add_argument(
        "--log-every",
        type=float,
        default=5.0,
        help="Раз в сколько секунд печатать полную сводку по степеням; 0 отключает",
    )

    hypo_group = parser.add_argument_group("Параметры режима hypo")
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
        help="Второй генератор в формате generator-spec; если не задан, используется случайный G2",
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
        ok = check_dimensions()
        sys.exit(0 if ok else 1)

    elif args.mode == "beldiev":
        result = run_beldiev(
            args.n,
            args.d,
            D_max=args.Dmax,
            log_every_s=args.log_every,
        )
        sys.exit(0 if result.success else 1)

    elif args.mode == "hypo":
        try:
            first_gen = _resolve_hypothesis_generator(args.n, args.g1, args.Da, args.Db)
            second_gen = parse_generator_spec(args.g2, args.n) if args.g2 else None
        except ValueError as exc:
            parser.error(str(exc))

        successes = run_hypothesis(
            n=args.n,
            d=args.d,
            first_gen=first_gen,
            second_gen=second_gen,
            num_trials=args.trials,
            E_max_degree=args.E_max_deg,
            seed=args.seed,
            log_every_s=args.log_every,
        )
        sys.exit(0 if successes > 0 else 1)


if __name__ == "__main__":
    main()
