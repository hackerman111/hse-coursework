#!/usr/bin/env python3
"""
numeric_enumerate.py -- finite search over G1 = P(x,y)dx + Q(x,y)dy.

Important:
  Enumerating "all polynomials" is infinite unless coefficients are restricted.
  This script therefore enumerates a finite search space defined by:
    1. bounds on monomial exponents in x and y,
    2. a finite set of allowed nonzero coefficients,
    3. a bound on the number of nonzero monomials in G1.

Examples:
  python numeric_enumerate.py --count-only
  python numeric_enumerate.py --coeffs 1,-1 --max-terms 2 --limit 20
  python numeric_enumerate.py --coeffs 1 --max-terms 1 --match fail
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib.numeric.enumeration import (
    build_basis_terms,
    build_generator,
    candidate_to_spec,
    estimate_candidate_count,
    format_coefficient,
    iter_candidates,
    parse_coefficients,
)
from lib.numeric.specs import generator_max_degree, parse_generator_spec
from lib.numeric.structure import StructureConstantCache
from lib.numeric.workflows import check_numeric_generation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-x-deg", type=int, default=7)
    parser.add_argument("--max-y-deg", type=int, default=7)
    parser.add_argument(
        "--include-constants",
        action="store_true",
        help="include constant dx and dy terms in the search basis",
    )
    parser.add_argument(
        "--coeffs",
        type=str,
        default="1,-1",
        help="comma-separated list of allowed nonzero coefficients",
    )
    parser.add_argument(
        "--min-terms",
        type=int,
        default=1,
        help="minimum number of nonzero monomials in G1",
    )
    parser.add_argument(
        "--max-terms",
        type=int,
        default=3,
        help="maximum number of nonzero monomials in G1",
    )
    parser.add_argument(
        "--g2",
        type=str,
        default="dx + x*dy",
        help="fixed second generator",
    )
    parser.add_argument(
        "--d",
        type=int,
        default=2,
        help="check Lie generation up to degree d",
    )
    parser.add_argument(
        "--Dmax",
        type=int,
        default=None,
        help="optional truncation degree override",
    )
    parser.add_argument(
        "--match",
        choices=["pass", "fail", "all"],
        default="fail",
        help="which candidates to print",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="maximum number of matching candidates to print; 0 means no limit",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help="print progress every N processed candidates; 0 disables progress",
    )
    parser.add_argument(
        "--count-only",
        action="store_true",
        help="only print basis size and search-space size, do not run the checker",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.max_x_deg < 0 or args.max_y_deg < 0:
        raise ValueError("max degrees must be nonnegative")
    if args.min_terms < 1:
        raise ValueError("--min-terms must be at least 1")
    if args.max_terms < args.min_terms:
        raise ValueError("--max-terms must be >= --min-terms")

    coeffs = parse_coefficients(args.coeffs)
    basis_terms = build_basis_terms(
        max_x_deg=args.max_x_deg,
        max_y_deg=args.max_y_deg,
        include_constants=args.include_constants,
    )

    if args.max_terms > len(basis_terms):
        raise ValueError("--max-terms exceeds the size of the search basis")

    total_candidates = estimate_candidate_count(
        basis_size=len(basis_terms),
        coeff_count=len(coeffs),
        min_terms=args.min_terms,
        max_terms=args.max_terms,
    )

    print("=== Search Space ===")
    print(f"basis terms      = {len(basis_terms)}")
    print(f"coefficients     = {', '.join(format_coefficient(c) for c in coeffs)}")
    print(f"nonzero terms    = {args.min_terms}..{args.max_terms}")
    print(f"candidate count  = {total_candidates}")

    if args.count_only:
        return 0

    g2 = parse_generator_spec(args.g2, 2)
    dmax = args.Dmax
    if dmax is None:
        dmax = max(args.d, args.max_x_deg + args.max_y_deg - 1, generator_max_degree(g2))

    print(f"G2              = {args.g2}")
    print(f"d, Dmax         = {args.d}, {dmax}")
    print()

    structure_constants = StructureConstantCache(n=2, D_max=dmax)

    processed = 0
    passed = 0
    failed = 0
    printed = 0

    for candidate in iter_candidates(
        basis_terms=basis_terms,
        coeffs=coeffs,
        min_terms=args.min_terms,
        max_terms=args.max_terms,
    ):
        processed += 1
        g1 = build_generator(candidate)
        result = check_numeric_generation(
            n=2,
            d=args.d,
            generators=[g1, g2],
            D_max=dmax,
            verbose=False,
            structure_constants=structure_constants,
        )

        status = "pass" if result.success else "fail"
        if result.success:
            passed += 1
        else:
            failed += 1

        should_print = args.match == "all" or args.match == status
        if should_print:
            spec = candidate_to_spec(candidate)
            degree_summary = ", ".join(
                f"{degree.degree}:{degree.rank}/{degree.target_dim}"
                for degree in result.degrees
            )
            print(f"[{status.upper()}] {spec}")
            print(f"  degrees = {degree_summary}")
            printed += 1
            if args.limit and printed >= args.limit:
                break

        if args.progress_every and processed % args.progress_every == 0:
            print(
                f"... processed={processed}, pass={passed}, fail={failed}, printed={printed}"
            )

    print()
    print("=== Summary ===")
    print(f"processed = {processed}")
    print(f"pass      = {passed}")
    print(f"fail      = {failed}")
    print(f"printed   = {printed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
