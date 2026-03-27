---
name: numeric-lie-generation
description: Use this skill when you need to check Lie-generation numerically in this repository, parse generator specs like `x^a*y^b*dx`, apply the project criterion "all derivations up to degree 2 imply the whole Lie algebra", run `numeric_check.py` or `numeric_enumerate.py`, or extend the public numeric API under `lib/numeric`.
---

# Numeric Lie Generation

This skill is for the numeric layer in this repository. Prefer the public API in `lib/numeric` and the thin CLIs `numeric_check.py` and `numeric_enumerate.py`. Do not reimplement generator parsing, Beldiev recipes, random generator creation, or the degree-2 criterion inside ad hoc scripts.

## Use This Skill For

- Checking whether a fixed set of generators fills `W_n^(<= d)` numerically.
- Applying the project criterion: if all derivations up to degree `2` are generated, treat the set as generating the whole Lie algebra.
- Parsing or formatting generator specs like `dx + x*dy` or `x^5*y^5*dx`.
- Enumerating finite search spaces of candidate generators in dimension `2`.
- Refactoring or extending the numeric API without pushing more logic into CLI files.

## Preferred Entry Points

- CLI quick check: `python numeric_check.py --mode criterion --g1 "..." --g2 "..."`
- CLI truncated check: `python numeric_check.py --mode hypo --n 2 --d 2 --g1 "..." --g2 "..."`
- CLI enumeration: `python numeric_enumerate.py ...`
- Python API:
  - `lib.numeric.check_numeric_generation`
  - `lib.numeric.check_full_lie_generation`
  - `lib.numeric.parse_generator_spec`
  - `lib.numeric.make_beldiev_generators`
  - `lib.numeric.StructureConstantCache`

## Workflow

1. For a yes/no answer about the whole Lie algebra, use the degree-2 criterion first.
2. For raw layer-by-layer ranks, use `check_numeric_generation(...)` or `--mode hypo`.
3. If you need reusable structure constants across many runs, create one `StructureConstantCache` and pass it into repeated checks.
4. If you change the numeric API, keep `numeric_check.py` and `numeric_enumerate.py` as thin wrappers over `lib.numeric`.

## References

- Public API and examples: `references/api.md`
