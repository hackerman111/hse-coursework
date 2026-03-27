# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Research platform for investigating **Beldiev's Problem 3.2**: is the Lie algebra W_n of polynomial vector fields on K^n always 1.5-generated? That is, for any nonzero D in W_n, does there exist E in W_n such that Lie(D, E) = W_n?

Two known positive results are implemented:
- **Beldiev's theorem**: W_n is 2-generated (explicit formulas in `lib/generator/beldiev.py`)
- **Andrist's theorem**: W_n is 3-generated (explicit formulas in `lib/generator/andrist.py`)

The Lie bracket is: [f*d/dx_i, g*d/dx_j] = f*(dg/dx_i)*d/dx_j - g*(df/dx_j)*d/dx_i

**Language note**: code comments, docstrings, and `.agent/AGENT.md` are in Russian. Commit messages are also in Russian.

## Commands

```bash
# Run all numeric tests (no SageMath needed)
pytest tests/test_numeric.py -v
pytest tests/test_numeric_enumerate.py -v

# Run a single test
pytest tests/test_numeric.py::NumericSolverTests::test_beldiev_generators_fill_w2_up_to_degree_7 -s

# Symbolic tests (require SageMath)
pytest tests/test_gen/ tests/test_check/ -v

# CLI: verify Beldiev's generators numerically
python numeric_check.py --mode beldiev --n 2 --d 7

# CLI: test 1.5-generation hypothesis with fixed generators
python numeric_check.py --mode hypo --n 2 --d 8 --g1 "x^5*y^5*dx" --g2 "dx + x*dy"

# CLI: check full Lie generation criterion (generation up to degree 2 implies full generation)
python numeric_check.py --mode criterion --n 2 --g1 "x^5*y^5*dx" --g2 "dx + x*dy"

# CLI: enumerate candidate generators over a finite search space
python numeric_enumerate.py --coeffs 1,-1 --max-terms 2 --limit 20
python numeric_enumerate.py --count-only
```

### Generator spec format

The CLI and `parse_generator_spec()` accept human-readable strings like:
- `"x^5*y^5*dx"` = x^5 y^5 d/dx
- `"dx + x*dy"` = d/dx + x*d/dy
- `"3*x^2*y*dx - y^3*dy"` = 3x^2 y d/dx - y^3 d/dy

Variables: `x,y,z,u,v,w` (positional, mapped to indices 0..5). Partial derivatives: `dx,dy,dz,...`

## Architecture

Two parallel computation layers that do not share code:

### Symbolic layer (`lib/derivation/`, `lib/solver/`)
Requires SageMath. `LieDerivation` wraps a Sage derivation and adds `.bracket()`, `.leading_term`, `.degree()`. `LieBasisSolver` in `lib/solver/basis.py` runs a Grobner-like algorithm: maintains a reduced basis, generates pairwise brackets, terminates when all d/dx_i are found.

### Numeric layer (`lib/numeric/`)
NumPy only. Works with finite-dimensional truncations W_n^(<=d). Much faster than symbolic, but only proves generation up to degree d.

Key modules in the numeric layer:
- `solver.py` — `NumericLieGenerationChecker`, the main entry point. Runs bracket closure with SVD-based rank tracking
- `bracket.py` — precomputed structure constants for the Lie bracket
- `subspace.py` — `OrthonormalBasis`, SVD-based rank tracker
- `indexing.py` — multiindex enumeration and dimension formulas (`dim_Wn_k`)
- `specs.py` — generator builders and parser (`parse_generator_spec`, `make_monomial_generator`, `make_random_E`)
- `recipes.py` — `make_beldiev_generators` (numeric duplicate of `lib/generator/beldiev.py`)
- `criteria.py` — project-level generation criterion: if all derivations up to degree `FULL_LIE_GENERATION_DEGREE=2` are generated, treat as full Lie generation
- `workflows.py` — high-level orchestration (`run_beldiev`, `run_hypothesis`, `check_full_lie_generation`)
- `enumeration.py` — finite search over candidate generators
- `results.py` — `DegreeStatus`, `NumericResult` data classes
- `types.py` — `Generator` type alias: `Dict[int, Dict[Tuple[int,...], float]]` mapping variable index -> {multiindex: coefficient}

### Generator implementations (`lib/generator/`)
`get_Beldiev(ring)` and `get_Andristy(ring)` return explicit generator lists. `check(generators, max_iter=1000)` is a convenience wrapper around the symbolic solver.

### Logging (`lib/utils/logger.py`)
`LibraryLogger` — Rich console + JSONL logging with periodic status snapshots.

## Key data representation

A numeric generator is `Dict[int, Dict[Tuple[int,...], float]]`:
```python
# x^2*y*dx + 3*dy  in W_2
{0: {(2,1): 1.0}, 1: {(0,0): 3.0}}
# component 0 = d/dx direction, (2,1) = x^2*y^1, coefficient 1.0
```

## Known Architectural Issues

Per `.agent/AGENT.md` (the authoritative development guide, written in Russian):
- **Incomplete termination criterion**: solver checks for d/dx_i only, not full W_n generation
- **No degree-bound truncation** in symbolic solver (only `max_iter`)
- **Duplication**: generator formulas exist independently in both symbolic and numeric layers; bracket computation is not shared
- **No unified pipeline** between symbolic and numeric paths

The research direction is documented in `hypo.md` (counterexample candidates) and `Modal.md` (planned 5-layer refactoring). `.agent/AGENT.md` contains the full mathematical context and priority roadmap — consult it before making architectural changes.
