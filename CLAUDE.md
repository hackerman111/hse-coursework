# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Research platform for investigating **Beldiev's Problem 3.2**: is the Lie algebra W_n of polynomial vector fields on K^n always 1.5-generated? That is, for any nonzero D ∈ W_n, does there exist E ∈ W_n such that Lie(D, E) = W_n?

Two known positive results are implemented:
- **Beldiev's theorem**: W_n is 2-generated (explicit formulas in `lib/generator/beldiev.py`)
- **Andrist's theorem**: W_n is 3-generated (explicit formulas in `lib/generator/andrist.py`)

The Lie bracket is: [f∂_{x_i}, g∂_{x_j}] = f·∂_{x_i}(g)·∂_{x_j} − g·∂_{x_j}(f)·∂_{x_i}

## Commands

```bash
# Run all tests
pytest tests/

# Run a single test file
pytest tests/test_numeric.py -v

# Run a single test
pytest tests/test_numeric.py::NumericSolverTests::test_beldiev_generators_fill_w2_up_to_degree_7 -s

# Symbolic tests require SageMath
pytest tests/test_gen/ tests/test_check/ -v

# CLI: verify Beldiev's generators numerically
python numeric_check.py --mode beldiev --n 2 --d 7

# CLI: test 1.5-generation hypothesis with fixed generators
python numeric_check.py --mode hypo --n 2 --d 8 --g1 "x^5*y^5*dx" --g2 "dx + x*dy"

# CLI: enumerate candidate generators
python numeric_enumerate.py --coeffs 1,-1 --max-terms 2 --limit 20
```

## Architecture

The library has two parallel computation layers:

### Symbolic layer (`lib/derivation/`, `lib/solver/`)
Uses SageMath. `LieDerivation` wraps a Sage derivation and adds `.bracket()`, `.leading_term`, `.degree()`. `LieBasisSolver` in `lib/solver/basis.py` runs a Gröbner-like algorithm: it maintains a reduced basis, generates pairwise brackets, and terminates when all ∂_{x_i} are found.

### Numeric layer (`lib/numeric/`)
Uses NumPy only. Works with finite-dimensional truncations W_n^(≤d). `NumericLieGenerationChecker` in `lib/numeric/solver.py` is the main entry point. Structure constants for the bracket are precomputed in `lib/numeric/bracket.py`. Rank is tracked via SVD in `lib/numeric/subspace.py`. Much faster than symbolic, but only proves generation up to degree d.

### Generator implementations (`lib/generator/`)
`get_Beldiev(ring)` and `get_Andristy(ring)` return explicit generator lists. `check(generators, max_iter=1000)` is a convenience wrapper around the symbolic solver.

### Key classes
- `LieDerivation` (`lib/derivation/base.py`) — symbolic derivation with Lie bracket
- `LieBasisSolver` (`lib/solver/basis.py`) — symbolic closure algorithm
- `NumericLieGenerationChecker` (`lib/numeric/solver.py`) — numeric truncated checker
- `OrthonormalBasis` (`lib/numeric/subspace.py`) — SVD-based rank tracker
- `LibraryLogger` (`lib/utils/logger.py`) — Rich console + JSONL logging

## Known Architectural Issues

Per `.agent/AGENT.md` (the authoritative development guide, written in Russian):
- **Incomplete termination criterion**: solver checks for ∂_{x_i} only, not full W_n generation
- **No degree-bound truncation** in symbolic solver (only `max_iter`)
- **Duplication**: generator formulas exist independently in both symbolic and numeric layers; bracket computation is not shared
- **No unified pipeline** between symbolic and numeric paths

The research direction is documented in `hypo.md` (counterexample candidates) and `Modal.md` (planned 5-layer refactoring). `.agent/AGENT.md` contains the full mathematical context and priority roadmap — consult it before making architectural changes.
