# Fix Phantom Brackets in Numeric Solver — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the numeric Lie-generation solver so it does not produce false positives for non-homogeneous generators by replacing the graded-pool bracket loop with a full-element-pair bracket algorithm.

**Architecture:** The current solver decomposes generators into homogeneous degree components and brackets all vectors in `B_p` against all vectors in `B_q` regardless of provenance, creating "phantom" brackets (e.g., `[pi_{-1}(G2), pi_0(G2)]` when `[G2, G2] = 0`). The fix maintains the graded `OrthonormalBasis` infrastructure for rank checking, but replaces the bracket loop: instead of bracketing basis pools, we track full elements (each a `Dict[int, np.ndarray]` mapping degree → coordinate vector) and bracket element pairs. Only newly-discovered elements are bracketed against existing ones (Buchberger-style new-vs-existing optimization). The helper `bracket_full_elements()` computes `[f, g]` by summing `bracket_vectors(f_p, p, g_q, q)` over all valid `(p, q)` degree pairs.

**Tech Stack:** Python 3.14, NumPy, pytest

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `lib/numeric/bracket.py` | Modify | Add `bracket_full_elements()` function |
| `lib/numeric/solver.py` | Modify | Replace `__init__` seeding + `run()` bracket loop with full-element algorithm |
| `lib/numeric/components.py` | Modify | Add docstring warning about decomposition semantics |
| `tests/test_numeric.py` | Modify | Add negative regression tests + bracket unit tests |

No new files created. No changes to `subspace.py`, `structure.py`, `indexing.py`, `criteria.py`, `workflows.py`, `specs.py`, `types.py`, `results.py`, or `__init__.py`. The public API is unchanged.

---

## Context for Implementers

### The Bug

`NumericLieGenerationChecker` decomposes each generator into homogeneous degree components (via `decompose_generator`) and seeds independent subspace bases `B[k]`. The bracket loop then computes `[B_p, B_q]` for ALL pairs of vectors across degree subspaces, regardless of which original element they came from. This produces "phantom" brackets.

**Example:** G2 = `d/dx + x*d/dy` is decomposed into B[-1] = {d/dx} and B[0] = {x*d/dy}. The solver computes `[d/dx, x*d/dy] = d/dy` and adds it to B[-1]. But `[G2, G2] = 0`, so `d/dy` is a phantom — it doesn't exist in `Lie(G1, G2)`.

### The Fix (Option C from review)

**Keep** the graded `OrthonormalBasis` structure for rank checking.
**Replace** the bracket source: instead of bracketing degree-pool bases, bracket tracked full elements.

A "full element" is `Dict[int, np.ndarray]` — the degree-decomposed form of one element of the generated subalgebra. The bracket of two full elements `f` and `g` produces a new full element where the degree-`r` component is `sum([f_p, g_q] for p+q=r)`.

**New-vs-existing optimization:** Maintain two lists — `processed` (already bracketed against everything) and `new` (not yet). Each iteration: bracket every `new` element against every `processed` element and against other `new` elements; move `new` into `processed`; any elements with new rank contributions become the next batch of `new`.

### Key Files to Read Before Implementing

- `lib/numeric/solver.py` — the main file being modified (lines 51-205)
- `lib/numeric/bracket.py:106-134` — `bracket_vectors()` which is the low-level building block
- `lib/numeric/components.py` — `decompose_generator()` still used for initial decomposition
- `lib/numeric/structure.py` — `StructureConstantCache` (unchanged, provides `get(p, q)`)
- `lib/numeric/subspace.py` — `OrthonormalBasis` (unchanged, used for rank tracking)
- `lib/numeric/indexing.py` — `dim_Wn_k()` (unchanged)
- `tests/test_numeric.py` — existing test suite (all tests must continue passing)

### Ground Truth for Existing Tests

All existing test cases with non-homogeneous generators have been verified symbolically:
- `y^3*dx + x^3*dy`, `dx + x*dy` → PASS (symbolic confirms generation)
- `x^3*dx + y^3*dy`, `dx + x*dy` → PASS (symbolic confirms generation)
- `x^4*dx + 2*x^2*y^2*dx + y^4*dx`, `dx + x*dy` → PASS (symbolic confirms generation)
- `x^2*dx`, `dx + x^3*dx` (n=1) → PASS (symbolic confirms generation)

The known false-positive case:
- `2*x*y^2*dx - x^3*y*dx`, `dx + x*dy` → **FAIL** (symbolic confirms non-generation, numeric incorrectly says PASS)

---

### Task 1: Add `bracket_full_elements()` to `bracket.py` + unit tests

**Files:**
- Modify: `lib/numeric/bracket.py` (add function after line 134)
- Modify: `tests/test_numeric.py` (add test class at end)

- [ ] **Step 1: Write the failing tests for `bracket_full_elements`**

Add to the end of `tests/test_numeric.py`, before the `if __name__` block:

```python
from lib.numeric.bracket import bracket_vectors, bracket_full_elements


class BracketVectorsTests(unittest.TestCase):
    """Unit tests for low-level bracket computation."""

    def test_partial_x_bracket_x_partial_y_gives_partial_y(self):
        """[d/dx, x*d/dy] = d/dy in W_2."""
        # d/dx in W_2^(-1): dim=2, index 0 = (0,0)*d/dz0
        u = np.array([1.0, 0.0])
        # x*d/dy in W_2^(0): dim=4, index 2 = z0 * d/dz1
        v = np.array([0.0, 0.0, 1.0, 0.0])
        result = bracket_vectors(u, -1, v, 0, 2)
        expected = np.array([0.0, 1.0])  # d/dy in W_2^(-1)
        np.testing.assert_allclose(result, expected)

    def test_commuting_diagonal_fields_give_zero(self):
        """[x*d/dx, y*d/dy] = 0 in W_2."""
        # x*d/dx in W_2^(0): dim=4, index 0 = z0*d/dz0
        u = np.array([1.0, 0.0, 0.0, 0.0])
        # y*d/dy in W_2^(0): dim=4, index 3 = z1*d/dz1
        v = np.array([0.0, 0.0, 0.0, 1.0])
        result = bracket_vectors(u, 0, v, 0, 2)
        np.testing.assert_allclose(result, np.zeros(6), atol=1e-14)

    def test_antisymmetry(self):
        """[u, v] = -[v, u] for bracket_vectors."""
        u = np.array([1.0, 0.0, 0.0, 0.0])  # z0*d/dz0 in W_2^(0)
        v = np.array([0.0, 1.0, 0.0, 0.0])  # z1*d/dz0 in W_2^(0)
        from lib.numeric.bracket import precompute_structure_constants
        sc_pq = precompute_structure_constants(2, 0, 0)
        fwd = bracket_vectors(u, 0, v, 0, 2, sc=sc_pq)
        rev = bracket_vectors(v, 0, u, 0, 2, sc=sc_pq)
        np.testing.assert_allclose(fwd, -rev, atol=1e-14)


class BracketFullElementsTests(unittest.TestCase):
    """Unit tests for full-element bracket."""

    def test_self_bracket_is_zero(self):
        """[G, G] must be zero for any element."""
        from lib.numeric.structure import StructureConstantCache
        sc = StructureConstantCache(n=2, D_max=3)
        # G = d/dx + x*d/dy (degrees -1 and 0)
        elem = {
            -1: np.array([1.0, 0.0]),               # d/dx
            0: np.array([0.0, 0.0, 1.0, 0.0]),      # x*d/dy
        }
        result = bracket_full_elements(elem, elem, 2, sc)
        for k, vec in result.items():
            np.testing.assert_allclose(
                vec, np.zeros_like(vec), atol=1e-14,
                err_msg=f"self-bracket nonzero at degree {k}",
            )

    def test_bracket_matches_bracket_vectors_for_homogeneous(self):
        """For single-degree elements, bracket_full_elements matches bracket_vectors."""
        from lib.numeric.structure import StructureConstantCache
        sc = StructureConstantCache(n=2, D_max=3)
        f = {-1: np.array([1.0, 0.0])}              # d/dx
        g = {0: np.array([0.0, 0.0, 1.0, 0.0])}     # x*d/dy
        result = bracket_full_elements(f, g, 2, sc)
        expected = bracket_vectors(
            np.array([1.0, 0.0]), -1,
            np.array([0.0, 0.0, 1.0, 0.0]), 0,
            2,
        )
        self.assertIn(-1, result)
        np.testing.assert_allclose(result[-1], expected)

    def test_cross_bracket_of_non_homogeneous_elements(self):
        """[G1, G2] where both span multiple degrees."""
        from lib.numeric.structure import StructureConstantCache
        sc = StructureConstantCache(n=2, D_max=5)
        # G1 = 2*x*y^2*d/dx - x^3*y*d/dx (degrees 2, 3)
        from lib.numeric.components import decompose_generator
        g1_spec = {0: {(1, 2): 2.0, (3, 1): -1.0}}
        g2_spec = {0: {(0, 0): 1.0}, 1: {(1, 0): 1.0}}
        elem1 = decompose_generator(g1_spec, 2, 5)
        elem2 = decompose_generator(g2_spec, 2, 5)
        result = bracket_full_elements(elem1, elem2, 2, sc)
        # Just verify it returns something and doesn't crash.
        # The result should have components at degrees 1,2,3 (from cross-terms).
        has_nonzero = any(np.linalg.norm(v) > 1e-14 for v in result.values())
        self.assertTrue(has_nonzero, "bracket of G1, G2 should be nonzero")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `PYTHONPATH=. python -m pytest tests/test_numeric.py::BracketVectorsTests tests/test_numeric.py::BracketFullElementsTests -v`
Expected: `BracketVectorsTests` pass (they use existing `bracket_vectors`), `BracketFullElementsTests` fail with `ImportError: cannot import name 'bracket_full_elements'`

- [ ] **Step 3: Implement `bracket_full_elements` in `bracket.py`**

Add after the `bracket_vectors` function (after line 134 of `lib/numeric/bracket.py`):

```python
def bracket_full_elements(
    f: Dict[int, np.ndarray],
    g: Dict[int, np.ndarray],
    n: int,
    structure_cache: "StructureConstantCache",
) -> Dict[int, np.ndarray]:
    """
    Compute the Lie bracket [f, g] of two full (non-homogeneous) elements.

    Each element is a dict {degree: coordinate_vector_in_Wn_k}.
    The bracket is computed as:
        [f, g]_r = sum_{p+q=r} [f_p, g_q]
    for each target degree r.

    This correctly handles non-homogeneous elements by summing ALL
    cross-degree contributions, unlike bracketing graded pools which
    can produce phantom brackets.

    Args:
        f: first element, {degree: np.ndarray}
        g: second element, {degree: np.ndarray}
        n: number of variables
        structure_cache: precomputed structure constants

    Returns:
        Dict[int, np.ndarray] — the bracket [f, g], keyed by degree.
        Only degrees with nonzero result are included.
    """
    result: Dict[int, np.ndarray] = {}

    for p, f_p in f.items():
        for q, g_q in g.items():
            r = p + q
            sc = structure_cache.get(p, q) if p <= q else structure_cache.get(q, p)
            if sc is None:
                continue

            if p <= q:
                contribution = bracket_vectors(f_p, p, g_q, q, n, sc=sc)
            else:
                # [f_p, g_q] = -[g_q, f_p] when using sc for (q, p)
                contribution = -bracket_vectors(g_q, q, f_p, p, n, sc=sc)

            N_r = dim_Wn_k(n, r)
            if r not in result:
                result[r] = np.zeros(N_r, dtype=np.float64)
            result[r] += contribution

    # Filter out zero components
    return {k: v for k, v in result.items() if np.linalg.norm(v) > 1e-15}
```

Also add the necessary import at the top of `bracket.py` (after line 7):

```python
from typing import Dict
```

Note: The `StructureConstantCache` type is used in the signature as a string annotation to avoid circular imports. It already works because of `from __future__ import annotations` at the top of the file.

- [ ] **Step 4: Run tests to verify they pass**

Run: `PYTHONPATH=. python -m pytest tests/test_numeric.py::BracketVectorsTests tests/test_numeric.py::BracketFullElementsTests -v`
Expected: all 6 tests PASS

- [ ] **Step 5: Commit**

```bash
git add lib/numeric/bracket.py tests/test_numeric.py
git commit -m "feat: bracket_full_elements для скобок полных (неоднородных) элементов"
```

---

### Task 2: Rewrite `NumericLieGenerationChecker.__init__` and `run()` with full-element algorithm

**Files:**
- Modify: `lib/numeric/solver.py:51-205` (rewrite `__init__` seeding and `run()` loop)

This is the core fix. The changes are:

1. In `__init__`: still decompose generators for initial rank seeding, but ALSO store them as full elements.
2. In `run()`: replace the graded-pool bracket loop with full-element-pair brackets using new-vs-existing optimization.

- [ ] **Step 1: Write the failing regression test**

Add to `tests/test_numeric.py` inside the `NumericSolverTests` class:

```python
    def test_non_generating_pair_returns_fail(self):
        """Generators that don't generate W_2 must return success=False.

        2*x*y^2*dx - x^3*y*dx and dx + x*dy do NOT generate W_2.
        Confirmed by symbolic solver (0/2 partial derivatives found).
        The old algorithm falsely reported PASS due to phantom brackets
        from graded decomposition of non-homogeneous generators.
        """
        first_gen = parse_generator_spec("2*x*y^2*dx - x^3*y*dx", 2)
        second_gen = parse_generator_spec("dx + x*dy", 2)
        checker = NumericLieGenerationChecker(
            2,
            7,
            [first_gen, second_gen],
            D_max=7,
            verbose=False,
        )

        result = checker.run()

        self.assertFalse(result.success)

    def test_trivially_incomplete_diagonal_pair_returns_fail(self):
        """x*dx and y*dy can only generate diagonal fields — never cross-terms."""
        first_gen = parse_generator_spec("x*dx", 2)
        second_gen = parse_generator_spec("y*dy", 2)
        checker = NumericLieGenerationChecker(
            2,
            2,
            [first_gen, second_gen],
            D_max=3,
            verbose=False,
        )

        result = checker.run()

        self.assertFalse(result.success)
```

- [ ] **Step 2: Run the new regression test to verify it fails (current bug)**

Run: `PYTHONPATH=. python -m pytest tests/test_numeric.py::NumericSolverTests::test_non_generating_pair_returns_fail -v`
Expected: FAIL with `AssertionError: True is not false` (the bug — solver incorrectly returns success=True)

- [ ] **Step 3: Rewrite `__init__` to track full elements**

Replace lines 79-90 of `lib/numeric/solver.py` (from `# Инициализация подпространств` through the generator seeding loop) with:

```python
        # Инициализация подпространств B[k] для k = -1 … D_max
        self.bases: Dict[int, OrthonormalBasis] = {}
        for k in range(-1, self.D_max + 1):
            N_k = dim_Wn_k(n, k)
            self.bases[k] = OrthonormalBasis(N_k)

        # Разложить генераторы в полные элементы и засеять подпространства
        self._elements: List[Dict[int, np.ndarray]] = []
        for gen in generators:
            elem = decompose_generator(gen, n, self.D_max)
            self._elements.append(elem)
            for k, vec in elem.items():
                if k in self.bases:
                    self.bases[k].try_add(vec)
```

- [ ] **Step 4: Rewrite `run()` with full-element bracket loop**

Replace the `run` method (lines 106-205 of `lib/numeric/solver.py`) with:

```python
    def run(self) -> NumericResult:
        """Запустить итеративное замыкание и вернуть результат."""
        from .bracket import bracket_full_elements

        t0 = time.perf_counter()
        iteration = 0
        next_detail_log_at = (
            t0 + self.log_every_s
            if self.verbose and self.logger.enabled and self.log_every_s is not None
            else None
        )

        if self.verbose and self.logger.enabled:
            self._print_status(iteration)

        # Алгоритм «новые-против-существующих» (Buchberger-style):
        # processed — элементы, уже прошедшие через все скобки;
        # pending   — элементы, которые ещё нужно прокоммутировать.
        processed: List[Dict[int, np.ndarray]] = []
        pending: List[Dict[int, np.ndarray]] = list(self._elements)

        while pending:
            iteration += 1
            next_pending: List[Dict[int, np.ndarray]] = []

            # Скобки: каждый pending × (все processed + остальные pending)
            for i, elem_new in enumerate(pending):
                # С каждым processed
                for elem_old in processed:
                    next_detail_log_at = self._maybe_print_detailed_status(
                        iteration, t0, next_detail_log_at,
                    )
                    bracket = bracket_full_elements(
                        elem_new, elem_old, self.n, self.structure_constants,
                    )
                    self._try_add_bracket(bracket, next_pending)

                # С другими pending (только i < j, чтобы не дублировать)
                for j in range(i + 1, len(pending)):
                    next_detail_log_at = self._maybe_print_detailed_status(
                        iteration, t0, next_detail_log_at,
                    )
                    bracket = bracket_full_elements(
                        elem_new, pending[j], self.n, self.structure_constants,
                    )
                    self._try_add_bracket(bracket, next_pending)

            processed.extend(pending)
            pending = next_pending

            if self.verbose and self.logger.enabled:
                self._print_status(iteration)
            next_detail_log_at = self._maybe_print_detailed_status(
                iteration, t0, next_detail_log_at,
            )

        elapsed = time.perf_counter() - t0

        # Собираем результат
        degrees: List[DegreeStatus] = []
        success = True
        for k in range(-1, self.d + 1):
            N_k = dim_Wn_k(self.n, k)
            rank = self.bases[k].rank
            ds = DegreeStatus(degree=k, rank=rank, target_dim=N_k)
            degrees.append(ds)
            if not ds.is_full:
                success = False

        return NumericResult(
            n=self.n,
            d=self.d,
            D_max=self.D_max,
            success=success,
            degrees=degrees,
            iterations=iteration,
            elapsed_s=elapsed,
        )

    def _try_add_bracket(
        self,
        bracket: Dict[int, np.ndarray],
        next_pending: List[Dict[int, np.ndarray]],
    ) -> None:
        """Попробовать добавить компоненты скобки в базисы.

        Если хотя бы одна компонента увеличила ранг соответствующего
        подпространства, весь элемент добавляется в next_pending
        для дальнейших итераций.
        """
        if not bracket:
            return
        added_any = False
        for k, vec in bracket.items():
            if k not in self.bases:
                continue
            norm = np.linalg.norm(vec)
            if norm < 1e-15:
                continue
            normalized = vec / norm
            if self.bases[k].try_add(normalized):
                added_any = True
        if added_any:
            next_pending.append(bracket)
```

- [ ] **Step 5: Run the full test suite**

Run: `PYTHONPATH=. python -m pytest tests/test_numeric.py -v`
Expected: ALL tests pass, including the two new regression tests.

If any of the existing "assertTrue(result.success)" tests fail, it means the new algorithm needs more iterations to converge. In that case, increase D_max for those specific tests. However, based on the symbolic verification, all should pass.

- [ ] **Step 6: Verify the false-positive case is fixed via CLI**

Run: `python numeric_check.py --mode hypo --n 2 --d 7 --g1 '2*x*y^2*dx - x^3*y*dx' --g2 'dx + x*dy'`
Expected: output should contain `FAIL` and exit code should be 1.

Also verify a known-good case still works:
Run: `python numeric_check.py --mode beldiev --n 2 --d 7`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add lib/numeric/solver.py tests/test_numeric.py
git commit -m "fix: исправить фантомные скобки в численном солвере

Заменён алгоритм скобок: вместо попарного перемножения степенных
пулов B_p × B_q (что создавало фантомные скобки для неоднородных
генераторов) — скобки вычисляются между полными элементами с
суммированием по всем парам степеней.

Применена оптимизация new-vs-existing (Buchberger-style):
каждый новый элемент коммутируется только с уже обработанными."
```

---

### Task 3: Add docstring warning to `decompose_generator` + clean up dead code

**Files:**
- Modify: `lib/numeric/components.py:13-37` (add docstring warning)
- Modify: `lib/numeric/solver.py` (remove now-unused `batch_bracket` import)

- [ ] **Step 1: Add docstring warning to `decompose_generator`**

Replace the docstring of `decompose_generator` in `lib/numeric/components.py` (lines 13-23):

```python
def decompose_generator(
    generator: Generator,
    n: int,
    D_max: int,
) -> dict[int, np.ndarray]:
    """
    Split a sparse generator into homogeneous numeric components.

    ``generator`` has the form ``{axis: {alpha: coeff}}`` and the result maps
    each degree ``k`` to a coordinate vector in ``W_n^(k)``.

    WARNING: The returned components are degree-projections of a single element,
    NOT independent elements of a Lie subalgebra.  Do NOT bracket components
    from different degrees as if they were separate subalgebra members — this
    produces phantom brackets.  Use ``bracket_full_elements`` from ``bracket.py``
    to correctly compute Lie brackets of non-homogeneous elements.
    """
```

- [ ] **Step 2: Remove unused `batch_bracket` import from `solver.py`**

In `lib/numeric/solver.py`, line 19, change:

```python
from .bracket import batch_bracket
```

to:

```python
# batch_bracket no longer used here; full-element brackets are computed
# via bracket_full_elements, imported lazily inside run().
```

If there are no other references to `batch_bracket` in `solver.py`, this is safe. (`batch_bracket` is still used by the public API functions `try_add_to_basis` and `check_degree_fullness` at the bottom of solver.py — check this. If they don't use it, remove the import entirely. If they do, keep it.)

Actually, looking at lines 280-333 of `solver.py`, the public primitives `try_add_to_basis` and `check_degree_fullness` do NOT use `batch_bracket` — they operate on `OrthonormalBasis` directly. So the import can be safely removed.

Remove line 19 (`from .bracket import batch_bracket`) entirely.

- [ ] **Step 3: Run the full test suite**

Run: `PYTHONPATH=. python -m pytest tests/test_numeric.py -v`
Expected: all tests pass (no behavioral change, just docs + cleanup)

- [ ] **Step 4: Commit**

```bash
git add lib/numeric/components.py lib/numeric/solver.py
git commit -m "docs: предупреждение о фантомных скобках в decompose_generator"
```

---

### Task 4: Run the full existing test suite as a final regression check

**Files:**
- No files modified — verification only.

- [ ] **Step 1: Run all numeric tests**

Run: `PYTHONPATH=. python -m pytest tests/test_numeric.py tests/test_numeric_enumerate.py -v`
Expected: all tests pass

- [ ] **Step 2: Run CLI smoke tests**

Run each command and verify output:

```bash
# Beldiev generators (homogeneous, known PASS)
python numeric_check.py --mode beldiev --n 2 --d 7
# Expected: PASS

# Known PASS with non-homogeneous generators
python numeric_check.py --mode hypo --n 2 --d 2 --g1 'x^3*dx + y^3*dy' --g2 'dx + x*dy' --Dmax 3
# Expected: PASS

# Known FAIL (the original false-positive case)
python numeric_check.py --mode hypo --n 2 --d 7 --g1 '2*x*y^2*dx - x^3*y*dx' --g2 'dx + x*dy'
# Expected: FAIL

# Criterion mode with known PASS
python numeric_check.py --mode criterion --n 2 --g1 'x^3*dx + y^3*dy' --g2 'dx + x*dy' --Dmax 3
# Expected: PASS

# Trivially incomplete diagonal
python numeric_check.py --mode hypo --n 2 --d 2 --g1 'x*dx' --g2 'y*dy' --Dmax 3
# Expected: FAIL
```

- [ ] **Step 3: Run symbolic tests (if SageMath available)**

Run: `PYTHONPATH=. python -m pytest tests/test_solver_stop_condition.py tests/test_symbolic_check.py -v`
Expected: all pass (these test the symbolic layer, which is unaffected)
