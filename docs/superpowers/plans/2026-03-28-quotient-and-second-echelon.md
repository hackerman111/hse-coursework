# Quotient And Second Echelon Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote quotient derivations to production code, add the second-echelon spec/backend API, and restore the intended layer boundaries around `DerivationSpec`.

**Architecture:** Introduce a small symbolic `PolynomialOps` adapter for coefficient inspection, keep `DerivationSpec` as a DTO with structural helpers only, and move backend-shaped helpers into backend modules while preserving the existing public surface where compatibility matters.

**Tech Stack:** Python, pytest, NumPy, optional SageMath integration

---

### Task 1: Write red tests for quotient adapter behavior

**Files:**
- Modify: `tests/test_backends.py`
- Modify: `tests/test_core_spec.py`

- [ ] **Step 1: Write the failing tests**

Add tests that assert:

```python
from lib.derivation import LieDerivationFactory, QuotientLieDerivation


def test_factory_uses_quotient_wrapper_for_quotient_codomain():
    fake = FakeSageDerivation(...)
    wrapped = LieDerivationFactory.create(fake)
    assert isinstance(wrapped, QuotientLieDerivation)
```

and that quotient leading-term/degree inspection uses lifted representatives.

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_backends.py -k quotient -q`
Expected: FAIL because the adapter-based behavior and test doubles are not yet implemented.

- [ ] **Step 3: Write minimal implementation**

Add symbolic polynomial-ops adapters and route `LieDerivation` / factory through them.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_backends.py -k quotient -q`
Expected: PASS

### Task 2: Write red tests for second-echelon spec helpers

**Files:**
- Modify: `tests/test_core_spec.py`

- [ ] **Step 1: Write the failing tests**

Add tests for:

```python
def test_zero_spec_returns_empty_terms():
    assert zero_spec(3).terms == {}


def test_homogeneous_components_groups_by_degree():
    ...


def test_leading_term_spec_returns_none_for_zero():
    ...
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_core_spec.py -k "zero_spec or homogeneous or leading_term_spec or normalize_spec" -q`
Expected: FAIL because the new helpers are absent.

- [ ] **Step 3: Write minimal implementation**

Implement only the helpers required for the failing assertions.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_core_spec.py -k "zero_spec or homogeneous or leading_term_spec or normalize_spec" -q`
Expected: PASS

### Task 3: Write red tests for backend-layer ownership and numeric reverse conversion

**Files:**
- Modify: `tests/test_backends.py`
- Modify: `tests/test_core_spec.py`

- [ ] **Step 1: Write the failing tests**

Add tests that assert:

```python
from lib.backends.numeric import from_numeric_components, to_generator
from lib.backends.spec_ops import lie_bracket, ad


def test_numeric_roundtrip():
    recovered = from_numeric_components(2, to_numeric_components(spec, D_max=5))
    assert recovered == spec
```

and compatibility imports remain usable.

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_backends.py tests/test_core_spec.py -k "roundtrip or to_generator or lie_bracket or ad" -q`
Expected: FAIL because the functions still live in the old module layout or reverse conversion is missing.

- [ ] **Step 3: Write minimal implementation**

Move the helper implementations to backend-oriented modules, re-export compatibility shims, and add reverse numeric decoding.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_backends.py tests/test_core_spec.py -k "roundtrip or to_generator or lie_bracket or ad" -q`
Expected: PASS

### Task 4: Full verification

**Files:**
- Modify as needed based on prior tasks

- [ ] **Step 1: Run focused verification**

Run: `pytest tests/test_core_spec.py tests/test_backends.py tests/test_solver_api.py -q`
Expected: PASS

- [ ] **Step 2: Run any broader targeted verification if new failures appear**

Run: `pytest -q`
Expected: PASS or a clearly reported existing unrelated failure set.
