# Quotient And Second Echelon Design

## Scope

This change promotes quotient support from a manual-test convenience into production code and completes the agreed second-echelon API around `DerivationSpec`.

It also restores the intended layer boundary by moving three non-DTO helpers out of `lib/core/spec.py`:

- `lie_bracket`
- `ad`
- `to_generator`

## Design

### 1. Symbolic quotient support

`LieDerivation` stops knowing how to inspect polynomial coefficients directly.
Instead it receives a small `PolynomialOps` adapter responsible for:

- converting an algebra element into a polynomial representative,
- extracting the leading monomial/coefficient,
- extracting degree.

Two implementations are added:

- `BasePolynomialOps` for ordinary polynomial rings,
- `QuotientPolynomialOps` for quotient rings, using `.lift()` when available.

`LieDerivationFactory` selects the correct adapter by inspecting `codomain()`.
`QuotientLieDerivation` remains public for compatibility, but becomes a thin specialization over the same adapter-based logic instead of owning bespoke behavior.

### 2. DTO layer cleanup

`DerivationSpec` remains a DTO plus small structural helpers.
The following helpers move out of `lib/core/spec.py`:

- `lie_bracket`, `ad` -> `lib/backends/spec_ops.py`
- `to_generator` -> `lib/backends/numeric.py`

Public imports stay available from `lib` and `lib.core` where compatibility is practical, but the implementation ownership moves to backend-oriented modules.

### 3. Second-echelon spec API

The following functions are added to `lib/core/spec.py`:

- `zero_spec`
- `scale_spec`
- `normalize_spec`
- `homogeneous_component`
- `homogeneous_components`
- `leading_term_spec`

`leading_term_spec` is a pure spec-level inspection helper:

- input: `DerivationSpec`
- output: `LeadingTerm | None`
- zero derivation -> `None`
- default order: `lex`
- tie-breaker: larger `axis` wins after equal monomial

`LeadingTerm` stores `axis`, `alpha`, `coeff`, and `degree`.

### 4. Numeric reverse adapter

`lib/backends/numeric.py` gains `from_numeric_components(n, components)`.
It reconstructs a `DerivationSpec` by decoding each coordinate vector back into basis terms using the existing numeric indexing order.

## Testing

The implementation is driven by tests for:

- quotient leading term and degree via fake ring/algebra objects,
- factory selection of quotient-aware derivations,
- second-echelon spec helpers,
- numeric round-trip `spec -> to_numeric_components -> from_numeric_components -> spec`,
- compatibility imports after moving the three helper functions.
