# Numeric API

Use these public modules instead of reaching into private helpers or copying logic out of CLI scripts.

## Main functions

- `lib.numeric.check_numeric_generation(n, d, generators, D_max=None, ...)`
  - Runs the truncated numeric closure algorithm and returns `NumericResult`.
- `lib.numeric.check_full_lie_generation(n, generators, required_degree=2, ...)`
  - Runs the project criterion for full Lie-generation and returns `FullLieGenerationCheck`.
- `lib.numeric.evaluate_full_lie_generation(result, required_degree=2)`
  - Applies the criterion to an existing `NumericResult`.
- `lib.numeric.run_hypothesis(...)`
  - Multi-trial workflow for a fixed or random second generator.

## Generator specs

- `lib.numeric.parse_generator_spec("dx + x*dy", n=2)`
- `lib.numeric.make_monomial_D(n, a, b, target_var=0)`
- `lib.numeric.make_general_monomial_generator(n, exponents, target_var=0)`
- `lib.numeric.make_random_generator(n, max_degree, rng=None, sparsity=0.5)`
- `lib.numeric.describe_generator(generator, n)`

Representation:

```python
{axis: {alpha_tuple: coefficient}}
```

Example:

```python
from lib.numeric import check_full_lie_generation, parse_generator_spec

g1 = parse_generator_spec("x^3*dx + y^3*dy", 2)
g2 = parse_generator_spec("dx + x*dy", 2)
check = check_full_lie_generation(n=2, generators=[g1, g2], D_max=3)
print(check.criterion_result.satisfied)
print(check.numeric_result.summary())
```

## Structure constants

- `lib.numeric.StructureConstantCache(n, D_max)`
  - Reuse this across repeated calls, especially in enumeration workloads.

Example:

```python
from lib.numeric import StructureConstantCache, check_numeric_generation, parse_generator_spec

cache = StructureConstantCache(n=2, D_max=3)
g1 = parse_generator_spec("x^3*dx + y^3*dy", 2)
g2 = parse_generator_spec("dx + x*dy", 2)
result = check_numeric_generation(
    n=2,
    d=2,
    generators=[g1, g2],
    D_max=3,
    structure_constants=cache,
    verbose=False,
)
```

## CLI modes

- `python numeric_check.py --mode criterion --n 2 --g1 "..." --g2 "..."`
- `python numeric_check.py --mode hypo --n 2 --d 2 --g1 "..." --g2 "..."`
- `python numeric_check.py --mode beldiev --n 2 --d 7`
- `python numeric_enumerate.py --coeffs 1,-1 --max-terms 2 --limit 20`

## File map

- `lib/numeric/solver.py`: closure algorithm only
- `lib/numeric/results.py`: result dataclasses
- `lib/numeric/specs.py`: parser, formatting, generator builders
- `lib/numeric/recipes.py`: named generator families
- `lib/numeric/criteria.py`: degree-2 full-generation criterion
- `lib/numeric/workflows.py`: high-level workflows for agents and CLI wrappers
- `lib/numeric/enumeration.py`: reusable finite-search helpers
