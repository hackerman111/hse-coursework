# hse-coursework

Исследовательский код для проверки конечной порождённости алгебр Ли полиномиальных векторных полей. Репозиторий содержит два вычислительных слоя:

- `numeric` — быстрые усечённые проверки в `W_n^(<= d)` на NumPy без зависимости от SageMath.
- `symbolic` — точные символьные проверки через `LieBasisSolver` и SageMath.

Основная идея проекта: задавать генераторы полиномиальных дифференцирований, запускать проверку замыкания по скобке Ли и анализировать, порождают ли они нужные слои `W_n`.

## Возможности

- Numeric solver для проверки порождения до степени `d`.
- Project criterion: если порождены все дифференцирования до степени `2`, набор считается порождающим всю алгебру Ли.
- Symbolic solver для точной проверки через базис Ли.
- Единый верхнеуровневый Python API через `lib`.
- CLI-скрипты для быстрых проверок и конечного перебора кандидатов.
- Готовые семейства генераторов Бельдиева и Андриста.

## Требования

### Numeric layer

- Python 3
- `numpy`
- `pytest` для тестов

### Symbolic layer

- Всё из numeric layer
- SageMath

Numeric API можно использовать без SageMath. Symbolic API и `symbolic_check.py` требуют установленный SageMath.

## Быстрый старт

Работа ведётся из корня репозитория.

### Numeric CLI

```bash
python numeric_check.py --mode beldiev --n 2 --d 7
python numeric_check.py --mode hypo --n 2 --d 2 --g1 "x^5*y^5*dx" --g2 "dx + x*dy"
python numeric_check.py --mode criterion --n 2 --g1 "x^5*y^5*dx" --g2 "dx + x*dy"
python numeric_enumerate.py --coeffs 1,-1 --max-terms 2 --limit 20
```

### Symbolic CLI

```bash
python symbolic_check.py --mode hypo --n 2 --g1 "dx" --g2 "x*dy"
python symbolic_check.py --mode criterion --n 2 --g1 "dy" --g2 "y^8*dx + x^4*y^4*dy"
```

Если `sage.all` доступен только через `sage -python`, запускайте так:

```bash
sage -python symbolic_check.py --mode hypo --n 2 --g1 "dx" --g2 "x*dy"
```

## Форматы данных

### Canonical spec

Каноническое представление библиотеки:

```python
DerivationSpec(n=2, terms={0: {(1, 0): 1.0}, 1: {(0, 0): 1.0}})
```

Это означает:

```text
x * d/dz0 + 1 * d/dz1
```

### Numeric generator

Разреженный numeric-формат:

```python
{0: {(2, 1): 1.0}, 1: {(0, 0): 3.0}}
```

Здесь ключ верхнего уровня — индекс компоненты `d/dz_i`, а вложенный словарь задаёт моном `z^alpha -> coeff`.

### Строковый generator-spec

Поддерживаются записи вида:

```text
dx
x*dy
x^5*y^5*dx
2*z0^3*d1 - 0.5*z1*dz0
```

Допустимые обозначения переменных: `x,y,z,u,v,w` и `z0,z1,...`.

## Python API

### Верхний уровень `lib`

`lib` теперь служит единой публичной точкой входа:

```python
from lib import (
    check_generation,
    check_numeric_generation,
    check_full_lie_generation,
    parse_generator_spec,
    make_beldiev_generators,
)
```

Через него же доступны lazy-loaded symbolic entry points:

```python
from lib import check_symbolic_generation, check_symbolic_degree2_criterion
```

### Numeric API

Основные функции:

- `check_numeric_generation(n, d, generators, ...) -> NumericResult`
- `check_full_lie_generation(n, generators, ...) -> FullLieGenerationCheck`
- `evaluate_full_lie_generation(result, required_degree=2) -> LieGenerationCriterionResult`
- `run_beldiev(n, d, ...) -> NumericResult`
- `run_hypothesis(n, d, first_gen, ...) -> HypothesisRunResult`

Пример:

```python
from lib import check_numeric_generation, parse_generator_spec

g1 = parse_generator_spec("x^3*dx + y^3*dy", 2)
g2 = parse_generator_spec("dx + x*dy", 2)

result = check_numeric_generation(
    n=2,
    d=2,
    generators=[g1, g2],
    D_max=3,
    verbose=False,
)

print(result.success)
print(result.summary())
```

Если нужен переиспользуемый кеш структурных констант:

```python
from lib import StructureConstantCache, check_numeric_generation, parse_generator_spec

cache = StructureConstantCache(n=2, D_max=3)
g1 = parse_generator_spec("x^3*dx + y^3*dy", 2)
g2 = parse_generator_spec("dx + x*dy", 2)

result = check_numeric_generation(
    n=2,
    d=2,
    generators=[g1, g2],
    D_max=3,
    verbose=False,
    structure_constants=cache,
)
```

### Symbolic API

Library-level symbolic entry points:

- `check_symbolic_generation(generators, max_iter=1000) -> SymbolicResult`
- `check_symbolic_degree2_criterion(generators, max_iter=1000) -> SymbolicCriterionResult`
- `make_symbolic_algebra(n)`
- `spec_to_symbolic_derivation(spec, algebra)`

`generators` можно передавать либо как `list[DerivationSpec]`, либо как `list[LieDerivation]`.

Пример:

```python
from lib import (
    check_symbolic_generation,
    check_symbolic_degree2_criterion,
    partial_spec,
)

generators = [partial_spec(2, axis=0), partial_spec(2, axis=1)]

full_result = check_symbolic_generation(generators, max_iter=100)
criterion_result = check_symbolic_degree2_criterion(generators, max_iter=100)

print(full_result.summary())
print(criterion_result.summary())
```

### Unified solver API

Если нужен один dispatching entry point:

```python
from lib import check_generation
from lib.generators.beldiev import beldiev_specs

result = check_generation(beldiev_specs(2), mode="numeric", d=2)
print(result.success)
```

Поддерживаемые режимы:

- `mode="numeric"`
- `mode="symbolic"`
- `mode="auto"`

## CLI-скрипты

### `numeric_check.py`

Режимы:

- `beldiev` — проверка генераторов Бельдиева
- `hypo` — проверка пары генераторов до степени `d`
- `criterion` — project criterion по степеням `<= 2`
- `dims` — проверка табличных размерностей

### `numeric_enumerate.py`

Конечный перебор кандидатов для

```text
G1 = P(x,y)dx + Q(x,y)dy
```

Ограничения задаются через:

- максимальные степени по `x` и `y`
- конечный набор коэффициентов
- ограничение на число ненулевых мономов

### `symbolic_check.py`

Режимы:

- `hypo` — полная символьная проверка пары генераторов
- `criterion` — ранняя остановка, когда найдены все цели степени `<= 2`

## Структура репозитория

```text
lib/
  core/           канонические типы и результаты
  backends/       конвертация между spec/numeric/sage-представлениями
  numeric/        numeric solver, workflows, criteria, enumeration
  solver/         symbolic solver и stop criteria
  generators/     спецификации Бельдиева и Андриста
  symbolic.py     high-level symbolic API

numeric_check.py
numeric_enumerate.py
symbolic_check.py
tests/
```

## Разработка

### Запуск тестов

Быстрый целевой прогон:

```bash
python -m pytest tests/test_public_api.py tests/test_api_doctests.py tests/test_solver_api.py tests/test_symbolic_check.py -q
```

Numeric regression:

```bash
python -m pytest tests/test_numeric.py tests/test_numeric_enumerate.py -q
```

Solver and parser:

```bash
python -m pytest tests/test_solver_api.py tests/test_io_parse.py tests/test_backends.py -q
```

### Важные замечания

- Numeric layer не требует SageMath.
- Symbolic layer требует SageMath и может использовать временный writable `DOT_SAGE`, если домашний каталог read-only.
- `symbolic_check.py` теперь является тонкой CLI-обёрткой над `lib.symbolic`, а не отдельной реализацией solver-логики.

## Что использовать в новом коде

- Для numeric-проверок используйте `lib` или `lib.numeric`.
- Для symbolic-проверок используйте `lib.symbolic` или lazy exports из `lib`.
- Для стабильного API избегайте обращения к приватным helper-функциям из CLI-скриптов.
