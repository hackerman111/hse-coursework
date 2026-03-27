# План Повышения Модульности Библиотеки

## 1. Вывод после анализа кода

Сейчас библиотека уже состоит из трёх естественных направлений:

1. `lib.derivation`: символическое представление дифференцирований через Sage.
2. `lib.solver`: символический алгоритм замыкания по скобке Ли.
3. `lib.numeric`: численный алгоритм для проверки порождения до степени `d`.

Проблема не в отсутствии модулей, а в том, что **одни и те же математические сущности представлены по-разному** и поэтому каждый новый метод тянет за собой новую локальную реализацию.

Главный архитектурный дефект:

- математический объект "дифференцирование" живёт сразу в трёх формах:
  - `LieDerivation` поверх Sage;
  - словарь вида `{i: {alpha: coeff}}` в `numeric`;
  - строковые спецификации в `numeric_check.py` и `numeric_enumerate.py`.

Из-за этого библиотека сейчас не строится из кирпичиков, а из готовых вертикальных сценариев.

## 2. Где уже есть дублирование

### 2.1. Формулы генераторов дублируются

- `lib/generator/beldiev.py`
- `lib/numeric/solver.py::make_beldiev_generators`

Одна и та же математическая формула реализована дважды для двух представлений.

### 2.2. Конструкторы генераторов дублируются

- `numeric_check.py::make_monomial_D`
- `numeric_check.py::make_general_monomial_generator`
- `numeric_check.py::make_random_E`
- `lib/numeric/solver.py::make_monomial_generator`
- `lib/numeric/solver.py::make_random_generator`

Это должно быть одним слоем builders, а не кодом в CLI и solver.

### 2.3. Разбор и форматирование генераторов живут вне библиотеки

- `numeric_check.py::parse_generator_spec`
- `numeric_check.py::_describe_generator`
- `numeric_check.py::_describe_random_generator`
- `numeric_enumerate.py::candidate_to_spec`

То есть важная часть пользовательского API сейчас вообще не библиотечная.

### 2.4. Поведение `LieDerivation` и `QuotientLieDerivation` частично дублируется

Методы:

- `leading_term`
- `degree`

различаются только тем, как извлекается коэффициент из quotient ring.

Это признак отсутствия отдельного слоя "полиномиальные операции / адаптер коэффициентов".

### 2.5. Solver и numeric solver решают похожую задачу, но без общего контракта

Оба алгоритма делают одно и то же на разном backend:

- строят замыкание по скобке Ли;
- поддерживают базис подпространства;
- останавливаются по критерию полноты;
- формируют отчёт.

Но общий интерфейс у них отсутствует, поэтому логика повторяется на уровне идей, а не кода.

### 2.6. `lib/numeric/solver.py` перегружен ответственностями

В одном файле сейчас смешаны:

- модель результата;
- основной алгоритм;
- логирование прогресса;
- декомпозиция генератора по степеням;
- фабрики генераторов.

Это делает расширение дорогим.

## 3. Цель переработки

Нужна архитектура, в которой:

1. есть **одно каноническое представление** дифференцирования;
2. все готовые методы собираются из малых примитивов;
3. symbolic и numeric слои используют **одни и те же рецепты**, но разные backend-адаптеры;
4. CLI только вызывает библиотеку и ничего не "знает" о внутренней математике;
5. пользователь-математик может брать кирпичики и строить свои методы без копирования кода.

## 4. Базовый архитектурный принцип

Нужно разделить библиотеку на 5 слоёв.

### 4.1. Слой модели

Единая абстракция: `DerivationSpec`.

Предлагаемое каноническое представление:

```python
@dataclass(frozen=True)
class DerivationSpec:
    n: int
    terms: dict[int, dict[tuple[int, ...], Scalar]]

    def __post_init__(self):
        # Валидация: target_var в диапазоне 0..n-1
        for target_var in self.terms:
            if not (0 <= target_var < self.n):
                raise ValueError(f"target_var={target_var} вне диапазона 0..{self.n - 1}")
            for alpha, coeff in self.terms[target_var].items():
                if len(alpha) != self.n:
                    raise ValueError(f"len(alpha)={len(alpha)} != n={self.n}")
        # Нормализация: убрать нулевые коэффициенты
        normalized = {}
        for target_var, mono_dict in self.terms.items():
            cleaned = {a: c for a, c in mono_dict.items() if abs(c) > 1e-15}
            if cleaned:
                normalized[target_var] = cleaned
        object.__setattr__(self, "terms", normalized)
```

Смысл:

- ключ верхнего уровня: номер компоненты `∂_{z_i}`;
- внутренний словарь: моном `z^alpha -> coeff`.

Это уже фактически используется в numeric-части, значит перенос дешёвый.

Важно:

- `DerivationSpec` не зависит ни от Sage, ни от `numpy`;
- все theorem-generators должны сначала строить именно `DerivationSpec`;
- только потом спецификация преобразуется в Sage-объект или в плотные численные векторы;
- `__post_init__` гарантирует инварианты при создании: валидность `target_var`, длина мультииндексов `len(alpha) == n`, отсутствие нулевых коэффициентов.

### 4.2. Слой кирпичиков

Это минимальные математические операции, из которых собираются более сложные методы.

Принцип: начинаем с минимального набора, который покрывает **текущие** use-cases (рецепты Бельдиева/Андриста, CLI-парсинг, numeric solver). Остальные кирпичики добавляются по мере появления реальных потребностей (YAGNI).

**Нужны на первом этапе** (есть текущие use-cases):

- создать мономиальное поле (`monomial_spec`) — для рецептов генераторов;
- создать частную производную (`partial_spec`) — для U = ∂_{z_n};
- сложить поля (`combine_specs`) — для сборки V из слагаемых;
- вычислить степень (`spec_degree`) — для определения D_max;
- случайный генератор (`random_spec`) — для hypothesis testing.

**Добавить позже** (нет текущих use-cases, но могут понадобиться):

- создать нулевое поле;
- умножить поле на скаляр;
- выделить однородную компоненту;
- разложить поле по степеням;
- вычислить старший член;
- итерировать `ad_X`;
- строить линейную оболочку и редукцию к базису.

### 4.3. Слой backend-адаптеров

Два адаптера:

1. `SageBackend`
2. `NumericBackend`

Они должны реализовывать один набор преобразований:

- `from_spec(...)`
- `to_spec(...)`
- `bracket(...)`
- `degree(...)`
- `leading_term(...)`
- `components(...)`

Принципиальное решение: **`DerivationSpec` — это data transfer object, а не вычислительный объект.** Скобка Ли вычисляется бэкендом, а не на уровне spec. Причины:

1. Символическая скобка сейчас вычисляется через Sage, numeric — через structure constants + матрицы. Функция `lie_bracket_spec` на чистом Python потребовала бы **третьей реализации** скобки (мини-CAS для полиномиальной арифметики), что противоречит цели убрать дублирование.
2. Если в будущем понадобится «чистый Python backend» без Sage, его можно добавить как третий бэкенд, а не как метод spec.

Важно: не надо делать один "магический" универсальный solver для всего. Нужно сделать **общую модель и общие примитивы**, а движки оставить раздельными.

### 4.4. Слой алгоритмов

Здесь лежат уже не конкретные теоремы, а универсальные процедуры:

- замыкание по скобке Ли;
- поддержка базиса;
- проверка полноты;
- перечисление кандидатов;
- случайная генерация;
- вычисление структурных констант;
- работа с градуировкой.

### 4.5. Слой готовых рецептов

Это уже библиотека методов "из статьи":

- генераторы Бельдиева;
- генераторы Андриста;
- якобиановы дифференцирования;
- Вайтценбёковы дифференцирования;
- пользовательские конструктивные рецепты.

Каждый такой рецепт должен быть просто функцией, которая **собирает `DerivationSpec` из кирпичиков**.

## 5. Целевая структура пакетов

Текущая библиотека — ~1500 строк кода + ~600 строк CLI. Для этого масштаба структура должна быть компактной: каждый модуль содержит реальную логику, а не обёртку вокруг одной функции. Правило: если модуль вырос до 300+ строк — разделить в тот момент.

Предлагаемая раскладка (12 модулей, 6 пакетов):

```text
lib/
  core/
    spec.py            # DerivationSpec + валидация + кирпичики (monomial_spec, partial_spec, combine_specs, spec_degree, random_spec)
    results.py         # DegreeStatus, NumericResult, SymbolicResult
  backends/
    sage.py            # to_sage, from_sage + LieDerivation + QuotientLieDerivation
    numeric.py         # to_numeric_components, _decompose_generator
  solver/
    symbolic.py        # LieBasisSolver (бывший solver/basis.py + operations.py)
    numeric.py         # NumericLieGenerationChecker
    enumeration.py     # BasisTerm, iter_candidates, build_basis_terms, estimate_candidate_count
  generators/
    beldiev.py         # beldiev_specs() + beldiev() wrapper
    andrist.py         # andrist_specs() + andrist() wrapper
  io/
    parse.py           # parse_generator_spec, spec_to_string, _describe_generator
  utils/
    logger.py
```

Что убрано по сравнению с исходным планом и почему:

- `core/grading.py`, `core/terms.py` — в текущем коде это 5–10 строк, не оправдывают отдельных файлов; влить в `spec.py`.
- `builders/basic.py`, `builders/families.py`, `builders/recipes.py` — кирпичики и рецепты естественнее живут в `core/spec.py` и `generators/`; отдельный пакет `builders` избыточен.
- `algorithms/basis_symbolic.py` + `algorithms/symbolic_closure.py` — это один `LieBasisSolver`, нет смысла разделять на два файла. То же для numeric.
- `generators/jacobian.py`, `generators/weitzenbock.py` — это конструкторы из `lib/derivation/constructors.py`, каждый 10 строк; оставить в `backends/sage.py`.
- `io/format.py` — форматирование генераторов (~20 строк) влезает в `io/parse.py`.

Принцип раскладки:

- `core` ничего не знает о Sage и CLI;
- `backends` переводит spec в конкретное представление;
- `solver` содержит алгоритмы замыкания и перечисления;
- `generators` хранит готовые рецепты;
- `io` отвечает за строковые спецификации и форматирование.

## 6. Какие кирпичики надо сделать публичными

Ниже список функций, которые после переработки имеет смысл дать пользователю.

### 6.1. Низкоуровневые кирпичики

Функции, из которых математик собирает свои методы. Разделены на два эшелона.

**Первый эшелон** (реализовать сразу — есть текущие use-cases):

```python
def monomial_spec(n: int, axis: int, alpha: tuple[int, ...], coeff: Scalar = 1) -> DerivationSpec
def partial_spec(n: int, axis: int) -> DerivationSpec
def combine_specs(*specs: DerivationSpec) -> DerivationSpec
def spec_degree(spec: DerivationSpec) -> int
def random_spec(n: int, max_degree: int, sparsity: float = 0.5, seed: int | None = None) -> DerivationSpec
```

**Второй эшелон** (добавить когда появятся use-cases):

```python
def zero_spec(n: int) -> DerivationSpec
def scale_spec(spec: DerivationSpec, scalar: Scalar) -> DerivationSpec
def normalize_spec(spec: DerivationSpec, eps: float = 1e-15) -> DerivationSpec
def homogeneous_component(spec: DerivationSpec, degree: int) -> DerivationSpec
def homogeneous_components(spec: DerivationSpec) -> dict[int, DerivationSpec]
def leading_term_spec(spec: DerivationSpec, order: str = "lex") -> LeadingTerm | None
```

Скобка Ли (`lie_bracket_spec`) и итерирование `ad_X` **не входят** в API spec-уровня.
Скобка вычисляется бэкендом (см. §4.3). Если в будущем понадобится чистый Python backend — он добавляется как `PurePythonBackend`, а не как метод spec.

### 6.2. Конструкторы поверх кирпичиков

```python
def from_mapping(algebra, mapping) -> LieDerivation
def from_linear(algebra, matrix) -> LieDerivation
def from_weitzenbock(algebra, matrix) -> LieDerivation
def from_jacobian(algebra, polynomials) -> LieDerivation

def spec_from_mapping(n: int, mapping: dict[int, dict[tuple[int, ...], Scalar]]) -> DerivationSpec
def spec_from_string(text: str, n: int) -> DerivationSpec
def spec_to_string(spec: DerivationSpec, variables: list[str] | None = None) -> str
def random_spec(n: int, max_degree: int, sparsity: float = 0.5, seed: int | None = None) -> DerivationSpec
```

Смысл:

- один и тот же пользовательский генератор должен легко попадать и в symbolic, и в numeric pipeline;
- строковый парсер должен быть библиотечной функцией, а не CLI-скриптом.

### 6.3. Функции для работы с backend

```python
def to_sage(spec: DerivationSpec, algebra) -> LieDerivation
def from_sage(derivation) -> DerivationSpec
def to_numeric_components(spec: DerivationSpec, D_max: int) -> dict[int, np.ndarray]
def from_numeric_components(n: int, components: dict[int, np.ndarray]) -> DerivationSpec
```

Это убирает пропасть между `LieDerivation` и numeric dict-ами.

### 6.4. Функции алгоритмического уровня

```python
def symbolic_closure(generators: list[LieDerivation], max_iter: int = 1000) -> SymbolicClosureResult
def numeric_closure(specs: list[DerivationSpec], n: int, d: int, D_max: int | None = None) -> NumericResult
def check_generation(
    generators,
    *,
    mode: str = "symbolic",
    d: int | None = None,
    D_max: int | None = None,
    max_iter: int = 1000,
)
def enumerate_candidates(basis_terms, coeffs, min_terms: int, max_terms: int)
def estimate_search_space(basis_size: int, coeff_count: int, min_terms: int, max_terms: int) -> int
```

Ключевая идея:

- `check_generation(...)` становится единой точкой входа;
- старый `generator.check(...)` превращается в thin wrapper.

### 6.5. Готовые математические рецепты

```python
def beldiev_specs(n: int) -> list[DerivationSpec]
def andrist_specs(n: int) -> list[DerivationSpec]
def beldiev(algebra) -> list[LieDerivation]
def andrist(algebra) -> list[LieDerivation]
```

Здесь рецепт один, а способов материализации два:

- либо получить `DerivationSpec`;
- либо сразу Sage-объекты.

Именно так убирается текущее дублирование Бельдиева.

## 7. Как математик будет из этого собирать свои методы

После рефакторинга пользователь сможет собирать конструкции так:

```python
u = partial_spec(n=2, axis=1)
v = combine_specs(
    monomial_spec(2, axis=0, alpha=(0, 8)),
    monomial_spec(2, axis=1, alpha=(4, 4)),
)

result = numeric_closure([u, v], n=2, d=7, D_max=7)
```

Или так:

```python
g1 = spec_from_string("y^3*dx + x^3*dy", n=2)
g2 = spec_from_string("dx + x*dy", n=2)

result = check_generation([g1, g2], mode="numeric", d=2, D_max=3)
```

Или так:

```python
ring = PolynomialRing(QQ, "z1, z2, z3")
gens = beldiev(ring)
result = check_generation(gens, mode="symbolic", max_iter=500)
```

То есть библиотека станет не "набором готовых скриптов", а конструктором.

## 8. Что именно надо вынести из текущих файлов

### 8.1. Из `numeric_check.py`

В библиотеку (с указанием транзитивных зависимостей):

- `parse_generator_spec` → в `lib/io/parse.py`
  - тянет за собой: `_parse_derivative_factor`, `_parse_variable_factor`, `_alias_index`, `_validate_axis`, константу `VARIABLE_ALIASES` (5 функций + 1 константа)
- `_describe_generator`, `_describe_random_generator` → в `lib/io/parse.py` (formatting)
  - тянет за собой: `_generator_term_count`, `_generator_max_degree`
- `make_monomial_D`, `make_general_monomial_generator` → заменяются на `monomial_spec` в `lib/core/spec.py`
- `make_random_E` → заменяется на `random_spec` в `lib/core/spec.py`
- `check_dimensions` → в `lib/core/results.py` или убрать (это тестовая утилита)

В CLI должны остаться только:

- `parse_args`
- выбор режима;
- печать результата;
- вызов библиотечных функций.

### 8.2. Из `numeric_enumerate.py`

В библиотеку:

- `BasisTerm`
- `build_basis_terms`
- `parse_coefficients`
- `estimate_candidate_count`
- `build_generator`
- `candidate_to_spec`
- `iter_candidates`

CLI должен стать thin wrapper над `lib.algorithms.enumeration`.

### 8.3. Из `lib/numeric/solver.py`

В отдельные модули:

- `DegreeStatus`, `NumericResult` -> `core/results.py`
- `_decompose_generator` -> `builders/` или `backends/numeric_backend.py`
- `make_beldiev_generators`, `make_monomial_generator`, `make_random_generator` -> `builders/families.py` и `generators/beldiev.py`
- логика статусов -> отдельный reporter/helper

Сам `solver.py` должен содержать только алгоритм замыкания.

### 8.4. Из `lib/derivation/base.py` и `lib/derivation/quotient.py`

`QuotientLieDerivation` используется только в `tests/manual_test_quotient/` — это не production-код. Дублирование `leading_term` и `degree` реально, но масштаб не оправдывает отдельный Protocol.

**Текущее решение (минимальное):** добавить в `LieDerivation` overridable метод `_extract_poly(element)`, который по умолчанию возвращает `element` as-is, а в `QuotientLieDerivation` делает `.lift()`. Это устраняет дублирование одной строкой вместо нового типа.

**Отложенное решение:** если quotient-логика станет production-кодом, тогда ввести `PolynomialOps` Protocol. Не раньше.

## 9. Как избежать дублирования в будущем

Нужно принять три жёстких правила.

### Правило 1. Любая новая математическая формула сначала пишется как `DerivationSpec`

Не как Sage-объект. Не как `numpy`-вектор. Не как CLI-строка.

Только после этого добавляются адаптеры:

- `to_sage`
- `to_numeric_components`

### Правило 2. Любой новый CLI использует только библиотечный API

Если в CLI появился новый математический код, значит архитектура снова начала расползаться.

### Правило 3. Теорема, алгоритм и представление должны жить в разных местах

Пример:

- теорема Бельдиева -> `generators/beldiev.py`
- замыкание по скобке -> `algorithms/`
- преобразование в Sage -> `backends/sage_backend.py`

Это минимизирует пересечения.

## 10. Что оставить совместимым

Чтобы не ломать существующий код, старые функции можно сохранить как обёртки:

```python
get_Beldiev(algebra)  -> return beldiev(algebra)
get_Andristy(algebra) -> return andrist(algebra)
check(generators, max_iter=1000) -> return check_generation(..., mode="symbolic")
Derivation(algebra, mapping) -> return from_mapping(algebra, mapping)
```

То есть публичный API можно расширить без резкого breaking change.

## 11. Минимальный план рефакторинга

### Этап 0. Тестовая инфраструктура (ПЕРЕД любым рефакторингом)

Рефакторинг такого масштаба без тестовой стратегии — рецепт для регрессий.

1. **Зафиксировать baseline:** убедиться, что `pytest tests/ -v` проходит полностью. Каждый последующий этап должен заканчиваться прогоном всех тестов без регрессий.

2. **Round-trip тесты для DerivationSpec** (добавить в этап 1):
   - `spec → to_sage → from_sage → spec` — проверяет корректность Sage-адаптера;
   - `spec → to_numeric_components → from_numeric_components → spec` — проверяет numeric-адаптер;
   - Оба round-trip должны давать идентичный spec (с точностью до eps для float).

3. **Round-trip тесты для parse/format:**
   - `parse_generator_spec(spec_to_string(spec)) == spec` — проверяет, что парсер и форматтер согласованы.
   - Сейчас парсер тестируется только косвенно через CLI.

4. **Кросс-бэкендный математический тест:**
   - Для набора тестовых генераторов: вычислить скобку Ли через Sage и через numeric, сравнить результаты.
   - Это самый ценный тест — он гарантирует, что два бэкенда реализуют одну математику.

### Этап 1. Ввести каноническую модель

- добавить `DerivationSpec` с валидацией в `__post_init__`;
- добавить `spec_from_string`, `spec_to_string`;
- добавить `monomial_spec`, `partial_spec`, `combine_specs`, `spec_degree`, `random_spec`;
- написать round-trip тесты (из этапа 0).

### Этап 2. Перенести theorem-generators на specs

- `beldiev_specs`
- `andrist_specs`

Старые `get_Beldiev`, `get_Andristy` должны стать thin wrapper.

### Этап 3. Вынести builders из CLI и numeric solver

- случайные генераторы;
- мономиальные генераторы;
- перечисление кандидатов;
- парсер и форматтер генераторов (с транзитивными зависимостями, см. §8.1).

### Этап 4. Разделить алгоритмы и результаты

- вынести result/reporting;
- очистить `lib/numeric/solver.py`;
- сделать единый `check_generation`;
- добавить возможность переиспользовать structure constants между вызовами checker (см. §14).

### Этап 5. Упростить symbolic layer

- добавить `_extract_poly` в `LieDerivation` (overridable);
- переписать `QuotientLieDerivation` через `_extract_poly` вместо дублирования `leading_term` и `degree`.

## 12. Итоговая архитектурная формула

После переработки библиотека должна читатьcя так:

```text
Spec -> Builder -> Backend -> Algorithm -> Result
```

Где:

- `Spec` — единое математическое описание;
- `Builder` — кирпичики и рецепты;
- `Backend` — Sage или numeric;
- `Algorithm` — symbolic closure, numeric closure, enumeration;
- `Result` — единый формат ответа и логирования.

Это даст три ключевых эффекта:

1. математик сможет собирать новые методы из маленьких функций;
2. формулы Бельдиева, Андриста и будущих конструкций больше не будут дублироваться;
3. numeric и symbolic части перестанут расходиться как две разные библиотеки.

## 13. Короткий список того, что должен вызывать человек после переработки

Минимальный публичный API, разделённый на два эшелона.

**Первый эшелон** (реализовать сразу):

```python
# Кирпичики (core/spec.py)
monomial_spec
partial_spec
combine_specs
spec_degree
random_spec

# I/O (io/parse.py)
spec_from_string
spec_to_string

# Бэкенды (backends/)
to_sage
from_sage
to_numeric_components

# Алгоритмы (solver/)
symbolic_closure
numeric_closure
check_generation

# Рецепты (generators/)
beldiev_specs
andrist_specs
beldiev
andrist
```

**Второй эшелон** (добавить по мере появления use-cases):

```python
zero_spec
scale_spec
homogeneous_component
homogeneous_components
leading_term_spec
from_numeric_components
```

Скобка Ли и `iterate_ad` **не входят** в API spec-уровня — они живут в бэкендах.

## 14. Performance: кеширование structure constants

### Проблема

Самый тяжёлый этап в numeric solver — `precompute_structure_constants` в `__init__`. Для `n=2, D_max=10` это O(D_max²) пар, каждая с O(n · |multiindices|²) операциями.

В `numeric_enumerate.py` structure constants **одинаковые** для всех кандидатов (одни и те же `n`, `D_max`), но пересчитываются при каждом создании `NumericLieGenerationChecker`. Это чистый waste при тысячах вызовов.

### Решение

1. **Вынести structure constants в отдельный объект:**

```python
class StructureConstants:
    def __init__(self, n: int, D_max: int):
        self._cache = {}
        for p in range(-1, D_max + 1):
            for q in range(p, D_max + 1):
                r = p + q
                if -1 <= r <= D_max:
                    self._cache[(p, q)] = precompute_structure_constants(n, p, q)
```

2. **Принимать в `NumericLieGenerationChecker` как опциональный параметр:**

```python
def __init__(self, ..., precomputed_sc: StructureConstants | None = None):
    self._sc_cache = precomputed_sc or StructureConstants(n, D_max)
```

3. **В `numeric_enumerate.py` — создавать один раз перед циклом:**

```python
sc = StructureConstants(n=2, D_max=dmax)
for candidate in iter_candidates(...):
    checker = NumericLieGenerationChecker(..., precomputed_sc=sc)
```

### Когда реализовать

Этап 4 (разделение алгоритмов и результатов) — естественный момент для этого изменения, поскольку `solver.py` и так перестраивается.
