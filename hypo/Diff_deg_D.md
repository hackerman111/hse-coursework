Вот полная формализация. Структурная схема взаимодействия фаз — сразу, затем псевдокод каждой фазы с точными определениями.---

## Алгоритм GENCHECK — полная формализация

### Обозначения и предварительные определения

Зафиксируем: $K$ — алгебраически замкнутое поле, $\mathrm{char}\,K = 0$, $R = K[z_1,\ldots,z_n]$, $W_n = \mathrm{Der}(R)$.

Степень поля $X = \sum_k f_k \partial_k$: $\deg X = \max_k \deg f_k$. Мономиальная база целевого пространства:
$$\mathcal{B}_d = \bigl\{z^\alpha \partial_k : |\alpha| \leq d,\; 1 \leq k \leq n\bigr\}, \quad N_d = |\mathcal{B}_d| = n\tbinom{n+d}{n}.$$

Рабочий базис при степени $D$: $\mathcal{B}_D$ аналогично, $N_D = n\binom{n+D}{n}$.

Представление элемента $X \in W_n^{\leq D}$ как вектора: $\mathbf{v}(X) \in K^{N_D}$ — координаты в $\mathcal{B}_D$.

`REF` — операция приведения к ступенчатому виду (row echelon form).

---

### Фаза 0: инициализация

**Вход:** $G = \{G_1,\ldots,G_m\} \subset W_n$, целевая степень $d$, размерность $n$.

```
D       ← max(d, max_{i} deg G_i)
N_d     ← n · C(n+d, n)        // целевая размерность
N_D     ← n · C(n+D, n)        // рабочая размерность

L       ← пустая матрица размера 0 × N_D   // ступенчатый базис подпространства L ⊆ W_n^{≤D}
M       ← ∅                                // множество мономов, заведомо ∈ L (моном-трекер)
P       ← ∅                                // множество i: ∂/∂z_i ∈ L
Q       ← min-heap, ключ = deg(·)          // приоритетная очередь

для каждого G_i:
    если REF_INSERT(L, v(G_i)):             // вставить, если линейно независим
        Q.push(G_i, key=deg(G_i))
        MONOMIAL_TRACK(G_i, M, P)
```

`REF_INSERT(L, v)` добавляет вектор $v$ в ступенчатую матрицу $L$, возвращает True, если $v \notin \mathrm{span}(L)$, и обновляет $L$. Стоимость: $O(N_D \cdot \mathrm{rank}\, L)$.

`MONOMIAL_TRACK(X, M, P)`: если $X = c \cdot z^\alpha \partial_k$ (моном), добавить $z^\alpha \partial_k$ в $M$; если $|\alpha| = 0$, добавить $k$ в $P$.

---

### Фаза 1: пилинг (extraction of $\partial/\partial z_i$)

Ключевое тождество: $[\partial/\partial z_i,\, X] = \sum_k \frac{\partial f_k}{\partial z_i} \partial_k$ для $X = \sum_k f_k \partial_k$.

```
changed ← True
пока changed:
    changed ← False

    для каждого X ∈ elements(Q) в порядке возрастания deg:
        для каждого i ∈ P:
            cur ← X
            пока deg(cur) > 0:
                cur ← [∂/∂z_i, cur]     // дифференцирование коэффициентов по z_i
                cur ← π_D(cur)           // усечение до степени D
                если cur = 0: break
                если REF_INSERT(L, v(cur)):
                    Q.push(cur, key=deg(cur))
                    MONOMIAL_TRACK(cur, M, P)
                    если P изменился: changed ← True
                    если SUFFICIENT_CONDITION(M, n):  // H2
                        вернуть True
    
    // Дополнительный источник: скобки пар генераторов
    для каждой пары (G_i, G_j), i < j:
        Y ← π_D([G_i, G_j])
        если REF_INSERT(L, v(Y)):
            Q.push(Y, key=deg(Y))
            MONOMIAL_TRACK(Y, M, P)
            // Применить пилинг к Y немедленно
            PEEL(Y, P, L, M, Q)
```

**Инвариант после Ф1:** $P = \{1,\ldots,n\}$ при «хороших» генераторах; алгоритм продолжает к Ф2 только если $|P| > 0$.

---

### Фаза 2: построение $\mathfrak{gl}_n$ (элементы степени 1)

Используется: $[z_i \partial_n,\, z_n \partial_j] = z_i \partial_j$ для $i \neq j$ (кросс-терм).

```
// Извлечь элементы W_n^(1) пилингом из более высоких степеней
для каждого X ∈ elements(L) с deg(X) ≥ 2:
    для каждого i ∈ P:
        Y ← [∂/∂z_i, X]
        если deg(Y) = 1 и REF_INSERT(L, v(Y)):
            Q.push(Y, key=1)
            MONOMIAL_TRACK(Y, M, P)

// Замкнуть W_n^(1) ≅ gl_n под скобкой (W_n^(1) — подалгебра Ли)
changed ← True
пока changed:
    changed ← False
    для каждых A, B ∈ elements(L) с deg(A) = deg(B) = 1:
        C ← [A, B]
        если deg(C) = 1 и REF_INSERT(L, v(C)):
            Q.push(C, key=1); changed ← True

// Кросс-термы
для i, j ∈ {1,…,n} с i ≠ j:
    если z_i·∂_n ∈ L и z_n·∂_j ∈ L:
        добавить z_i·∂_j в M    // [z_i∂_n, z_n∂_j] = z_i∂_j

// Проверка полноты gl_n
если dim(L ∩ W_n^(1)) < n²:
    перейти к Ф4   // Ф3 недоступна без полного gl_n

если SUFFICIENT_CONDITION(M, n): вернуть True  // H2
```

---

### Фаза 3: индуктивный подъём по степеням

Для каждой степени $p$ от 2 до $d$ последовательно выполняем четыре подфазы.

#### Ф3a: повышение степени (power raising)

Основная рекурсия (Бельдиев, шаг 5):

$$\bigl[z_k^s \partial_k,\; z_k^2 \partial_k\bigr] = (2-s)\, z_k^{s+1} \partial_k, \quad 2 - s \neq 0 \text{ при } s \geq 3.$$

```
для каждого k = 1,…,n:
    для s = 3, 4, …, p:
        если z_k^s ∂_k ∈ M и z_k^2 ∂_k ∈ M:
            new ← z_k^{s+1} ∂_k           // = [z_k^s ∂_k, z_k^2 ∂_k] / (2-s)
            добавить new в M
```

#### Ф3b: смена направления (direction change)

Тождество (Бельдиев, шаг 6):

$$\bigl[z_k^s \partial_k,\; z_k \partial_j\bigr] = z_k^s \partial_j.$$

```
для каждых k ≠ j:
    для s = 0,1,…,p:
        если z_k^s ∂_k ∈ M и z_k ∂_j ∈ M:
            добавить z_k^s ∂_j в M
```

#### Ф3c: умножение мономов (variable mixing)

Тождество (Бельдиев, шаг 7), где $p_k = \deg_{z_k}(f)$:

$$\bigl[z_\ell^s \partial_k,\; f \cdot z_k \partial_k\bigr] = (p_k + 1)\, z_\ell^s \cdot f\, \partial_k.$$

```
для каждого k = 1,…,n:
    для каждого монома f·∂_k ∈ M, зависящего от z_1,…,z_{ℓ-1}:
        для s = 1,…,p - deg(f):
            если z_ℓ^s ∂_k ∈ M и f·z_k ∂_k ∈ M:
                добавить z_ℓ^s · f · ∂_k в M   // = [z_ℓ^s ∂_k, f·z_k ∂_k] / (p_k+1)
```

#### Ф3d: действие $\mathfrak{gl}_n$ (заполнение $\mathfrak{gl}_n$-модуля)

Каждая однородная компонента $W_n^{(p)} \cong S^p(V) \otimes V^*$ — неприводимый $\mathfrak{gl}_n$-модуль. Одного элемента в «общем положении» достаточно.

```
для каждого A = z_i ∂_j ∈ L (степень 1):
    для каждого X ∈ elements(L) с deg(X) = p:
        Y ← [A, X]         // deg(Y) = p, поскольку deg(A) = 1
        если REF_INSERT(L, v(Y)):
            MONOMIAL_TRACK(Y, M, P)

// Послойная проверка (H3)
если dim(L ∩ W_n^(p)) = n · C(n+p-1, p):
    удалить W_n^(p) из активных элементов Q  // пропустить в дальнейшем

если SUFFICIENT_CONDITION(M, n): вернуть True
```

---

### Фаза 4: запасное насыщение (fallback BFS)

Запускается только если Ф3 не достигла $N_d$.

```
пока dim(L|_{≤d}) < N_d и Q не пуст:
    X ← Q.pop_min()     // элемент наименьшей степени
    для каждого Y ∈ elements(L):
        если (X, Y) ∉ visited:
            visited.add((X, Y))
            если deg(X) + deg(Y) - 1 ≤ D:
                Z ← π_D([X, Y])
                если REF_INSERT(L, v(Z)):
                    Q.push(Z, key=deg(Z))
                    MONOMIAL_TRACK(Z, M, P)
                    PEEL(Z, P, L, M, Q)   // немедленный пилинг
                    если SUFFICIENT_CONDITION(M, n): вернуть True

    если Q.size() стабилизировался (2 итерации без роста dim):
        вернуть False   // стабилизация — поколение не достигается
```

---

### Фаза 5: верификация

```
вернуть dim(L|_{≤d}) ≥ N_d
```

---

### Достаточное условие (эвристика H2)

```python
def SUFFICIENT_CONDITION(M, n):
    """
    Теорема: если ∀i: ∂/∂z_i ∈ L  и  z_i∂/∂z_j ∈ L  и  z_i²∂/∂z_i ∈ L,
    то L ⊇ W_n^{≤d} для любого d ≥ 0.
    Проверяем 2n + n² элементов вместо N_d.
    """
    return (all((1, k) in M         for k in range(n)) and   # ∂/∂z_k
            all((e_i, k) in M       for i in range(n)         # z_i ∂/∂z_j
                                    for k in range(n)) and
            all((2*e_k, k) in M     for k in range(n)))       # z_k² ∂/∂z_k
```

Проверка занимает $O(n^2)$ операций и выполняется после каждой вставки в $M$. При достижении этого условия возвращаем `True` немедленно, не вычисляя оставшиеся $N_d - (2n + n^2)$ мономов.

---

### Сводка сложности

| Фаза | Стоимость | Условие применимости |
|---|---|---|
| Ф0 | $O(m \cdot N_D^2)$ | всегда |
| Ф1 (пилинг) | $O(D \cdot n \cdot m \cdot N_D^2)$ | всегда |
| Ф2 ($\mathfrak{gl}_n$) | $O(n^4 \cdot N_D)$ | после Ф1 |
| Ф3a–3d | $O(d \cdot n^2 \cdot N_d)$ | при $\mathfrak{gl}_n \subset L$ |
| Ф4 (BFS) | $O(N_D^3)$ | запасная |
| H2 (проверка) | $O(n^2)$ | после каждой вставки |

Типичный случай (генераторы вида Андрист/Бельдиев): алгоритм завершается на Ф1–Ф3 с суммарной стоимостью $O(d \cdot n^2 \cdot N_d)$, значительно не достигая $O(N_D^3)$ общего BFS.
