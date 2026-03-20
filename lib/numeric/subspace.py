"""
Ортонормальный базис подпространства с инкрементальным обновлением.

Используем QR-декомпозицию для гарантированной численной устойчивости:
- Все добавления проходят через полный QR-пересчёт.
- Ранг определяется по диагонали R с относительным порогом.
"""

from __future__ import annotations

import numpy as np


class OrthonormalBasis:
    """
    Инкрементально поддерживаемый ортонормальный базис подпространства R^N.

    Строки матрицы `rows` образуют ортонормальный набор.
    """

    __slots__ = ("dim", "rows", "_rank", "_eps")

    def __init__(self, dim: int, eps: float = 1e-12) -> None:
        self.dim = dim
        self.rows = np.empty((0, dim), dtype=np.float64)
        self._rank = 0
        self._eps = eps

    # ------------------------------------------------------------------
    # Публичный интерфейс
    # ------------------------------------------------------------------

    @property
    def rank(self) -> int:
        return self._rank

    @property
    def all_rows(self) -> np.ndarray:
        return self.rows

    def is_full(self, target_dim: int) -> bool:
        return self._rank >= target_dim

    # ------------------------------------------------------------------
    # Добавление векторов
    # ------------------------------------------------------------------

    def try_add(self, v: np.ndarray) -> bool:
        """
        Попробовать добавить вектор v. Возвращает True если он линейно
        независим от текущих строк (ранг вырос).
        """
        if self._rank >= self.dim:
            return False
        v = np.asarray(v, dtype=np.float64).ravel()
        return self.try_add_many(v.reshape(1, -1)) > 0

    def try_add_many(self, M: np.ndarray) -> int:
        """
        Добавить строки из матрицы M, возвращает число реально добавленных
        (линейно независимых) строк.

        Всегда использует стабильный QR-пересчёт.
        """
        M = np.asarray(M, dtype=np.float64)
        if M.ndim == 1:
            M = M.reshape(1, -1)
        if M.shape[0] == 0 or self._rank >= self.dim:
            return 0
        if M.shape[1] != self.dim:
            raise ValueError(
                f"expected row dimension {self.dim}, got {M.shape[1]}"
            )

        # Собираем текущие строки + новые кандидаты
        if self._rank > 0:
            combined = np.vstack([self.rows, M])
        else:
            combined = M

        # QR с перестановкой столбцов для устойчивого ранга
        # scipy.linalg.qr с pivoting даёт лучшую оценку ранга,
        # но для переносимости используем numpy SVD (золотой стандарт для ранга)
        old_rank = self._rank
        new_rank = self._compute_rank_and_basis(combined)

        return max(0, new_rank - old_rank)

    def _compute_rank_and_basis(self, combined: np.ndarray) -> int:
        """
        Вычислить ранг объединённого набора строк и обновить
        ортонормальный базис через SVD.

        SVD — наиболее устойчивый способ определения числового ранга.
        """
        # SVD: combined^T = U @ diag(s) @ V^T
        # Ранг = число сингулярных значений > threshold
        U, s, Vt = np.linalg.svd(combined.T, full_matrices=False)

        # Relative threshold (как в numpy.linalg.matrix_rank)
        tol = max(combined.shape) * s[0] * np.finfo(np.float64).eps if len(s) > 0 else 0.0
        tol = max(tol, self._eps)

        new_rank = int(np.sum(s > tol))
        new_rank = min(new_rank, self.dim)

        if new_rank > self._rank:
            # Ортонормальный базис = первые new_rank столбцов U, транспонированные
            self.rows = U[:, :new_rank].T.copy()
            self._rank = new_rank

        return new_rank
