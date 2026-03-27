"""
Типы результатов проверки порождения алгебры Ли W_n.

Этот модуль является каноническим местом хранения DegreeStatus, NumericResult
и SymbolicResult. Модуль lib/numeric/results.py теперь является re-export shim.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class DegreeStatus:
    """Статус одного однородного слоя степени d."""

    degree: int
    rank: int
    target_dim: int

    @property
    def is_full(self) -> bool:
        return self.rank >= self.target_dim

    def __repr__(self) -> str:
        tag = "FULL" if self.is_full else f"{self.rank}/{self.target_dim}"
        return f"Degree {self.degree:3d}: {tag}"


@dataclass
class NumericResult:
    """Результат численной проверки порождения до степени d."""

    n: int
    d: int
    D_max: int
    success: bool
    degrees: list[DegreeStatus] = field(default_factory=list)
    iterations: int = 0
    elapsed_s: float = 0.0

    def degree_status(self, degree: int) -> DegreeStatus | None:
        """Статус конкретной степени, если она была проверена."""
        for status in self.degrees:
            if status.degree == degree:
                return status
        return None

    def covers_degree(self, degree: int) -> bool:
        """Была ли проверена данная степень."""
        return degree <= self.d and self.degree_status(degree) is not None

    def missing_degrees(self, max_degree: int | None = None) -> list[DegreeStatus]:
        """
        Список неполных степеней.

        При задании max_degree — только степени до него.
        """
        if max_degree is None:
            return [status for status in self.degrees if not status.is_full]
        return [
            status
            for status in self.degrees
            if status.degree <= max_degree and not status.is_full
        ]

    def is_full_up_to(self, degree: int) -> bool:
        """Все проверенные степени до degree полны."""
        if not self.covers_degree(degree):
            return False
        return not self.missing_degrees(max_degree=degree)

    def summary(self) -> str:
        lines = [
            "=== Numeric Lie Generation Check ===",
            f"n = {self.n}, d = {self.d}, D_max = {self.D_max}",
            f"Итераций: {self.iterations}, Время: {self.elapsed_s:.3f}s",
            (
                "Результат: PASS — все степени ≤ d полны"
                if self.success
                else "Результат: FAIL — есть неполные степени"
            ),
            "",
        ]
        for status in self.degrees:
            lines.append(str(status))
        return "\n".join(lines)


@dataclass
class SymbolicResult:
    """Результат символической проверки порождения через LieBasisSolver."""

    success: bool
    iterations: int = 0
    basis_size: int = 0
    targets_found: dict[int, bool] = field(default_factory=dict)

    def summary(self) -> str:
        lines = [
            "=== Symbolic Lie Generation Check ===",
            f"Итераций: {self.iterations}, Базис: {self.basis_size} элементов",
            (
                "Результат: PASS — все ∂/∂z_i найдены"
                if self.success
                else "Результат: FAIL — не все целевые производные найдены"
            ),
        ]
        if self.targets_found:
            found = [i for i, v in self.targets_found.items() if v]
            missing = [i for i, v in self.targets_found.items() if not v]
            if found:
                lines.append(f"Найдено: ∂/∂z_{{{', '.join(map(str, found))}}}")
            if missing:
                lines.append(f"Не найдено: ∂/∂z_{{{', '.join(map(str, missing))}}}")
        return "\n".join(lines)
