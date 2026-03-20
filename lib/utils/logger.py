"""
Переиспользуемый логгер для текстового прогресса библиотечных методов.

Класс отделяет форматирование человекочитаемых сообщений от конкретных
CLI-скриптов, чтобы будущие методы библиотеки могли использовать единый
механизм вывода без прямых вызовов print().
"""

from __future__ import annotations

import sys
import time
from typing import Callable, Optional, Sequence, TextIO


class LibraryLogger:
    """
    Простой потоковый логгер для консольных сценариев библиотеки.

    Он не привязан к конкретному модулю и подходит для переиспользования
    в любых методах, которым нужен человекочитаемый прогресс.
    """

    def __init__(self, enabled: bool = True, stream: TextIO | None = None) -> None:
        self.enabled = enabled
        self.stream = stream if stream is not None else sys.stdout

    def line(self, message: str = "", *, indent: int = 0) -> None:
        """Вывести строку или многострочный блок с отступом."""
        if not self.enabled:
            return

        prefix = " " * indent
        if not message:
            self.stream.write("\n")
        else:
            for raw_line in str(message).splitlines():
                self.stream.write(f"{prefix}{raw_line}\n")
        self.stream.flush()

    def info(self, message: str, *, indent: int = 2) -> None:
        """Вывести обычное информационное сообщение."""
        self.line(message, indent=indent)

    def kv(
        self,
        label: str,
        value: object,
        *,
        indent: int = 2,
        separator: str = " = ",
    ) -> None:
        """Вывести пару ключ-значение."""
        self.line(f"{label}{separator}{value}", indent=indent)

    def rule(self, label: str | None = None, *, indent: int = 2, char: str = "-") -> None:
        """Вывести разделитель секций."""
        text = char * 8 if label is None else f"{char * 4} {label} {char * 4}"
        self.line(text, indent=indent)

    def banner(self, title: str, details: Sequence[str] = ()) -> None:
        """Вывести короткий баннер из нескольких строк."""
        if not self.enabled:
            return

        content = [title, *details]
        width = max(len(line) for line in content)
        border = f"+-{'-' * width}-+"

        self.line(border, indent=0)
        for line in content:
            self.line(f"| {line.ljust(width)} |", indent=0)
        self.line(border, indent=0)
        self.line()

    def maybe_emit_periodic(
        self,
        *,
        next_log_at: float | None,
        interval_s: float | None,
        message_factory: Callable[[], str],
        now: float | None = None,
        surround_with_blank_lines: bool = True,
    ) -> Optional[float]:
        """
        Периодически выводить блок текста и вернуть время следующего лога.
        """
        if not self.enabled or next_log_at is None or interval_s is None:
            return next_log_at

        current_time = time.perf_counter() if now is None else now
        if current_time < next_log_at:
            return next_log_at

        if surround_with_blank_lines:
            self.line()
        self.line(message_factory(), indent=0)
        if surround_with_blank_lines:
            self.line()
        return current_time + interval_s
