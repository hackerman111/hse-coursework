"""Тесты CLI symbolic_check.py."""

import subprocess
import sys
import unittest


class TestSymbolicCheckHypo(unittest.TestCase):
    """Тест режима hypo через subprocess (требует SageMath)."""

    def _run(self, *args):
        result = subprocess.run(
            [sys.executable, "symbolic_check.py", *args],
            capture_output=True,
            text=True,
            timeout=120,
        )
        return result

    def test_hypo_with_known_generators_n2(self):
        """dx и x*dy порождают W_2 — солвер должен вернуть success."""
        result = self._run(
            "--mode", "hypo",
            "--n", "2",
            "--g1", "dx",
            "--g2", "x*dy",
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("PASS", result.stdout)

    def test_hypo_with_incomplete_generators(self):
        """x*dx и y*dy не порождают W_2 — солвер должен вернуть fail."""
        result = self._run(
            "--mode", "hypo",
            "--n", "2",
            "--g1", "x*dx",
            "--g2", "y*dy",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("FAIL", result.stdout)

    def test_hypo_random_g2(self):
        """Без --g2 генерируем случайный — скрипт не должен падать."""
        result = self._run(
            "--mode", "hypo",
            "--n", "2",
            "--g1", "dx",
            "--seed", "42",
            "--trials", "1",
        )
        # Может быть 0 или 1, главное — не crash
        self.assertIn("Symbolic", result.stdout)


class TestSymbolicCheckCriterion(unittest.TestCase):
    """Тест режима criterion через subprocess."""

    def _run(self, *args):
        result = subprocess.run(
            [sys.executable, "symbolic_check.py", *args],
            capture_output=True,
            text=True,
            timeout=120,
        )
        return result

    def test_criterion_with_full_generators(self):
        """Генераторы Бельдиева для W_2 порождают все дифференцирования степени <= 2."""
        result = self._run(
            "--mode", "criterion",
            "--n", "2",
            "--g1", "dy",
            "--g2", "y^8*dx + x^4*y^4*dy",
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("PASS", result.stdout)

    def test_criterion_with_incomplete_generators(self):
        """x*dx и y*dy не порождают все степени <= 2."""
        result = self._run(
            "--mode", "criterion",
            "--n", "2",
            "--g1", "x*dx",
            "--g2", "y*dy",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("FAIL", result.stdout)


if __name__ == "__main__":
    unittest.main()
