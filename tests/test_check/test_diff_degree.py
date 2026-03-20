import io
import os
import pathlib
import subprocess
import sys
import unittest

from sage.all import PolynomialRing, QQ

from lib.derivation import Derivation
from lib.generator import check_diff_degree, get_Beldiev
from lib.solver import (
    DiffDegreeLieGenerationChecker,
    dim_Wn_coeff_k,
    dim_Wn_coeff_leq_d,
)
from lib.utils import LibraryLogger

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]


class DiffDegreeDimensionsTests(unittest.TestCase):
    def test_dimension_examples(self):
        expected = [
            (2, 5, 42),
            (2, 10, 132),
            (3, 5, 168),
            (3, 10, 858),
            (5, 5, 1260),
            (5, 10, 15015),
        ]

        for n, d, total_dim in expected:
            with self.subTest(n=n, d=d):
                self.assertEqual(dim_Wn_coeff_leq_d(n, d), total_dim)

        self.assertEqual(dim_Wn_coeff_k(2, 0), 2)
        self.assertEqual(dim_Wn_coeff_k(2, 2), 6)


class DiffDegreeCheckerTests(unittest.TestCase):
    def test_dx_and_x_squared_dx_generate_w1_up_to_degree_2(self):
        ring = PolynomialRing(QQ, "x")
        (x,) = ring.gens()
        dx = Derivation(ring, {x: 1})
        x2dx = Derivation(ring, {x: x**2})

        checker = DiffDegreeLieGenerationChecker(
            [dx, x2dx],
            2,
            verbose=False,
            max_iter=200,
        )
        result = checker.run()

        self.assertTrue(result.success)
        self.assertTrue(result.sufficient_condition)
        self.assertEqual(result.partial_axes, (0,))

    def test_single_linear_field_fails(self):
        ring = PolynomialRing(QQ, "x")
        (x,) = ring.gens()
        xdx = Derivation(ring, {x: x})

        result = check_diff_degree([xdx], 2, verbose=False, max_iter=100)

        self.assertFalse(result.success)
        self.assertEqual(result.partial_axes, ())

    def test_beldiev_generators_work_symbolically(self):
        ring = PolynomialRing(QQ, "x,y")
        generators = get_Beldiev(ring)

        result = check_diff_degree(
            generators,
            2,
            work_degree=8,
            verbose=False,
            max_iter=2000,
        )

        self.assertTrue(result.success)
        self.assertGreaterEqual(result.target_rank, 1)

    def test_checker_uses_library_logger(self):
        ring = PolynomialRing(QQ, "x")
        (x,) = ring.gens()
        dx = Derivation(ring, {x: 1})
        x2dx = Derivation(ring, {x: x**2})
        stream = io.StringIO()
        logger = LibraryLogger(stream=stream)

        checker = DiffDegreeLieGenerationChecker(
            [dx, x2dx],
            2,
            verbose=True,
            logger=logger,
            log_every_s=0,
            max_iter=200,
        )
        checker.run()

        output = stream.getvalue()
        self.assertIn("F0", output)
        self.assertIn("F1", output)


class DiffDegreeCliTests(unittest.TestCase):
    def _run_cli(self, *args: str):
        env = os.environ.copy()
        env["HOME"] = "/tmp"
        return subprocess.run(
            [sys.executable, "diff_deg_check.py", *args],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )

    def test_dims_mode(self):
        proc = self._run_cli("--mode", "dims")

        self.assertEqual(proc.returncode, 0)
        self.assertIn("dim W_2[coeff<=5] = 42", proc.stdout)

    def test_beldiev_mode_succeeds_for_n2(self):
        proc = self._run_cli(
            "--mode",
            "beldiev",
            "--n",
            "2",
            "--d",
            "2",
            "--Dmax",
            "8",
            "--log-every",
            "0",
        )

        self.assertEqual(proc.returncode, 0)
        self.assertIn("mode = beldiev, n = 2, d = 2", proc.stdout)
        self.assertIn("status = PASS", proc.stdout)

    def test_andrist_mode_succeeds_for_n2(self):
        proc = self._run_cli(
            "--mode",
            "andrist",
            "--n",
            "2",
            "--d",
            "2",
            "--Dmax",
            "8",
            "--log-every",
            "0",
        )

        self.assertEqual(proc.returncode, 0)
        self.assertIn("mode = andrist, n = 2, d = 2", proc.stdout)
        self.assertIn("status = PASS", proc.stdout)

    def test_andrist_mode_succeeds_for_n3(self):
        proc = self._run_cli(
            "--mode",
            "andrist",
            "--n",
            "3",
            "--d",
            "2",
            "--Dmax",
            "12",
            "--log-every",
            "0",
        )

        self.assertEqual(proc.returncode, 0)
        self.assertIn("mode = andrist, n = 3, d = 2", proc.stdout)
        self.assertIn("status = PASS", proc.stdout)

    def test_hypothesis_mode_accepts_fixed_second_generator(self):
        proc = self._run_cli(
            "--mode",
            "hypo",
            "--n",
            "1",
            "--d",
            "2",
            "--g1",
            "x^2*dx",
            "--g2",
            "dx",
            "--log-every",
            "0",
        )

        self.assertEqual(proc.returncode, 0)
        self.assertIn("G2 mode = fixed", proc.stdout)
        self.assertIn("Гипотеза ПОДТВЕРЖДЕНА", proc.stdout)


if __name__ == "__main__":
    unittest.main()
