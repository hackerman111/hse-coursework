import io
import pathlib
import subprocess
import sys
import unittest

import numpy as np

from lib.numeric import (
    StructureConstantCache,
    check_full_lie_generation,
    parse_generator_spec,
)
from lib.numeric.criteria import evaluate_full_lie_generation
from lib.numeric.indexing import dim_Wn_leq_d
from lib.numeric.recipes import make_beldiev_generators
from lib.numeric.solver import NumericLieGenerationChecker
from lib.numeric.subspace import OrthonormalBasis
from lib.utils import LibraryLogger

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]


class OrthonormalBasisTests(unittest.TestCase):
    def test_try_add_many_reports_rank_growth(self):
        basis = OrthonormalBasis(3)

        self.assertEqual(basis.try_add_many(np.array([[1.0, 0.0, 0.0]])), 1)
        self.assertEqual(basis.try_add_many(np.array([[1.0, 0.0, 0.0]])), 0)
        self.assertEqual(
            basis.try_add_many(np.array([[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])),
            1,
        )
        self.assertEqual(basis.rank, 2)


class LibraryLoggerTests(unittest.TestCase):
    def test_logger_writes_banner_rule_and_key_value_to_custom_stream(self):
        stream = io.StringIO()
        logger = LibraryLogger(stream=stream)

        logger.banner("Numeric Check", ["n = 2, d = 3"])
        logger.rule("Trial 1/2")
        logger.kv("G1", "dx")

        output = stream.getvalue()
        self.assertIn("Numeric Check", output)
        self.assertIn("Trial 1/2", output)
        self.assertIn("G1 = dx", output)


class NumericSolverTests(unittest.TestCase):
    def test_beldiev_generators_fill_w2_up_to_degree_7(self):
        checker = NumericLieGenerationChecker(
            2,
            7,
            make_beldiev_generators(2),
            D_max=7,
            verbose=False,
        )

        result = checker.run()

        self.assertTrue(result.success)
        self.assertGreater(result.iterations, 1)
        self.assertTrue(all(ds.is_full for ds in result.degrees))

    def test_hypothesis_candidate_y3x3_is_not_counterexample_up_to_degree_2(self):
        first_gen = parse_generator_spec("y^3*dx + x^3*dy", 2)
        second_gen = parse_generator_spec("dx + x*dy", 2)
        checker = NumericLieGenerationChecker(
            2,
            2,
            [first_gen, second_gen],
            D_max=3,
            verbose=False,
        )

        result = checker.run()

        self.assertTrue(result.success)
        self.assertTrue(all(ds.is_full for ds in result.degrees))

    def test_hypothesis_candidate_x3y3_diagonal_is_not_counterexample_up_to_degree_2(self):
        first_gen = parse_generator_spec("x^3*dx + y^3*dy", 2)
        second_gen = parse_generator_spec("dx + x*dy", 2)
        checker = NumericLieGenerationChecker(
            2,
            2,
            [first_gen, second_gen],
            D_max=3,
            verbose=False,
        )

        result = checker.run()

        self.assertTrue(result.success)
        self.assertTrue(all(ds.is_full for ds in result.degrees))

    def test_hypothesis_candidate_radial_quartic_dx_is_not_counterexample_up_to_degree_2(self):
        first_gen = parse_generator_spec("x^4*dx + 2*x^2*y^2*dx + y^4*dx", 2)
        second_gen = parse_generator_spec("dx + x*dy", 2)
        checker = NumericLieGenerationChecker(
            2,
            2,
            [first_gen, second_gen],
            D_max=4,
            verbose=False,
        )

        result = checker.run()

        self.assertTrue(result.success)
        self.assertTrue(all(ds.is_full for ds in result.degrees))

    def test_structure_constant_cache_can_be_reused_explicitly(self):
        first_gen = parse_generator_spec("x^3*dx + y^3*dy", 2)
        second_gen = parse_generator_spec("dx + x*dy", 2)
        structure_constants = StructureConstantCache(n=2, D_max=3)

        first = NumericLieGenerationChecker(
            2,
            2,
            [first_gen, second_gen],
            D_max=3,
            verbose=False,
            structure_constants=structure_constants,
        ).run()
        second = NumericLieGenerationChecker(
            2,
            2,
            [first_gen, second_gen],
            D_max=3,
            verbose=False,
            structure_constants=structure_constants,
        ).run()

        self.assertTrue(first.success)
        self.assertEqual(first.degrees, second.degrees)

    def test_solver_uses_injected_library_logger(self):
        stream = io.StringIO()
        logger = LibraryLogger(stream=stream)
        checker = NumericLieGenerationChecker(
            2,
            1,
            make_beldiev_generators(2),
            D_max=1,
            verbose=True,
            logger=logger,
        )

        checker.run()

        output = stream.getvalue()
        self.assertIn("iter   0:", output)
        self.assertIn("iter   1:", output)

    def test_solver_respects_verbose_flag_with_injected_logger(self):
        stream = io.StringIO()
        logger = LibraryLogger(stream=stream)
        checker = NumericLieGenerationChecker(
            2,
            1,
            make_beldiev_generators(2),
            D_max=1,
            verbose=False,
            logger=logger,
        )

        checker.run()

        self.assertEqual(stream.getvalue(), "")

    def test_numeric_md_dimension_examples(self):
        expected = [
            (2, 5, 56),
            (2, 10, 156),
            (3, 5, 252),
            (3, 10, 1092),
            (5, 5, 2310),
            (5, 10, 21840),
        ]

        for n, d, total_dim in expected:
            with self.subTest(n=n, d=d):
                self.assertEqual(dim_Wn_leq_d(n, d), total_dim)

    def test_detailed_progress_summary_lists_all_degrees(self):
        checker = NumericLieGenerationChecker(
            2,
            2,
            make_beldiev_generators(2),
            D_max=3,
            verbose=False,
            log_every_s=1.0,
        )

        summary = checker._build_detailed_status(
            iteration=2,
            elapsed_s=1.25,
            active_degree=1,
            active_pair=(0, 1),
        )

        self.assertIn("=== Progress Summary ===", summary)
        self.assertIn("elapsed = 1.250s, iteration = 2", summary)
        self.assertIn("active degree = 1", summary)
        self.assertIn("active pair = (0, 1)", summary)
        self.assertIn("Degree  -1 [target ]", summary)
        self.assertIn("Degree   3 [closure]", summary)

    def test_degree_two_criterion_accepts_full_result(self):
        first_gen = parse_generator_spec("x^3*dx + y^3*dy", 2)
        second_gen = parse_generator_spec("dx + x*dy", 2)

        check = check_full_lie_generation(
            n=2,
            generators=[first_gen, second_gen],
            D_max=3,
        )

        self.assertTrue(check.success)
        self.assertTrue(check.criterion_result.satisfied)
        self.assertTrue(check.numeric_result.success)

    def test_degree_two_criterion_rejects_result_that_does_not_reach_degree_2(self):
        checker = NumericLieGenerationChecker(
            2,
            1,
            make_beldiev_generators(2),
            D_max=1,
            verbose=False,
        )

        criterion = evaluate_full_lie_generation(checker.run())

        self.assertFalse(criterion.satisfied)
        self.assertIn("only checks up to degree 1", criterion.reason)


class NumericCliTests(unittest.TestCase):
    def test_parse_generator_spec_supports_sum_and_aliases(self):
        gen = parse_generator_spec("dx + x^3*dx - 2*y*dy", 2)

        self.assertEqual(gen[0][(0, 0)], 1.0)
        self.assertEqual(gen[0][(3, 0)], 1.0)
        self.assertEqual(gen[1][(0, 1)], -2.0)

    def test_help_mentions_new_generator_and_logging_options(self):
        proc = subprocess.run(
            [sys.executable, "numeric_check.py", "--help"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(proc.returncode, 0)
        self.assertIn("--g1", proc.stdout)
        self.assertIn("--g2", proc.stdout)
        self.assertIn("--log-every", proc.stdout)
        self.assertIn("criterion", proc.stdout)
        self.assertIn("generator-spec", proc.stdout)

    def test_hypothesis_mode_accepts_fixed_second_generator(self):
        proc = subprocess.run(
            [
                sys.executable,
                "numeric_check.py",
                "--mode",
                "hypo",
                "--n",
                "1",
                "--d",
                "5",
                "--g1",
                "x^2*dx",
                "--g2",
                "dx + x^3*dx",
                "--log-every",
                "0",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(proc.returncode, 0)
        self.assertIn("G2 mode = fixed", proc.stdout)
        self.assertIn("Trial 1/1", proc.stdout)
        self.assertIn("Гипотеза ПОДТВЕРЖДЕНА", proc.stdout)

    def test_criterion_mode_accepts_fixed_generator_pair(self):
        proc = subprocess.run(
            [
                sys.executable,
                "numeric_check.py",
                "--mode",
                "criterion",
                "--n",
                "2",
                "--g1",
                "x^3*dx + y^3*dy",
                "--g2",
                "dx + x*dy",
                "--Dmax",
                "3",
                "--log-every",
                "0",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(proc.returncode, 0)
        self.assertIn("Проверка критерия полной порождаемости", proc.stdout)
        self.assertIn("PASS: full-lie-from-degree-bound", proc.stdout)


if __name__ == "__main__":
    unittest.main()
