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


from lib.numeric.bracket import bracket_vectors, bracket_full_elements


class BracketVectorsTests(unittest.TestCase):
    """Unit tests for low-level bracket computation."""

    def test_partial_x_bracket_x_partial_y_gives_partial_y(self):
        """[d/dx, x*d/dy] = d/dy in W_2."""
        u = np.array([1.0, 0.0])
        v = np.array([0.0, 0.0, 1.0, 0.0])
        result = bracket_vectors(u, -1, v, 0, 2)
        expected = np.array([0.0, 1.0])
        np.testing.assert_allclose(result, expected)

    def test_commuting_diagonal_fields_give_zero(self):
        """[x*d/dx, y*d/dy] = 0 in W_2."""
        u = np.array([1.0, 0.0, 0.0, 0.0])
        v = np.array([0.0, 0.0, 0.0, 1.0])
        result = bracket_vectors(u, 0, v, 0, 2)
        np.testing.assert_allclose(result, np.zeros(4), atol=1e-14)

    def test_antisymmetry(self):
        """[u, v] = -[v, u] for bracket_vectors."""
        u = np.array([1.0, 0.0, 0.0, 0.0])
        v = np.array([0.0, 1.0, 0.0, 0.0])
        from lib.numeric.bracket import precompute_structure_constants
        sc_pq = precompute_structure_constants(2, 0, 0)
        fwd = bracket_vectors(u, 0, v, 0, 2, sc=sc_pq)
        rev = bracket_vectors(v, 0, u, 0, 2, sc=sc_pq)
        np.testing.assert_allclose(fwd, -rev, atol=1e-14)


class BracketFullElementsTests(unittest.TestCase):
    """Unit tests for full-element bracket."""

    def test_self_bracket_is_zero(self):
        """[G, G] must be zero for any element."""
        from lib.numeric.structure import StructureConstantCache
        sc = StructureConstantCache(n=2, D_max=3)
        elem = {
            -1: np.array([1.0, 0.0]),
            0: np.array([0.0, 0.0, 1.0, 0.0]),
        }
        result = bracket_full_elements(elem, elem, 2, sc)
        for k, vec in result.items():
            np.testing.assert_allclose(
                vec, np.zeros_like(vec), atol=1e-14,
                err_msg=f"self-bracket nonzero at degree {k}",
            )

    def test_bracket_matches_bracket_vectors_for_homogeneous(self):
        """For single-degree elements, bracket_full_elements matches bracket_vectors."""
        from lib.numeric.structure import StructureConstantCache
        sc = StructureConstantCache(n=2, D_max=3)
        f = {-1: np.array([1.0, 0.0])}
        g = {0: np.array([0.0, 0.0, 1.0, 0.0])}
        result = bracket_full_elements(f, g, 2, sc)
        expected = bracket_vectors(
            np.array([1.0, 0.0]), -1,
            np.array([0.0, 0.0, 1.0, 0.0]), 0,
            2,
        )
        self.assertIn(-1, result)
        np.testing.assert_allclose(result[-1], expected)

    def test_cross_bracket_of_non_homogeneous_elements(self):
        """[G1, G2] where both span multiple degrees."""
        from lib.numeric.structure import StructureConstantCache
        sc = StructureConstantCache(n=2, D_max=5)
        from lib.numeric.components import decompose_generator
        g1_spec = {0: {(1, 2): 2.0, (3, 1): -1.0}}
        g2_spec = {0: {(0, 0): 1.0}, 1: {(1, 0): 1.0}}
        elem1 = decompose_generator(g1_spec, 2, 5)
        elem2 = decompose_generator(g2_spec, 2, 5)
        result = bracket_full_elements(elem1, elem2, 2, sc)
        has_nonzero = any(np.linalg.norm(v) > 1e-14 for v in result.values())
        self.assertTrue(has_nonzero, "bracket of G1, G2 should be nonzero")


if __name__ == "__main__":
    unittest.main()
