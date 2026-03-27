import pathlib
import subprocess
import sys
import unittest

from numeric_enumerate import build_basis_terms, estimate_candidate_count, parse_coefficients

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]


class NumericEnumerateTests(unittest.TestCase):
    def test_build_basis_terms_excludes_constants_by_default(self):
        terms = build_basis_terms(max_x_deg=1, max_y_deg=1, include_constants=False)

        self.assertEqual(len(terms), 6)

    def test_parse_coefficients_deduplicates_and_ignores_zero(self):
        coeffs = parse_coefficients("1, -1, 0, 1")

        self.assertEqual(coeffs, [1.0, -1.0])

    def test_estimate_candidate_count_matches_small_example(self):
        self.assertEqual(estimate_candidate_count(6, 2, 1, 2), 72)

    def test_count_only_cli_reports_search_space(self):
        proc = subprocess.run(
            [
                sys.executable,
                "numeric_enumerate.py",
                "--max-x-deg",
                "1",
                "--max-y-deg",
                "1",
                "--coeffs",
                "1,-1",
                "--max-terms",
                "2",
                "--count-only",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(proc.returncode, 0)
        self.assertIn("basis terms      = 6", proc.stdout)
        self.assertIn("candidate count  = 72", proc.stdout)


if __name__ == "__main__":
    unittest.main()
