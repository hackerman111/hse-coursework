"""Тест early-stop callback в LieBasisSolver."""

import unittest

from sage.all import PolynomialRing, QQ

from lib.derivation import LieDerivation
from lib.solver.basis import LieBasisSolver


class TestStopCondition(unittest.TestCase):
    def setUp(self):
        self.R = PolynomialRing(QQ, "x,y", order="degrevlex")
        self.x, self.y = self.R.gens()

    def _make(self, mapping):
        images = [self.R.zero(), self.R.zero()]
        for gen, image in mapping.items():
            idx = list(self.R.gens()).index(gen)
            images[idx] = self.R(image)
        return LieDerivation.from_mapping(self.R, images)

    def test_stop_condition_fires_before_full_completion(self):
        """Callback, который останавливает после первого найденного элемента базиса."""
        dx = self._make({self.x: 1})
        x_dy = self._make({self.y: self.x})

        call_count = [0]

        def stop_after_two(solver):
            call_count[0] += 1
            return len(solver.basis) >= 2

        solver = LieBasisSolver([dx, x_dy], max_iter=1000, stop_condition=stop_after_two)
        solver.run()

        self.assertGreaterEqual(call_count[0], 1)
        self.assertGreaterEqual(len(solver.basis), 2)
        self.assertLessEqual(solver.iter_count, solver.max_iter)

    def test_no_stop_condition_runs_to_completion(self):
        """Без callback солвер работает как раньше."""
        dx = self._make({self.x: 1})
        x_dy = self._make({self.y: self.x})

        solver = LieBasisSolver([dx, x_dy], max_iter=1000)
        success = solver.run()

        self.assertTrue(success)
        self.assertTrue(all(solver.targets_found.values()))


if __name__ == "__main__":
    unittest.main()
