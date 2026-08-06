from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "codes"))

from prepare_latest_bcc_qe_input import (  # noqa: E402
    generate_paramagnetic_spins_magnetic_sqs,
    magnetic_neighbor_shells,
    magnetic_shell_correlations,
    magnetic_sqs_objective,
)


def ideal_bcc_4x4x4() -> tuple[np.ndarray, np.ndarray]:
    positions = []
    for ix in range(4):
        for iy in range(4):
            for iz in range(4):
                positions.extend(
                    [
                        [ix / 4, iy / 4, iz / 4],
                        [(ix + 0.5) / 4, (iy + 0.5) / 4, (iz + 0.5) / 4],
                    ]
                )
    return np.asarray(positions), np.eye(3) * 10.2


class MagneticSqsTests(unittest.TestCase):
    def test_sqs_is_reproducible_zero_net_and_decorrelated(self) -> None:
        positions, cell = ideal_bcc_4x4x4()
        first = generate_paramagnetic_spins_magnetic_sqs(
            positions, cell, m_abs=0.35, seed=12345, optimization_steps=4000
        )
        second = generate_paramagnetic_spins_magnetic_sqs(
            positions, cell, m_abs=0.35, seed=12345, optimization_steps=4000
        )
        spins, theta, phi, net = first
        np.testing.assert_array_equal(spins, second[0])
        np.testing.assert_allclose(np.linalg.norm(spins, axis=1), 0.35, atol=1e-14)
        np.testing.assert_allclose(net, 0.0, atol=1e-13)
        self.assertTrue(np.all((theta >= 0.0) & (theta <= 180.0)))
        self.assertTrue(np.all((phi >= 0.0) & (phi < 360.0)))
        shells = magnetic_neighbor_shells(positions, cell, nshells=5)
        correlations = magnetic_shell_correlations(spins / 0.35, shells)
        self.assertLess(magnetic_sqs_objective(spins / 0.35, shells), 2.0e-4)
        self.assertLess(np.max(np.abs(correlations)), 0.03)


if __name__ == "__main__":
    unittest.main()
