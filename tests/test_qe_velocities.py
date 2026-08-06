from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "codes"))

from generate_qe_maxwell_velocities import (  # noqa: E402
    finalize_velocities_for_qe,
    parse_atomic_velocities,
    replace_or_append_atomic_velocities,
    sample_maxwell_velocities_au,
    temperature_from_velocities_au,
)


class QeVelocityTests(unittest.TestCase):
    def test_128_fe_atoms_reach_requested_temperatures(self) -> None:
        masses = np.full(128, 55.845)
        labels = ["Fe"] * 128
        for target in (2000.0, 4000.0, 6000.0):
            with self.subTest(target=target):
                sampled = sample_maxwell_velocities_au(masses, target, seed=4000)
                velocities, _, measured = finalize_velocities_for_qe(labels, sampled, masses, target)
                self.assertAlmostEqual(measured, target, places=8)
                np.testing.assert_allclose(np.sum(masses[:, None] * velocities, axis=0), 0.0, atol=1e-12)

    def test_doubling_velocities_quadruples_temperature(self) -> None:
        masses = np.full(128, 55.845)
        velocities = sample_maxwell_velocities_au(masses, 4000.0, seed=7)
        initial = temperature_from_velocities_au(velocities, masses, remove_com=False)
        doubled = temperature_from_velocities_au(2.0 * velocities, masses, remove_com=False)
        self.assertAlmostEqual(doubled, 4.0 * initial, places=10)

    def test_formatted_qe_input_round_trip_retains_temperature(self) -> None:
        masses = np.full(128, 55.845)
        labels = ["Fe"] * 128
        sampled = sample_maxwell_velocities_au(masses, 4000.0, seed=4000)
        _, block, _ = finalize_velocities_for_qe(labels, sampled, masses, 4000.0)
        base = "ATOMIC_SPECIES\nFe 55.845 Fe.UPF\n\nATOMIC_POSITIONS angstrom\n" + "\n".join(
            f"Fe {index} 0 0" for index in range(128)
        ) + "\n\nK_POINTS gamma\n"
        generated = replace_or_append_atomic_velocities(base, block)
        parsed_labels, parsed = parse_atomic_velocities(generated.splitlines())
        self.assertEqual(parsed_labels, labels)
        self.assertAlmostEqual(temperature_from_velocities_au(parsed, masses, remove_com=True), 4000.0, places=8)

    def test_random_seed_reproducibility(self) -> None:
        masses = np.full(128, 55.845)
        first = sample_maxwell_velocities_au(masses, 4000.0, seed=123)
        second = sample_maxwell_velocities_au(masses, 4000.0, seed=123)
        different = sample_maxwell_velocities_au(masses, 4000.0, seed=124)
        self.assertTrue(np.array_equal(first, second))
        self.assertFalse(np.array_equal(first, different))

    def test_mixed_masses_use_mass_weighted_com_removal(self) -> None:
        masses = np.array([1.0, 4.0, 16.0, 55.845])
        labels = ["H", "He", "O", "Fe"]
        sampled = sample_maxwell_velocities_au(masses, 3000.0, seed=9)
        velocities, _, measured = finalize_velocities_for_qe(labels, sampled, masses, 3000.0)
        np.testing.assert_allclose(np.sum(masses[:, None] * velocities, axis=0), 0.0, atol=1e-12)
        self.assertAlmostEqual(measured, 3000.0, places=8)

if __name__ == "__main__":
    unittest.main()
