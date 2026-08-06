from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "codes"))

import numpy as np

from data_compress import parse_forces_block, parse_qe_aimd_output  # noqa: E402


class QeForceParserTests(unittest.TestCase):
    def test_accepts_qe_overflowed_atom_type(self) -> None:
        lines = ["Forces acting on atoms (cartesian axes, Ry/au):\n"]
        lines.extend(f"atom {index} type {index} force = 0 0 0\n" for index in range(1, 99))
        lines.extend(
            [
                "atom   99 type 99   force =  0.1 0.2 0.3\n",
                "atom   99 type 99   force =  0.1 0.2 0.3\n",
                "atom  100 type **   force = -0.4 0.5 0.6\n",
            ]
        )
        forces, _ = parse_forces_block(lines, 0, 100)
        self.assertEqual(forces[98:].tolist(), [[0.1, 0.2, 0.3], [-0.4, 0.5, 0.6]])

    def test_accepts_d_exponents_and_arbitrary_type_token(self) -> None:
        lines = [
            "Forces acting on atoms (cartesian axes, Ry/au):\n",
            "atom 1 type *** force = 1.0D-02 -2.5d+00 3.0E-01\n",
        ]
        forces, _ = parse_forces_block(lines, 0, 1)
        np.testing.assert_allclose(forces, [[0.01, -2.5, 0.3]])

    def test_reports_missing_force_indices(self) -> None:
        lines = [
            "Forces acting on atoms (cartesian axes, Ry/au):\n",
            "atom 1 type 1 force = 0 0 0\n",
            "Total force = 0\n",
        ]
        with self.assertRaisesRegex(ValueError, r"missing \[2\]"):
            parse_forces_block(lines, 0, 2)

    def test_fixed_cell_md_frame_is_valid_without_cell_parameters_block(self) -> None:
        output = """number of atoms/cell = 2
lattice parameter (alat) = 5.0 a.u.
crystal axes: (cart. coord. in units of alat)
 a(1) = ( 1.0 0.0 0.0 )
 a(2) = ( 0.0 1.0 0.0 )
 a(3) = ( 0.0 0.0 1.0 )
Cartesian axes
site n.     atom                  positions (alat units)
 1 Fe tau( 1) = ( 0.0 0.0 0.0 )
 2 Fe tau( 2) = ( 0.5 0.5 0.5 )
Entering Dynamics: iteration = 1
 time = 0.001 pico-seconds
ATOMIC_POSITIONS (crystal)
Fe 0.01 0.02 0.03
Fe 0.51 0.52 0.53
Forces acting on atoms (cartesian axes, Ry/au):

atom 1 type 1 force = 1.0D-2 2.0D-2 3.0D-2
atom 2 type ** force = -1.0D-2 -2.0D-2 -3.0D-2
Total force = 0.0
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fe.out"
            path.write_text(output)
            data = parse_qe_aimd_output(path)

        self.assertEqual(data["nsteps"], 1)
        self.assertTrue(data["frame_valid"].tolist() == [True])
        self.assertEqual(data["cell_parameters_unit"].tolist(), ["alat"])
        np.testing.assert_allclose(data["cell_parameters"][0], np.eye(3))

    def test_rejects_archive_with_no_complete_force_frames(self) -> None:
        output = """number of atoms/cell = 1
lattice parameter (alat) = 5.0 a.u.
crystal axes: (cart. coord. in units of alat)
 a(1) = ( 1.0 0.0 0.0 )
 a(2) = ( 0.0 1.0 0.0 )
 a(3) = ( 0.0 0.0 1.0 )
Cartesian axes
site n.     atom                  positions (alat units)
 1 Fe tau( 1) = ( 0.0 0.0 0.0 )
Entering Dynamics: iteration = 1
ATOMIC_POSITIONS (crystal)
Fe 0.0 0.0 0.0
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fe.out"
            path.write_text(output)
            with self.assertRaisesRegex(ValueError, "No complete force frames"):
                parse_qe_aimd_output(path)

    def test_preserves_noncollinear_total_magnetization_vector(self) -> None:
        output = """number of atoms/cell = 1
lattice parameter (alat) = 5.0 a.u.
crystal axes: (cart. coord. in units of alat)
 a(1) = ( 1.0 0.0 0.0 )
 a(2) = ( 0.0 1.0 0.0 )
 a(3) = ( 0.0 0.0 1.0 )
Cartesian axes
site n.     atom                  positions (alat units)
 1 Fe tau( 1) = ( 0.0 0.0 0.0 )
Entering Dynamics: iteration = 1
ATOMIC_POSITIONS (crystal)
Fe 0.0 0.0 0.0
total magnetization = -1.25 2.50 -3.75 Bohr mag/cell
absolute magnetization = 4.50 Bohr mag/cell
Forces acting on atoms (cartesian axes, Ry/au):
atom 1 type 1 force = 0.0 0.0 0.0
Total force = 0.0
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fe.out"
            path.write_text(output)
            data = parse_qe_aimd_output(path)

        np.testing.assert_allclose(
            data["mag_total_vector_Bohr"], [[-1.25, 2.50, -3.75]]
        )
        self.assertEqual(data["mag_total_Bohr"].tolist(), [-1.25])

    def test_preserves_atom_resolved_noncollinear_magnetization(self) -> None:
        output = """number of atoms/cell = 2
lattice parameter (alat) = 5.0 a.u.
crystal axes: (cart. coord. in units of alat)
 a(1) = ( 1.0 0.0 0.0 )
 a(2) = ( 0.0 1.0 0.0 )
 a(3) = ( 0.0 0.0 1.0 )
Cartesian axes
site n.     atom                  positions (alat units)
 1 Fe tau( 1) = ( 0.0 0.0 0.0 )
 2 Fe tau( 2) = ( 0.5 0.5 0.5 )
Entering Dynamics: iteration = 1
ATOMIC_POSITIONS (crystal)
Fe 0.0 0.0 0.0
Fe 0.5 0.5 0.5
atom number 1 relative position : 0.0 0.0 0.0
charge : 13.0
magnetization : 0.10 -0.20 0.30
atom number 2 relative position : 0.5 0.5 0.5
charge : 13.0
magnetization : -0.40 0.50 -0.60
Forces acting on atoms (cartesian axes, Ry/au):
atom 1 type 1 force = 0.0 0.0 0.0
atom 2 type 2 force = 0.0 0.0 0.0
Total force = 0.0
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fe.out"
            path.write_text(output)
            data = parse_qe_aimd_output(path)

        np.testing.assert_allclose(
            data["local_magnetization_Bohr"],
            [[[0.10, -0.20, 0.30], [-0.40, 0.50, -0.60]]],
        )
        self.assertEqual(data["local_magnetization_frame_valid"].tolist(), [True])


if __name__ == "__main__":
    unittest.main()
