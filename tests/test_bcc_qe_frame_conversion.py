from __future__ import annotations

from pathlib import Path
import json
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "codes"))

from prepare_latest_bcc_qe_input import BOHR_TO_ANG, frame_cell_angstrom  # noqa: E402


class ArrayArchive:
    def __init__(self, **arrays) -> None:
        self._arrays = arrays
        self.files = list(arrays)

    def __getitem__(self, key: str):
        return self._arrays[key]


class BccQeFrameConversionTests(unittest.TestCase):
    def test_frame_cell_in_alat_uses_archive_metadata(self) -> None:
        alat_bohr = 14.4564
        archive = ArrayArchive(
            cell_parameters=np.eye(3, dtype=float)[None, :, :],
            cell_parameters_unit=np.asarray(["alat"]),
            metadata_json=np.asarray(json.dumps({"alat_bohr": alat_bohr})),
        )

        converted = frame_cell_angstrom(archive, 0, np.eye(3) * 99.0)

        np.testing.assert_allclose(converted, np.eye(3) * alat_bohr * BOHR_TO_ANG)

    def test_frame_cell_in_alat_requires_alat_metadata(self) -> None:
        archive = ArrayArchive(
            cell_parameters=np.eye(3, dtype=float)[None, :, :],
            cell_parameters_unit=np.asarray(["alat"]),
            metadata_json=np.asarray("{}"),
        )

        with self.assertRaisesRegex(ValueError, "alat_bohr"):
            frame_cell_angstrom(archive, 0, np.eye(3))


if __name__ == "__main__":
    unittest.main()
