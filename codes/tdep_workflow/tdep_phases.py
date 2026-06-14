#!/usr/bin/env python3
"""Phase-specific settings for the TDEP workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class PhaseSpec:
    key: str
    title: str
    long_name: str
    dataset_subdir: str
    conventional_atoms: int
    default_supercell: tuple[int, int, int]
    qpoint_bravais: str
    qpoint_segments: tuple[tuple[str, str], ...]

    def dataset_dir(self, repo_root: Path) -> Path:
        return repo_root / "dataset" / self.dataset_subdir

    def lattice_parameter_from_volume_per_atom(self, volume_per_atom: float | np.ndarray) -> float | np.ndarray:
        return np.cbrt(self.conventional_atoms * np.asarray(volume_per_atom, dtype=float))

    def conventional_cell_volume(self, volume_per_atom: float) -> float:
        return float(self.conventional_atoms) * float(volume_per_atom)

    def infer_supercell(self, cell_ang: np.ndarray, natoms: int) -> tuple[int, int, int]:
        cell = np.asarray(cell_ang, dtype=float)
        volume_per_atom = abs(float(np.linalg.det(cell))) / float(natoms)
        lattice_a = float(self.lattice_parameter_from_volume_per_atom(volume_per_atom))
        lengths = np.linalg.norm(cell, axis=1)
        counts = np.rint(lengths / lattice_a).astype(int)
        if np.any(counts < 1):
            raise ValueError(f"Could not infer a valid supercell for phase={self.key!r}: {counts.tolist()}")
        if not np.allclose(lengths, counts * lattice_a, rtol=1.0e-3, atol=1.0e-3):
            raise ValueError(
                f"Cell lengths {lengths.tolist()} are not commensurate with the inferred {self.key} lattice "
                f"parameter {lattice_a:.6f} A"
            )
        return tuple(int(value) for value in counts)

    def primitive_cell(self, cell_ang: np.ndarray, natoms: int) -> tuple[np.ndarray, np.ndarray]:
        volume_per_atom = abs(float(np.linalg.det(np.asarray(cell_ang, dtype=float)))) / float(natoms)
        lattice_a = float(self.lattice_parameter_from_volume_per_atom(volume_per_atom))
        if self.key == "bcc":
            cell = np.array(
                [
                    [-0.5, 0.5, 0.5],
                    [0.5, -0.5, 0.5],
                    [0.5, 0.5, -0.5],
                ],
                dtype=float,
            )
        elif self.key == "fcc":
            cell = np.array(
                [
                    [0.0, 0.5, 0.5],
                    [0.5, 0.0, 0.5],
                    [0.5, 0.5, 0.0],
                ],
                dtype=float,
            )
        else:
            raise ValueError(f"Unsupported phase: {self.key}")
        return lattice_a * cell, np.array([[0.0, 0.0, 0.0]], dtype=float)


PHASE_SPECS = {
    "bcc": PhaseSpec(
        key="bcc",
        title="BCC",
        long_name="body centered cubic",
        dataset_subdir="bcc",
        conventional_atoms=2,
        default_supercell=(4, 4, 4),
        qpoint_bravais="BCC",
        qpoint_segments=(("GM", "H"), ("H", "N"), ("N", "GM"), ("GM", "P"), ("P", "H")),
    ),
    "fcc": PhaseSpec(
        key="fcc",
        title="FCC",
        long_name="face centered cubic",
        dataset_subdir="fcc",
        conventional_atoms=4,
        default_supercell=(4, 4, 2),
        qpoint_bravais="FCC",
        qpoint_segments=(("GM", "X"), ("X", "U"), ("K", "GM"), ("GM", "L")),
    ),
}


def normalize_phase(phase: str) -> str:
    key = str(phase).strip().lower()
    if key not in PHASE_SPECS:
        choices = ", ".join(sorted(PHASE_SPECS))
        raise ValueError(f"Unsupported phase {phase!r}. Expected one of: {choices}")
    return key


def get_phase_spec(phase: str) -> PhaseSpec:
    return PHASE_SPECS[normalize_phase(phase)]


def write_qpoints_dispersion(path: Path, phase: str, points_per_path: int = 120) -> None:
    spec = get_phase_spec(phase)
    lines = [
        spec.qpoint_bravais,
        f"{points_per_path:5d}",
        f"{len(spec.qpoint_segments):5d}",
        *(f"{start:<3} {end}" for start, end in spec.qpoint_segments),
    ]
    path.write_text("\n".join(lines) + "\n")
