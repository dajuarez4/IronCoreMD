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
    supports_lattice_plots: bool = True

    def dataset_dir(self, repo_root: Path) -> Path:
        return repo_root / "dataset" / self.dataset_subdir

    def lattice_parameter_from_volume_per_atom(self, volume_per_atom: float | np.ndarray) -> float | np.ndarray:
        if self.key == "hcp":
            raise ValueError("HCP does not have a single lattice-parameter mapping from volume alone.")
        return np.cbrt(self.conventional_atoms * np.asarray(volume_per_atom, dtype=float))

    def conventional_cell_volume(self, volume_per_atom: float) -> float:
        return float(self.conventional_atoms) * float(volume_per_atom)

    def infer_supercell(self, cell_ang: np.ndarray, natoms: int) -> tuple[int, int, int]:
        if self.key == "hcp":
            if natoms % self.conventional_atoms != 0:
                raise ValueError(
                    f"Could not infer an HCP supercell for natoms={natoms}; expected a multiple of {self.conventional_atoms}."
                )
            ncell = natoms // self.conventional_atoms
            n = round(ncell ** (1.0 / 3.0))
            if n**3 == ncell:
                return (n, n, n)
            lengths = np.linalg.norm(np.asarray(cell_ang, dtype=float), axis=1)
            best: tuple[float, tuple[int, int, int]] | None = None
            for nx in range(1, ncell + 1):
                if ncell % nx != 0:
                    continue
                ny_nz = ncell // nx
                for ny in range(1, ny_nz + 1):
                    if ny_nz % ny != 0:
                        continue
                    nz = ny_nz // ny
                    if nx != ny:
                        continue
                    a_from_x = lengths[0] / nx
                    a_from_y = lengths[1] / ny
                    mismatch = abs(a_from_x - a_from_y) / max(1.0e-12, 0.5 * (a_from_x + a_from_y))
                    shape_penalty = abs(nx - nz)
                    candidate = (mismatch, float(shape_penalty), (int(nx), int(ny), int(nz)))
                    if best is None or candidate < best:
                        best = candidate
            if best is None or best[0] > 1.0e-3:
                raise ValueError(f"Could not infer a valid HCP supercell for natoms={natoms}.")
            return best[2]

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

    def primitive_cell(
        self,
        cell_ang: np.ndarray,
        natoms: int,
        supercell: tuple[int, int, int] | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        if self.key == "hcp":
            nx, _, nz = supercell if supercell is not None else self.infer_supercell(cell_ang, natoms)
            cell = np.asarray(cell_ang, dtype=float)
            lattice_a = float(np.linalg.norm(cell[0])) / float(nx)
            lattice_c = float(np.linalg.norm(cell[2])) / float(nz)
            primitive = np.array(
                [
                    [lattice_a, 0.0, 0.0],
                    [0.5 * lattice_a, 0.5 * np.sqrt(3.0) * lattice_a, 0.0],
                    [0.0, 0.0, lattice_c],
                ],
                dtype=float,
            )
            basis = np.array([[0.0, 0.0, 0.0], [2.0 / 3.0, 1.0 / 3.0, 0.5]], dtype=float)
            return primitive, basis

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
    "hcp": PhaseSpec(
        key="hcp",
        title="HCP",
        long_name="hexagonal close packed",
        dataset_subdir="hcp",
        conventional_atoms=2,
        default_supercell=(4, 4, 4),
        qpoint_bravais="HEX",
        qpoint_segments=(("GM", "M"), ("M", "K"), ("K", "GM"), ("GM", "A"), ("A", "L")),
        supports_lattice_plots=False,
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
