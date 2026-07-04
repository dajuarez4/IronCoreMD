#!/usr/bin/env python3
"""Prepare a QE MD input from the last frame of a BCC dataset NPZ.

The workflow is intentionally narrow:
1. Pick the newest BCC NPZ, or an explicit override.
2. Extract the requested frame, using the fixed input cell when per-frame cells
   are absent or invalid.
3. Generate zero-net noncollinear paramagnetic spins, using either paired-random
   or quasi-random directions on the sphere.
4. Generate Maxwell-Boltzmann velocities in QE atomic units.
5. Write sidecar text files plus a final QE input with spins and velocities.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import re
from pathlib import Path
import string

import numpy as np

from generate_qe_maxwell_velocities import (
    format_atomic_velocities_card,
    read_qe_species_sequence,
    remove_center_of_mass_drift_au,
    replace_or_append_atomic_velocities,
    rescale_to_temperature_au,
    sample_maxwell_velocities_au,
    temperature_from_velocities_au,
)

BOHR_TO_ANG = 0.529177210903
TEMPERATURE_RE = re.compile(r"_(\d+)(?:K)?(?:[-_].+)?$")
FAKE_SPECIES_LETTERS = string.ascii_uppercase
MAX_FAKE_SPECIES = len(FAKE_SPECIES_LETTERS) * 99
SPIN_MODE_RANDOM_ZERO_NET = "random_zero_net"
SPIN_MODE_QUASI_RANDOM_ZERO_NET = "quasi_random_zero_net"
QE_FLOAT_RE = r"[+\-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[EeDd][+\-]?\d+)?"
QE_ATOM_SPIN_BLOCK_RE = re.compile(
    rf"atom number\s+(\d+)\s+relative position\s*:\s*"
    rf"({QE_FLOAT_RE})\s+({QE_FLOAT_RE})\s+({QE_FLOAT_RE})\s+"
    rf"charge\s*:\s*({QE_FLOAT_RE})\s+\(integrated on a sphere of radius\s+({QE_FLOAT_RE})\)\s+"
    rf"magnetization\s*:\s*({QE_FLOAT_RE})\s+({QE_FLOAT_RE})\s+({QE_FLOAT_RE})\s+"
    rf"magnetization/charge:\s*({QE_FLOAT_RE})\s+({QE_FLOAT_RE})\s+({QE_FLOAT_RE})\s+"
    rf"polar coord\.: r, theta, phi \[deg\] :\s*({QE_FLOAT_RE})\s+({QE_FLOAT_RE})\s+({QE_FLOAT_RE})",
    re.MULTILINE,
)
SUPPORTED_SPIN_MODES = (
    SPIN_MODE_RANDOM_ZERO_NET,
    SPIN_MODE_QUASI_RANDOM_ZERO_NET,
)


@dataclass
class PreparationResult:
    npz_path: Path
    output_dir: Path
    frame_index: int
    natoms: int
    structure_tag: str
    target_temperature_k: float
    measured_velocity_temperature_k: float
    qe_input_base: Path
    qe_input_final: Path
    velocity_block: Path
    spin_vectors: Path
    spin_parameters: Path
    net_magnetization: tuple[float, float, float]
    spin_mode: str
    spin_source: str
    qe_spin_output_path: Path | None
    mean_starting_magnetization: float
    min_starting_magnetization: float
    max_starting_magnetization: float
    ntypes: int
    first_labels: tuple[str, ...]
    last_labels: tuple[str, ...]
    position_velocity_labels_match: bool

    def as_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["npz_path"] = str(self.npz_path)
        payload["output_dir"] = str(self.output_dir)
        payload["qe_input_base"] = str(self.qe_input_base)
        payload["qe_input_final"] = str(self.qe_input_final)
        payload["velocity_block"] = str(self.velocity_block)
        payload["spin_vectors"] = str(self.spin_vectors)
        payload["spin_parameters"] = str(self.spin_parameters)
        payload["qe_spin_output_path"] = None if self.qe_spin_output_path is None else str(self.qe_spin_output_path)
        payload["net_magnetization"] = [float(value) for value in self.net_magnetization]
        payload["first_labels"] = list(self.first_labels)
        payload["last_labels"] = list(self.last_labels)
        return payload


@dataclass
class QESpinGuess:
    qe_output_path: Path
    atom_numbers: np.ndarray
    charge_in_sphere: np.ndarray
    moment_cart_bohr: np.ndarray
    moment_over_charge_cart: np.ndarray
    local_moment_bohr: np.ndarray
    angle1: np.ndarray
    angle2: np.ndarray
    starting_magnetization: np.ndarray

    def summary_dict(self) -> dict[str, object]:
        return {
            "qe_output_path": str(self.qe_output_path),
            "natoms": int(len(self.atom_numbers)),
            "mean_starting_magnetization": float(np.mean(self.starting_magnetization)),
            "min_starting_magnetization": float(np.min(self.starting_magnetization)),
            "max_starting_magnetization": float(np.max(self.starting_magnetization)),
            "mean_local_moment_bohr": float(np.mean(self.local_moment_bohr)),
            "min_local_moment_bohr": float(np.min(self.local_moment_bohr)),
            "max_local_moment_bohr": float(np.max(self.local_moment_bohr)),
        }


def scalar_string(value) -> str:
    array = np.asarray(value)
    return str(array.item() if array.shape == () else value)


def qe_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def load_metadata(data) -> dict[str, object]:
    if "metadata_json" not in data.files:
        return {}
    value = data["metadata_json"]
    if isinstance(value, np.ndarray):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return json.loads(str(value))


def infer_latest_bcc_npz(dataset_dir: Path) -> Path:
    npz_files = sorted(dataset_dir.glob("*.npz"), key=lambda path: path.stat().st_mtime)
    if not npz_files:
        raise FileNotFoundError(f"No NPZ files found in {dataset_dir}")
    return npz_files[-1]


def parse_temperature_from_name(path: Path) -> float:
    match = TEMPERATURE_RE.search(path.stem)
    if not match:
        raise ValueError(f"Could not infer temperature from file name: {path.name}")
    return float(match.group(1))


def fixed_cell_angstrom(data) -> np.ndarray:
    cell = np.asarray(data["input_cell_parameters"], dtype=float)
    unit = scalar_string(data["input_cell_unit"]).lower()
    if unit in {"angstrom", "ang"}:
        return cell
    if unit in {"bohr", "a.u.", "au"}:
        return cell * BOHR_TO_ANG
    raise ValueError(f"Unsupported input_cell_unit={unit!r}")


def frame_cell_angstrom(data, frame_index: int, fallback_cell: np.ndarray) -> np.ndarray:
    if "cell_parameters" not in data.files:
        return fallback_cell

    frame_cells = np.asarray(data["cell_parameters"], dtype=float)
    if frame_cells.ndim != 3 or frame_index >= len(frame_cells):
        return fallback_cell

    cell = frame_cells[frame_index]
    if not np.isfinite(cell).all():
        return fallback_cell

    units = np.asarray(data["cell_parameters_unit"])
    unit = str(units[frame_index]).lower() if units.ndim > 0 else scalar_string(units).lower()
    if unit in {"none", ""}:
        return fallback_cell
    if unit in {"angstrom", "ang"}:
        return cell
    if unit in {"bohr", "a.u.", "au"}:
        return cell * BOHR_TO_ANG
    raise ValueError(f"Unsupported cell_parameters_unit={unit!r}")


def frame_positions_fractional(data, frame_index: int, cell_ang: np.ndarray) -> np.ndarray:
    positions = np.asarray(data["positions"], dtype=float)[frame_index]
    units = np.asarray(data["positions_unit"])
    unit = str(units[frame_index]).lower() if units.ndim > 0 else scalar_string(units).lower()

    if unit == "crystal":
        return positions
    if unit in {"angstrom", "ang"}:
        return positions @ np.linalg.inv(cell_ang)
    if unit in {"bohr", "a.u.", "au"}:
        return (positions * BOHR_TO_ANG) @ np.linalg.inv(cell_ang)
    if unit == "alat":
        metadata = load_metadata(data)
        alat_bohr = float(metadata["alat_bohr"])
        positions_ang = positions * (alat_bohr * BOHR_TO_ANG)
        return positions_ang @ np.linalg.inv(cell_ang)
    raise ValueError(f"Unsupported positions_unit={unit!r}")


def wrap_fractional(frac_pos: np.ndarray) -> np.ndarray:
    return np.asarray(frac_pos, dtype=float) - np.floor(np.asarray(frac_pos, dtype=float))


def infer_structure_tag(natoms: int) -> str:
    conventional_cells = natoms / 2.0
    n = round(conventional_cells ** (1.0 / 3.0))
    if 2 * n * n * n == natoms:
        return f"bcc_{n}x{n}x{n}"
    return f"bcc_{natoms}atoms"


def generate_unique_qe_species_labels(natoms: int) -> list[str]:
    if natoms < 1:
        raise ValueError("Need at least 1 atom to generate QE species labels.")
    if natoms > MAX_FAKE_SPECIES:
        raise ValueError(
            f"Cannot generate {natoms} unique 3-character QE labels with the current scheme; "
            f"maximum supported is {MAX_FAKE_SPECIES}."
        )

    labels: list[str] = []
    for atom_index in range(natoms):
        letter = FAKE_SPECIES_LETTERS[atom_index // 99]
        number = atom_index % 99 + 1
        labels.append(f"{letter}{number:02d}")
    return labels


def normalize_direction_vectors(
    vectors: np.ndarray,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    norms = np.linalg.norm(vectors, axis=1)
    bad = norms < 1.0e-14
    while np.any(bad):
        if rng is None:
            rng = np.random.default_rng()
        vectors[bad] = rng.normal(size=(np.sum(bad), 3))
        norms = np.linalg.norm(vectors, axis=1)
        bad = norms < 1.0e-14
    return vectors / norms[:, None]


def random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    quaternion = rng.normal(size=4)
    quaternion /= np.linalg.norm(quaternion)
    w, x, y, z = quaternion
    return np.array(
        [
            [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
            [2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
            [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)],
        ],
        dtype=float,
    )


def spins_from_directions(
    directions: np.ndarray,
    m_abs: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    spins_cart = m_abs * directions
    angle1 = np.degrees(np.arccos(np.clip(directions[:, 2], -1.0, 1.0)))
    angle2 = np.degrees(np.arctan2(directions[:, 1], directions[:, 0]))
    angle2 = np.where(angle2 < 0.0, angle2 + 360.0, angle2)
    net_m = np.sum(spins_cart, axis=0)
    return spins_cart, angle1, angle2, net_m


def generate_paramagnetic_spins_zero_net(
    natoms: int,
    m_abs: float = 0.35,
    seed: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if natoms < 2:
        raise ValueError("Need at least 2 atoms to build a zero-net paramagnetic spin pattern.")
    if natoms % 2 != 0:
        raise ValueError("The paired zero-net generator expects an even number of atoms.")

    rng = np.random.default_rng(seed)
    half = natoms // 2
    vectors = rng.normal(size=(half, 3))
    directions = normalize_direction_vectors(vectors, rng=rng)
    directions = np.vstack([directions, -directions])
    rng.shuffle(directions, axis=0)
    return spins_from_directions(directions, m_abs)


def generate_paramagnetic_spins_quasi_random_zero_net(
    natoms: int,
    m_abs: float = 0.35,
    seed: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if natoms < 2:
        raise ValueError("Need at least 2 atoms to build a zero-net paramagnetic spin pattern.")
    if natoms % 2 != 0:
        raise ValueError("The antipodal quasi-random generator expects an even number of atoms.")

    rng = np.random.default_rng(seed)
    half = natoms // 2
    indices = np.arange(half, dtype=float)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))

    z = 1.0 - 2.0 * (indices + 0.5) / half
    azimuth = golden_angle * indices
    radial = np.sqrt(np.clip(1.0 - z * z, 0.0, None))

    base_directions = np.column_stack(
        (
            radial * np.cos(azimuth),
            radial * np.sin(azimuth),
            z,
        )
    )
    directions = np.vstack([base_directions, -base_directions])
    directions = directions @ random_rotation_matrix(rng).T
    rng.shuffle(directions, axis=0)
    return spins_from_directions(directions, m_abs)


def generate_paramagnetic_spins(
    natoms: int,
    *,
    m_abs: float = 0.35,
    seed: int | None = None,
    mode: str = SPIN_MODE_RANDOM_ZERO_NET,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if mode == SPIN_MODE_RANDOM_ZERO_NET:
        return generate_paramagnetic_spins_zero_net(natoms, m_abs=m_abs, seed=seed)
    if mode == SPIN_MODE_QUASI_RANDOM_ZERO_NET:
        return generate_paramagnetic_spins_quasi_random_zero_net(natoms, m_abs=m_abs, seed=seed)
    raise ValueError(
        f"Unsupported spin mode {mode!r}. Supported modes: {', '.join(SUPPORTED_SPIN_MODES)}"
    )


def parse_qe_atom_resolved_spin_guess(
    qe_output_path: Path,
    *,
    natoms: int | None = None,
) -> QESpinGuess:
    qe_output_path = qe_output_path.resolve()
    text = qe_output_path.read_text()
    matches = list(QE_ATOM_SPIN_BLOCK_RE.finditer(text))
    if not matches:
        raise ValueError(
            "Could not find any atom-resolved magnetization blocks in the QE output. "
            "Expected blocks containing 'atom number', 'magnetization/charge', and 'polar coord.'."
        )

    atom_numbers = np.array([int(match.group(1)) for match in matches], dtype=int)
    if not np.array_equal(atom_numbers, np.arange(1, len(atom_numbers) + 1, dtype=int)):
        raise ValueError(
            "The atom-resolved QE spin blocks are incomplete or out of order. "
            f"Parsed atom numbers start as: {atom_numbers[:10].tolist()}"
        )
    if natoms is not None and len(atom_numbers) != natoms:
        raise ValueError(
            f"QE output contains {len(atom_numbers)} atom-resolved spin blocks, but the structure has {natoms} atoms."
        )

    charge_in_sphere = np.array([qe_float(match.group(5)) for match in matches], dtype=float)
    moment_cart_bohr = np.array(
        [[qe_float(match.group(group)) for group in (7, 8, 9)] for match in matches],
        dtype=float,
    )
    moment_over_charge_cart = np.array(
        [[qe_float(match.group(group)) for group in (10, 11, 12)] for match in matches],
        dtype=float,
    )
    local_moment_bohr = np.array([qe_float(match.group(13)) for match in matches], dtype=float)
    angle1 = np.array([qe_float(match.group(14)) for match in matches], dtype=float)
    angle2 = np.mod(np.array([qe_float(match.group(15)) for match in matches], dtype=float), 360.0)
    starting_magnetization = np.linalg.norm(moment_over_charge_cart, axis=1)

    return QESpinGuess(
        qe_output_path=qe_output_path,
        atom_numbers=atom_numbers,
        charge_in_sphere=charge_in_sphere,
        moment_cart_bohr=moment_cart_bohr,
        moment_over_charge_cart=moment_over_charge_cart,
        local_moment_bohr=local_moment_bohr,
        angle1=angle1,
        angle2=angle2,
        starting_magnetization=starting_magnetization,
    )


def write_spin_vectors_txt(spins_cart: np.ndarray, path: Path) -> None:
    lines = ["# mx my mz"]
    lines.extend(f"{mx:.10f} {my:.10f} {mz:.10f}" for mx, my, mz in spins_cart)
    path.write_text("\n".join(lines) + "\n")


def resolve_starting_magnetization_array(
    starting_magnetization: float | np.ndarray,
    natoms: int,
) -> np.ndarray:
    values = np.asarray(starting_magnetization, dtype=float)
    if values.ndim == 0:
        return np.full(natoms, float(values), dtype=float)
    values = values.reshape(-1)
    if values.shape != (natoms,):
        raise ValueError(
            f"Expected {natoms} starting_magnetization values, got shape {values.shape}."
        )
    return values


def write_qe_spin_parameters_txt(
    angle1: np.ndarray,
    angle2: np.ndarray,
    path: Path,
    starting_magnetization: float | np.ndarray,
) -> None:
    start_mags = resolve_starting_magnetization_array(starting_magnetization, len(angle1))
    lines = ["# QE noncollinear initial magnetic parameters"]
    for index, (theta, phi, start_mag) in enumerate(zip(angle1, angle2, start_mags), start=1):
        lines.append(f"starting_magnetization({index}) = {start_mag:.10f}")
        lines.append(f"angle1({index}) = {theta:.10f}")
        lines.append(f"angle2({index}) = {phi:.10f}")
    path.write_text("\n".join(lines) + "\n")


def find_card_start(lines: list[str], card_name: str) -> int:
    upper_name = card_name.upper()
    for index, line in enumerate(lines):
        if line.strip().upper().startswith(upper_name):
            return index
    raise ValueError(f"Could not find {card_name} card.")


def card_end_index(lines: list[str], start_index: int) -> int:
    index = start_index + 1
    while index < len(lines):
        stripped = lines[index].strip().upper()
        if stripped.startswith(
            (
                "&CONTROL",
                "&SYSTEM",
                "&ELECTRONS",
                "&IONS",
                "&CELL",
                "ATOMIC_SPECIES",
                "ATOMIC_POSITIONS",
                "ATOMIC_VELOCITIES",
                "K_POINTS",
                "CELL_PARAMETERS",
                "CONSTRAINTS",
                "OCCUPATIONS",
                "ATOMIC_FORCES",
                "SOLVENTS",
                "HUBBARD",
            )
        ):
            break
        index += 1
    return index


def extract_card_labels(qe_text: str, card_name: str, min_fields: int) -> list[str]:
    lines = qe_text.splitlines()
    start = find_card_start(lines, card_name)
    end = card_end_index(lines, start)
    labels: list[str] = []
    for line in lines[start + 1 : end]:
        stripped = line.strip()
        if not stripped or stripped.startswith(("!", "#")):
            continue
        fields = stripped.split()
        if len(fields) < min_fields:
            continue
        labels.append(fields[0])
    return labels


def validate_label_consistency(
    species_labels: list[str],
    position_labels: list[str],
    velocity_labels: list[str],
) -> None:
    natoms = len(position_labels)
    if len(species_labels) != natoms:
        raise ValueError(f"Expected {natoms} species labels, got {len(species_labels)}.")
    if len(set(species_labels)) != natoms:
        raise ValueError("Fake QE species labels are not unique.")

    species_set = set(species_labels)
    missing_positions = [label for label in position_labels if label not in species_set]
    if missing_positions:
        raise ValueError(
            "ATOMIC_POSITIONS uses labels missing from ATOMIC_SPECIES: "
            + ", ".join(missing_positions[:10])
        )

    missing_velocities = [label for label in velocity_labels if label not in species_set]
    if missing_velocities:
        raise ValueError(
            "ATOMIC_VELOCITIES uses labels missing from ATOMIC_SPECIES: "
            + ", ".join(missing_velocities[:10])
        )

    if len(velocity_labels) != natoms:
        raise ValueError(f"Expected {natoms} velocity labels, got {len(velocity_labels)}.")

    if position_labels != velocity_labels:
        for atom_index, (pos_label, vel_label) in enumerate(zip(position_labels, velocity_labels), start=1):
            if pos_label != vel_label:
                raise ValueError(
                    f"ATOMIC_POSITIONS and ATOMIC_VELOCITIES labels differ at atom {atom_index}: "
                    f"{pos_label!r} != {vel_label!r}"
                )
        raise ValueError("ATOMIC_POSITIONS and ATOMIC_VELOCITIES labels differ.")


def print_generation_diagnostic(
    natoms: int,
    species_labels: list[str],
    position_labels: list[str],
    velocity_labels: list[str],
) -> None:
    print(f"Number of atoms: {natoms}")
    print(f"Number of species: {len(species_labels)}")
    print(f"First 10 labels: {species_labels[:10]}")
    print(f"Last 10 labels: {species_labels[-10:]}")
    print(f"Position and velocity labels match exactly: {position_labels == velocity_labels}")


def build_noncollinear_md_input(
    frac_positions: np.ndarray,
    cell_ang: np.ndarray,
    *,
    species_labels: list[str],
    prefix: str,
    pseudo_dir: str,
    pseudo_file: str,
    temperature_k: float,
    dt_au: float,
    nstep: int,
    ecutwfc: float,
    ecutrho: float,
    degauss: float,
    k_grid: tuple[int, int, int],
    angle1: np.ndarray,
    angle2: np.ndarray,
    starting_magnetization: float | np.ndarray,
    constrained_magnetization: bool,
    lambda_value: float,
    mixing_beta: float,
    nosym: bool,
) -> str:
    natoms = len(frac_positions)
    if len(species_labels) != natoms:
        raise ValueError(f"Expected {natoms} species labels, got {len(species_labels)}.")
    start_mags = resolve_starting_magnetization_array(starting_magnetization, natoms)
    kx, ky, kz = k_grid

    lines = [
        "&control",
        "   calculation='md'",
        "   restart_mode='from_scratch'",
        f"   prefix='{prefix}'",
        f"   pseudo_dir='{pseudo_dir}'",
        "   outdir='./md'",
        f"   dt={dt_au:.4f}d0",
        "   tstress=.true.",
        "   tprnfor=.true.",
        f"   nstep={nstep}",
        "/",
        " &system",
        "    ibrav = 0",
        f"    nat = {natoms}, ntyp = {natoms},",
        f"    ecutwfc={ecutwfc:g}",
        f"    ecutrho={ecutrho:g}",
        f"    occupations='smearing',smearing='m-v',degauss={degauss:.2f}",
        f"    nosym={'.true.' if nosym else '.false.'}",
        "    noncolin=.true.",
        "    lspinorb=.false.",
    ]
    if constrained_magnetization:
        lines.append("    constrained_magnetization='atomic direction'")
        lines.append(f"    lambda={lambda_value:g}")
    for index, (theta, phi, start_mag) in enumerate(zip(angle1, angle2, start_mags), start=1):
        lines.append(f"    starting_magnetization({index})={start_mag:.10f}")
        lines.append(f"    angle1({index})={theta:.10f}")
        lines.append(f"    angle2({index})={phi:.10f}")

    lines.extend(
        [
            "/",
            " &electrons",
            "    conv_thr=1.0d-4",
            f"    mixing_beta={mixing_beta:.4f}d0",
            "/",
            " &ions",
            "    pot_extrapolation='second_order'",
            "    wfc_extrapolation='second_order'",
            "    ion_temperature='svr'",
            "    ion_velocities='from_input'",
            f"    tempw={temperature_k:.1f}",
            "    nraise=20",
            "/",
            "ATOMIC_SPECIES",
        ]
    )
    lines.extend(f"{label} 55.845 {pseudo_file}" for label in species_labels)
    lines.extend(["", "ATOMIC_POSITIONS (crystal)"])
    lines.extend(
        f"{label} {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}"
        for label, pos in zip(species_labels, frac_positions)
    )
    lines.extend(["", "CELL_PARAMETERS angstrom"])
    lines.extend(f"{row[0]:.12f} {row[1]:.12f} {row[2]:.12f}" for row in cell_ang)
    lines.extend(["", "K_POINTS automatic", f"{kx} {ky} {kz} 0 0 0", ""])
    return "\n".join(lines)


def prepare_bcc_qe_input(
    *,
    dataset_dir: Path,
    npz_path: Path | None = None,
    output_dir: Path,
    frame_index: int = -1,
    temperature_k: float | None = None,
    velocity_seed: int | None = None,
    spin_seed: int | None = None,
    spin_mode: str = SPIN_MODE_RANDOM_ZERO_NET,
    qe_spin_output_path: Path | None = None,
    m_abs: float = 0.35,
    pseudo_dir: str = ".",
    pseudo_file: str = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
    dt_au: float = 20.670,
    nstep: int = 400,
    ecutwfc: float = 71.0,
    ecutrho: float = 496.0,
    degauss: float = 0.02,
    k_grid: tuple[int, int, int] = (4, 4, 4),
    constrained_magnetization: bool = True,
    lambda_value: float = 0.2,
    mixing_beta: float = 0.01,
    remove_com_drift: bool = True,
    rescale_exact: bool = True,
    nosym: bool = True,
) -> PreparationResult:
    chosen_npz = npz_path.resolve() if npz_path is not None else infer_latest_bcc_npz(dataset_dir.resolve())
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    data = np.load(chosen_npz, allow_pickle=False)
    positions = np.asarray(data["positions"], dtype=float)
    natoms = int(positions.shape[1])
    resolved_index = frame_index if frame_index >= 0 else positions.shape[0] + frame_index
    if resolved_index < 0 or resolved_index >= positions.shape[0]:
        raise IndexError(f"Frame index {frame_index} is out of range for {chosen_npz.name}")

    fallback_cell = fixed_cell_angstrom(data)
    cell_ang = frame_cell_angstrom(data, resolved_index, fallback_cell)
    frac_positions = wrap_fractional(frame_positions_fractional(data, resolved_index, cell_ang))
    target_temperature = float(temperature_k) if temperature_k is not None else parse_temperature_from_name(chosen_npz)

    structure_tag = infer_structure_tag(natoms)
    prefix = f"Fe_{structure_tag}_latest_noncollinear_paramagnetic"
    base_input = output_dir / f"{prefix}_{int(round(target_temperature))}K_base.in"
    final_input = output_dir / f"{prefix}_{int(round(target_temperature))}K_with_velocities.in"
    velocity_block_path = output_dir / f"{prefix}_atomic_velocities_{int(round(target_temperature))}K.txt"
    spin_vectors_path = output_dir / f"{prefix}_spin_vectors.txt"
    spin_parameters_path = output_dir / f"{prefix}_qe_spin_parameters.txt"
    species_labels = generate_unique_qe_species_labels(natoms)

    qe_spin_output = None if qe_spin_output_path is None else qe_spin_output_path.resolve()
    if qe_spin_output is None:
        spins_cart, angle1, angle2, net_m = generate_paramagnetic_spins(
            natoms,
            m_abs=m_abs,
            seed=spin_seed,
            mode=spin_mode,
        )
        starting_magnetization = np.full(natoms, float(m_abs), dtype=float)
        spin_source = f"generated:{spin_mode}"
    else:
        spin_guess = parse_qe_atom_resolved_spin_guess(qe_spin_output, natoms=natoms)
        spins_cart = spin_guess.moment_cart_bohr
        angle1 = spin_guess.angle1
        angle2 = spin_guess.angle2
        starting_magnetization = spin_guess.starting_magnetization
        net_m = np.sum(spins_cart, axis=0)
        spin_source = f"qe_output:{qe_spin_output.name}"

    write_spin_vectors_txt(spins_cart, spin_vectors_path)
    write_qe_spin_parameters_txt(
        angle1,
        angle2,
        spin_parameters_path,
        starting_magnetization=starting_magnetization,
    )

    base_text = build_noncollinear_md_input(
        frac_positions,
        cell_ang,
        species_labels=species_labels,
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        pseudo_file=pseudo_file,
        temperature_k=target_temperature,
        dt_au=dt_au,
        nstep=nstep,
        ecutwfc=ecutwfc,
        ecutrho=ecutrho,
        degauss=degauss,
        k_grid=k_grid,
        angle1=angle1,
        angle2=angle2,
        starting_magnetization=starting_magnetization,
        constrained_magnetization=constrained_magnetization,
        lambda_value=lambda_value,
        mixing_beta=mixing_beta,
        nosym=nosym,
    )
    base_species_labels = extract_card_labels(base_text, "ATOMIC_SPECIES", min_fields=3)
    base_position_labels = extract_card_labels(base_text, "ATOMIC_POSITIONS", min_fields=4)
    validate_label_consistency(base_species_labels, base_position_labels, base_position_labels)
    base_input.write_text(base_text)

    labels, masses_amu = read_qe_species_sequence(base_input)
    velocities_au = sample_maxwell_velocities_au(masses_amu, target_temperature, velocity_seed)
    if remove_com_drift:
        velocities_au = remove_center_of_mass_drift_au(velocities_au, masses_amu)
    if rescale_exact:
        velocities_au = rescale_to_temperature_au(
            velocities_au,
            masses_amu,
            target_temperature,
            remove_com=remove_com_drift,
    )
    measured_temperature = temperature_from_velocities_au(
        velocities_au,
        masses_amu,
        remove_com=remove_com_drift,
    )
    velocity_block = format_atomic_velocities_card(labels, velocities_au)
    velocity_labels = extract_card_labels(velocity_block, "ATOMIC_VELOCITIES", min_fields=4)
    validate_label_consistency(base_species_labels, base_position_labels, velocity_labels)
    velocity_block_path.write_text(velocity_block)
    final_text = replace_or_append_atomic_velocities(base_text, velocity_block)
    final_species_labels = extract_card_labels(final_text, "ATOMIC_SPECIES", min_fields=3)
    final_position_labels = extract_card_labels(final_text, "ATOMIC_POSITIONS", min_fields=4)
    final_velocity_labels = extract_card_labels(final_text, "ATOMIC_VELOCITIES", min_fields=4)
    validate_label_consistency(final_species_labels, final_position_labels, final_velocity_labels)
    print_generation_diagnostic(natoms, final_species_labels, final_position_labels, final_velocity_labels)
    final_input.write_text(final_text)

    return PreparationResult(
        npz_path=chosen_npz,
        output_dir=output_dir,
        frame_index=resolved_index,
        natoms=natoms,
        ntypes=natoms,
        structure_tag=structure_tag,
        target_temperature_k=target_temperature,
        measured_velocity_temperature_k=float(measured_temperature),
        qe_input_base=base_input,
        qe_input_final=final_input,
        velocity_block=velocity_block_path,
        spin_vectors=spin_vectors_path,
        spin_parameters=spin_parameters_path,
        net_magnetization=(float(net_m[0]), float(net_m[1]), float(net_m[2])),
        spin_mode=spin_mode,
        spin_source=spin_source,
        qe_spin_output_path=qe_spin_output,
        mean_starting_magnetization=float(np.mean(starting_magnetization)),
        min_starting_magnetization=float(np.min(starting_magnetization)),
        max_starting_magnetization=float(np.max(starting_magnetization)),
        first_labels=tuple(final_species_labels[:10]),
        last_labels=tuple(final_species_labels[-10:]),
        position_velocity_labels_match=(final_position_labels == final_velocity_labels),
    )
