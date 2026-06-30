#!/usr/bin/env python3
"""Trace BCC/FCC/HCP character through an MD trajectory.

The script reads either:
- LAMMPS dump trajectories (`.lammpstrj`, `.dump`, `.traj`)
- QE-style NPZ archives used elsewhere in this repository

For each frame it computes a local bond-angle fingerprint around every atom and
compares it against ideal BCC, FCC, and HCP templates. A normalized
neighbor-shell distance fingerprint is added as a tie-breaker, which improves
FCC/HCP separation while keeping the score driven by local geometry. The result
is written as:
- `phase_trace.csv`: per-frame scores and phase fractions
- `phase_transitions.csv`: likely phase changes after temporal smoothing
- `phase_trace.png`: score/q_l/fraction plot versus MD step
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
import tempfile

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-codex"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


BOHR_TO_ANG = 0.529177210903
PHASES = ("bcc", "fcc", "hcp")
ANGLE_WEIGHT = 0.65
RADIAL_WEIGHT = 0.35
RADIAL_SIGMA = 0.08
PHASE_COLORS = {
    "bcc": "#1f6f5b",
    "fcc": "#c24b2a",
    "hcp": "#2f5597",
    "unknown": "#666666",
}


@dataclass(frozen=True)
class TrajectoryFrame:
    step: int
    positions: np.ndarray
    cell: np.ndarray


@dataclass(frozen=True)
class PhaseTemplate:
    name: str
    angle_neighbor_count: int
    radial_neighbor_count: int
    first_shell_count: int
    ideal_cosines: np.ndarray
    ideal_distances: np.ndarray


@dataclass(frozen=True)
class TransitionEvent:
    frame_index: int
    step: int
    from_phase: str
    to_phase: str


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze an MD trajectory and trace BCC/FCC/HCP character as a function of MD step "
            "using local bond-angle fingerprints."
        )
    )
    parser.add_argument("input", type=Path, help="Trajectory path: LAMMPS dump or QE-style NPZ.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for CSV/plot outputs. Default: <input_stem>_phase_trace next to the input file.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Analyze every Nth frame. Default: 1.",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional cap on the number of analyzed frames.",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=5,
        help="Moving-average window used before transition detection. Default: 5.",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=0.18,
        help="Similarity width in cosine-space RMSD. Lower values make classification stricter. Default: 0.18.",
    )
    parser.add_argument(
        "--atom-score-threshold",
        type=float,
        default=0.45,
        help="Atoms below this best-match score are labeled unknown in the phase fractions. Default: 0.45.",
    )
    parser.add_argument(
        "--atom-margin-threshold",
        type=float,
        default=0.03,
        help="Atoms are labeled unknown when the best and second-best scores differ by less than this margin. Default: 0.03.",
    )
    parser.add_argument(
        "--frame-score-threshold",
        type=float,
        default=0.50,
        help="Frames below this best smoothed score are labeled unknown for transition detection. Default: 0.50.",
    )
    parser.add_argument(
        "--frame-margin-threshold",
        type=float,
        default=0.02,
        help="Frames are labeled unknown when the best and second-best smoothed scores differ by less than this margin. Default: 0.02.",
    )
    parser.add_argument(
        "--transition-persistence",
        type=int,
        default=3,
        help="A new phase must persist for at least this many analyzed frames to count as a transition. Default: 3.",
    )
    parser.add_argument(
        "--ql-orders",
        type=int,
        nargs="+",
        default=[2, 4, 6],
        help="Steinhardt q_l orders to track. Default: 2 4 6.",
    )
    parser.add_argument(
        "--ql-neighbor-count",
        type=int,
        default=12,
        help="Number of nearest neighbors used to compute q_l. Default: 12.",
    )
    return parser.parse_args(argv)


def scalar_string(value) -> str:
    array = np.asarray(value)
    return str(array.item() if array.shape == () else value)


def load_metadata(data) -> dict[str, object]:
    if "metadata_json" not in data.files:
        return {}
    payload = data["metadata_json"]
    if isinstance(payload, np.ndarray):
        payload = payload.item()
    if isinstance(payload, bytes):
        payload = payload.decode("utf-8")
    return json.loads(str(payload))


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


def frame_positions_angstrom(data, frame_index: int, cell_ang: np.ndarray) -> np.ndarray:
    positions = np.asarray(data["positions"], dtype=float)[frame_index]
    units = np.asarray(data["positions_unit"])
    unit = str(units[frame_index]).lower() if units.ndim > 0 else scalar_string(units).lower()

    if unit == "crystal":
        return positions @ cell_ang
    if unit in {"angstrom", "ang"}:
        return positions
    if unit in {"bohr", "a.u.", "au"}:
        return positions * BOHR_TO_ANG
    if unit == "alat":
        metadata = load_metadata(data)
        alat_bohr = float(metadata["alat_bohr"])
        return positions * (alat_bohr * BOHR_TO_ANG)
    raise ValueError(f"Unsupported positions_unit={unit!r}")


def infer_npz_steps(data, nframes: int) -> np.ndarray:
    for key in ("steps", "step", "step_indices", "timesteps", "timestep", "md_steps"):
        if key in data.files:
            values = np.asarray(data[key]).reshape(-1)
            if len(values) == nframes:
                return np.rint(values).astype(int)
    return np.arange(nframes, dtype=int)


def wrap_positions(positions: np.ndarray, cell: np.ndarray) -> np.ndarray:
    inv_cell = np.linalg.inv(cell)
    frac = positions @ inv_cell
    frac -= np.floor(frac)
    return frac @ cell


def read_npz_frames(path: Path) -> list[TrajectoryFrame]:
    with np.load(path, allow_pickle=True) as data:
        positions_array = np.asarray(data["positions"])
        if positions_array.ndim != 3 or positions_array.shape[-1] != 3:
            raise ValueError(f"Unexpected positions array shape in {path}: {positions_array.shape}")

        nframes = positions_array.shape[0]
        steps = infer_npz_steps(data, nframes)
        base_cell = fixed_cell_angstrom(data)

        frames: list[TrajectoryFrame] = []
        for frame_index in range(nframes):
            cell = frame_cell_angstrom(data, frame_index, base_cell)
            positions = frame_positions_angstrom(data, frame_index, cell)
            positions = wrap_positions(np.asarray(positions, dtype=float), np.asarray(cell, dtype=float))
            frames.append(TrajectoryFrame(step=int(steps[frame_index]), positions=positions, cell=cell))
        return frames


def parse_box_bounds(header: str, lines: list[str]) -> tuple[np.ndarray, np.ndarray]:
    tokens = header.split()[3:]
    has_tilt = any(token in {"xy", "xz", "yz"} for token in tokens)

    if has_tilt:
        xlo_bound, xhi_bound, xy = map(float, lines[0].split()[:3])
        ylo_bound, yhi_bound, xz = map(float, lines[1].split()[:3])
        zlo_bound, zhi_bound, yz = map(float, lines[2].split()[:3])

        xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
        xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
        ylo = ylo_bound - min(0.0, yz)
        yhi = yhi_bound - max(0.0, yz)
        zlo = zlo_bound
        zhi = zhi_bound

        origin = np.array([xlo, ylo, zlo], dtype=float)
        cell = np.array(
            [
                [xhi - xlo, 0.0, 0.0],
                [xy, yhi - ylo, 0.0],
                [xz, yz, zhi - zlo],
            ],
            dtype=float,
        )
        return origin, cell

    bounds = [list(map(float, line.split()[:2])) for line in lines]
    xlo, xhi = bounds[0]
    ylo, yhi = bounds[1]
    zlo, zhi = bounds[2]
    origin = np.array([xlo, ylo, zlo], dtype=float)
    cell = np.array(
        [
            [xhi - xlo, 0.0, 0.0],
            [0.0, yhi - ylo, 0.0],
            [0.0, 0.0, zhi - zlo],
        ],
        dtype=float,
    )
    return origin, cell


def read_lammps_dump_frames(path: Path) -> list[TrajectoryFrame]:
    frames: list[TrajectoryFrame] = []
    with path.open("r", encoding="utf-8") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                raise ValueError(f"Expected 'ITEM: TIMESTEP' in {path}, got: {line.strip()}")

            step = int(handle.readline().strip())

            line = handle.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"Expected 'ITEM: NUMBER OF ATOMS' in {path}, got: {line.strip()}")
            natoms = int(handle.readline().strip())

            box_header = handle.readline().strip()
            if not box_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Expected 'ITEM: BOX BOUNDS' in {path}, got: {box_header}")
            box_lines = [handle.readline().strip() for _ in range(3)]
            origin, cell = parse_box_bounds(box_header, box_lines)

            atom_header = handle.readline().strip()
            if not atom_header.startswith("ITEM: ATOMS"):
                raise ValueError(f"Expected 'ITEM: ATOMS' in {path}, got: {atom_header}")
            columns = atom_header.split()[2:]
            column_map = {name: index for index, name in enumerate(columns)}

            atom_rows = [handle.readline().split() for _ in range(natoms)]
            atom_data = np.asarray(atom_rows, dtype=float)
            if atom_data.shape[0] != natoms:
                raise ValueError(f"Failed to read {natoms} atoms from {path}")

            if "id" in column_map:
                order = np.argsort(atom_data[:, column_map["id"]])
                atom_data = atom_data[order]

            if all(key in column_map for key in ("x", "y", "z")):
                positions = atom_data[:, [column_map["x"], column_map["y"], column_map["z"]]] - origin
            elif all(key in column_map for key in ("xu", "yu", "zu")):
                positions = atom_data[:, [column_map["xu"], column_map["yu"], column_map["zu"]]] - origin
            elif all(key in column_map for key in ("xs", "ys", "zs")):
                frac = atom_data[:, [column_map["xs"], column_map["ys"], column_map["zs"]]]
                positions = frac @ cell
            elif all(key in column_map for key in ("xsu", "ysu", "zsu")):
                frac = atom_data[:, [column_map["xsu"], column_map["ysu"], column_map["zsu"]]]
                positions = frac @ cell
            else:
                raise ValueError(
                    f"Could not find supported coordinates in {path}. "
                    "Need one of: x/y/z, xu/yu/zu, xs/ys/zs, xsu/ysu/zsu."
                )

            frames.append(TrajectoryFrame(step=step, positions=wrap_positions(positions, cell), cell=cell))
    return frames


def coerce_path(path: str | os.PathLike[str] | Path) -> Path:
    return Path(path).expanduser()


def read_frames(path: str | os.PathLike[str] | Path) -> list[TrajectoryFrame]:
    path = coerce_path(path)
    suffix = path.suffix.lower()
    if suffix == ".npz":
        return read_npz_frames(path)
    if suffix in {".lammpstrj", ".dump", ".traj"}:
        return read_lammps_dump_frames(path)
    raise ValueError(f"Unsupported input format for {path}. Use .npz or a LAMMPS dump file.")


def ideal_neighbor_vectors(cell: np.ndarray, basis: np.ndarray, neighbor_count: int) -> tuple[np.ndarray, np.ndarray]:
    vectors: list[np.ndarray] = []
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                shift = np.array([i, j, k], dtype=float)
                for basis_index, frac in enumerate(basis):
                    if basis_index == 0 and i == 0 and j == 0 and k == 0:
                        continue
                    vectors.append((shift + frac) @ cell)

    points = np.asarray(vectors, dtype=float)
    distances = np.linalg.norm(points, axis=1)
    order = np.argsort(distances)[:neighbor_count]
    return points[order], distances[order]


def ideal_neighbor_cosines(cell: np.ndarray, basis: np.ndarray, neighbor_count: int) -> np.ndarray:
    neighbors, distances = ideal_neighbor_vectors(cell, basis, neighbor_count)
    unit = neighbors / distances[:, None]
    cosine_matrix = unit @ unit.T
    tri_upper = np.triu_indices(neighbor_count, k=1)
    return np.sort(cosine_matrix[tri_upper])


def build_phase_templates() -> dict[str, PhaseTemplate]:
    hcp_c_over_a = math.sqrt(8.0 / 3.0)
    template_specs = {
        "bcc": (
            np.eye(3),
            np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], dtype=float),
            8,
            14,
            8,
        ),
        "fcc": (
            np.eye(3),
            np.array(
                [
                    [0.0, 0.0, 0.0],
                    [0.0, 0.5, 0.5],
                    [0.5, 0.0, 0.5],
                    [0.5, 0.5, 0.0],
                ],
                dtype=float,
            ),
            12,
            20,
            12,
        ),
        "hcp": (
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [-0.5, math.sqrt(3.0) / 2.0, 0.0],
                    [0.0, 0.0, hcp_c_over_a],
                ],
                dtype=float,
            ),
            np.array([[0.0, 0.0, 0.0], [2.0 / 3.0, 1.0 / 3.0, 0.5]], dtype=float),
            12,
            20,
            12,
        ),
    }

    templates: dict[str, PhaseTemplate] = {}
    for name, (cell, basis, angle_neighbor_count, radial_neighbor_count, first_shell_count) in template_specs.items():
        _, ideal_radial_distances = ideal_neighbor_vectors(cell, basis, radial_neighbor_count)
        templates[name] = PhaseTemplate(
            name=name,
            angle_neighbor_count=angle_neighbor_count,
            radial_neighbor_count=radial_neighbor_count,
            first_shell_count=first_shell_count,
            ideal_cosines=ideal_neighbor_cosines(cell, basis, angle_neighbor_count),
            ideal_distances=np.sort(ideal_radial_distances / np.mean(ideal_radial_distances[:first_shell_count])),
        )
    return templates


def moving_average(values: np.ndarray, window: int) -> np.ndarray:
    if values.ndim == 2 and values.shape[1] == 0:
        return values.copy()
    if window <= 1:
        return values.copy()
    kernel = np.ones(window, dtype=float) / float(window)
    return np.vstack([np.convolve(values[:, index], kernel, mode="same") for index in range(values.shape[1])]).T


def neighbor_vectors(positions: np.ndarray, cell: np.ndarray, max_neighbors: int) -> tuple[np.ndarray, np.ndarray]:
    natoms = positions.shape[0]
    if natoms - 1 < max_neighbors:
        raise ValueError(
            f"Trajectory frame has {natoms} atoms, but at least {max_neighbors + 1} are needed "
            "to evaluate BCC/FCC/HCP neighbor fingerprints."
        )

    inv_cell = np.linalg.inv(cell)
    frac = positions @ inv_cell
    frac -= np.floor(frac)

    delta_frac = frac[None, :, :] - frac[:, None, :]
    delta_frac -= np.rint(delta_frac)
    delta_cart = delta_frac @ cell
    distances = np.linalg.norm(delta_cart, axis=-1)
    np.fill_diagonal(distances, np.inf)

    nearest = np.argpartition(distances, kth=max_neighbors - 1, axis=1)[:, :max_neighbors]
    nearest_distances = np.take_along_axis(distances, nearest, axis=1)
    order = np.argsort(nearest_distances, axis=1)
    nearest = np.take_along_axis(nearest, order, axis=1)
    nearest_distances = np.take_along_axis(nearest_distances, order, axis=1)

    atom_indices = np.arange(natoms)[:, None]
    return delta_cart[atom_indices, nearest], nearest_distances


def sorted_pairwise_cosines(vectors: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(vectors, axis=2)
    unit = vectors / np.clip(norms[:, :, None], 1.0e-12, None)
    cosine_matrices = np.einsum("nik,njk->nij", unit, unit)
    tri_upper = np.triu_indices(vectors.shape[1], k=1)
    return np.sort(cosine_matrices[:, tri_upper[0], tri_upper[1]], axis=1)


def steinhardt_ql_for_vectors(vectors: np.ndarray, orders: tuple[int, ...]) -> np.ndarray:
    if not orders:
        return np.zeros((vectors.shape[0], 0), dtype=float)

    norms = np.linalg.norm(vectors, axis=2)
    unit = vectors / np.clip(norms[:, :, None], 1.0e-12, None)
    cosine_matrices = np.einsum("nik,njk->nij", unit, unit)

    ql_values = np.zeros((vectors.shape[0], len(orders)), dtype=float)
    for order_index, order in enumerate(orders):
        coeffs = np.zeros(order + 1, dtype=float)
        coeffs[order] = 1.0
        legendre_values = np.polynomial.legendre.legval(cosine_matrices, coeffs)
        ql_squared = np.mean(legendre_values, axis=(1, 2))
        ql_values[:, order_index] = np.sqrt(np.clip(ql_squared, 0.0, None))

    return ql_values


def phase_scores_and_custom_ql_for_frame(
    positions: np.ndarray,
    cell: np.ndarray,
    templates: dict[str, PhaseTemplate],
    sigma: float,
    ql_orders: tuple[int, ...],
    ql_neighbor_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    max_neighbors = max(max(template.radial_neighbor_count for template in templates.values()), ql_neighbor_count)
    local_vectors, local_distances = neighbor_vectors(positions, cell, max_neighbors=max_neighbors)

    cosine_cache: dict[int, np.ndarray] = {}
    scores = np.zeros((positions.shape[0], len(PHASES)), dtype=float)

    for phase_index, phase_name in enumerate(PHASES):
        template = templates[phase_name]
        angle_count = template.angle_neighbor_count
        radial_count = template.radial_neighbor_count

        if angle_count not in cosine_cache:
            cosine_cache[angle_count] = sorted_pairwise_cosines(local_vectors[:, :angle_count, :])
        angle_rmse = np.sqrt(np.mean((cosine_cache[angle_count] - template.ideal_cosines[None, :]) ** 2, axis=1))
        angle_score = 1.0 / (1.0 + (angle_rmse / sigma) ** 2)

        radial_distances = np.sort(
            local_distances[:, :radial_count] / np.mean(local_distances[:, : template.first_shell_count], axis=1, keepdims=True),
            axis=1,
        )
        radial_rmse = np.sqrt(np.mean((radial_distances - template.ideal_distances[None, :]) ** 2, axis=1))
        radial_score = 1.0 / (1.0 + (radial_rmse / RADIAL_SIGMA) ** 2)

        scores[:, phase_index] = ANGLE_WEIGHT * angle_score + RADIAL_WEIGHT * radial_score

    ql_values = steinhardt_ql_for_vectors(local_vectors[:, :ql_neighbor_count, :], ql_orders)
    return scores, ql_values


def dominant_labels(smoothed_scores: np.ndarray, threshold: float, margin_threshold: float) -> list[str]:
    best_indices = np.argmax(smoothed_scores, axis=1)
    best_scores = np.max(smoothed_scores, axis=1)
    sorted_scores = np.sort(smoothed_scores, axis=1)
    margins = sorted_scores[:, -1] - sorted_scores[:, -2]
    labels: list[str] = []
    for best_index, best_score, margin in zip(best_indices, best_scores, margins):
        labels.append(
            PHASES[int(best_index)] if float(best_score) >= threshold and float(margin) >= margin_threshold else "unknown"
        )
    return labels


def detect_transitions(
    steps: np.ndarray,
    labels: list[str],
    persistence: int,
) -> list[TransitionEvent]:
    if not labels:
        return []

    events: list[TransitionEvent] = []
    current_phase = labels[0]
    frame_index = 1

    while frame_index < len(labels):
        if labels[frame_index] == current_phase:
            frame_index += 1
            continue

        candidate = labels[frame_index]
        run_end = frame_index + 1
        while run_end < len(labels) and labels[run_end] == candidate:
            run_end += 1

        if run_end - frame_index >= persistence:
            events.append(
                TransitionEvent(
                    frame_index=frame_index,
                    step=int(steps[frame_index]),
                    from_phase=current_phase,
                    to_phase=candidate,
                )
            )
            current_phase = candidate

        frame_index = run_end

    return events


def default_output_dir(input_path: str | os.PathLike[str] | Path) -> Path:
    input_path = coerce_path(input_path)
    return input_path.parent / f"{input_path.stem}_phase_trace"


def analyze_frames(
    frames: list[TrajectoryFrame],
    templates: dict[str, PhaseTemplate],
    sigma: float,
    atom_score_threshold: float,
    atom_margin_threshold: float,
    ql_orders: tuple[int, ...] = (2, 4, 6),
    ql_neighbor_count: int = 12,
) -> dict[str, np.ndarray]:
    mean_scores: list[np.ndarray] = []
    mean_ql_values: list[np.ndarray] = []
    phase_fractions: list[np.ndarray] = []
    steps: list[int] = []

    for frame in frames:
        scores, ql_values = phase_scores_and_custom_ql_for_frame(
            frame.positions,
            frame.cell,
            templates=templates,
            sigma=sigma,
            ql_orders=ql_orders,
            ql_neighbor_count=ql_neighbor_count,
        )
        mean_scores.append(np.mean(scores, axis=0))
        mean_ql_values.append(np.mean(ql_values, axis=0))

        best_indices = np.argmax(scores, axis=1)
        best_scores = np.max(scores, axis=1)
        sorted_scores = np.sort(scores, axis=1)
        margins = sorted_scores[:, -1] - sorted_scores[:, -2]
        fractions = np.zeros(4, dtype=float)
        unknown_mask = (best_scores < atom_score_threshold) | (margins < atom_margin_threshold)
        fractions[3] = float(np.mean(unknown_mask))
        for phase_index in range(len(PHASES)):
            phase_mask = (best_indices == phase_index) & (~unknown_mask)
            fractions[phase_index] = float(np.mean(phase_mask))
        phase_fractions.append(fractions)
        steps.append(frame.step)

    return {
        "steps": np.asarray(steps, dtype=int),
        "mean_scores": np.asarray(mean_scores, dtype=float),
        "mean_ql_values": np.asarray(mean_ql_values, dtype=float),
        "phase_fractions": np.asarray(phase_fractions, dtype=float),
    }


def write_phase_trace_csv(
    path: Path,
    steps: np.ndarray,
    mean_scores: np.ndarray,
    smoothed_scores: np.ndarray,
    mean_ql_values: np.ndarray,
    smoothed_ql_values: np.ndarray,
    ql_orders: tuple[int, ...],
    phase_fractions: np.ndarray,
    labels: list[str],
    transition_steps: set[int],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        ql_headers = [f"mean_q{order}" for order in ql_orders] + [f"smoothed_q{order}" for order in ql_orders]
        writer.writerow(
            [
                "frame_index",
                "step",
                "mean_bcc_score",
                "mean_fcc_score",
                "mean_hcp_score",
                "smoothed_bcc_score",
                "smoothed_fcc_score",
                "smoothed_hcp_score",
                *ql_headers,
                "bcc_fraction",
                "fcc_fraction",
                "hcp_fraction",
                "unknown_fraction",
                "dominant_phase",
                "transition_here",
            ]
        )
        for frame_index in range(len(steps)):
            writer.writerow(
                [
                    frame_index,
                    int(steps[frame_index]),
                    *[f"{value:.8f}" for value in mean_scores[frame_index]],
                    *[f"{value:.8f}" for value in smoothed_scores[frame_index]],
                    *[f"{value:.8f}" for value in mean_ql_values[frame_index]],
                    *[f"{value:.8f}" for value in smoothed_ql_values[frame_index]],
                    *[f"{value:.8f}" for value in phase_fractions[frame_index]],
                    labels[frame_index],
                    "yes" if int(steps[frame_index]) in transition_steps else "no",
                ]
            )


def write_transition_csv(path: Path, events: list[TransitionEvent]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["frame_index", "step", "from_phase", "to_phase"])
        for event in events:
            writer.writerow([event.frame_index, event.step, event.from_phase, event.to_phase])


def plot_phase_trace(
    path: Path,
    input_name: str,
    steps: np.ndarray,
    mean_scores: np.ndarray,
    smoothed_scores: np.ndarray,
    mean_ql_values: np.ndarray,
    smoothed_ql_values: np.ndarray,
    ql_orders: tuple[int, ...],
    phase_fractions: np.ndarray,
    events: list[TransitionEvent],
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(11.5, 10.5), sharex=True, constrained_layout=True)

    score_ax, ql_ax, fraction_ax = axes
    for phase_index, phase_name in enumerate(PHASES):
        color = PHASE_COLORS[phase_name]
        score_ax.plot(steps, mean_scores[:, phase_index], color=color, alpha=0.22, linewidth=1.0)
        score_ax.plot(
            steps,
            smoothed_scores[:, phase_index],
            color=color,
            linewidth=2.2,
            label=f"{phase_name.upper()} score",
        )

    score_ax.set_ylabel("Structure score")
    score_ax.set_ylim(0.0, 1.02)
    score_ax.set_title(f"Local phase trace from {input_name}")
    score_ax.grid(alpha=0.24, linewidth=0.5)
    score_ax.legend(loc="upper right", ncol=3, frameon=False)

    ql_palette = ["#6d597a", "#e07a5f", "#3d405b", "#81b29a", "#f2cc8f", "#4d908e"]
    for ql_index, order in enumerate(ql_orders):
        color = ql_palette[ql_index % len(ql_palette)]
        ql_ax.plot(steps, mean_ql_values[:, ql_index], color=color, alpha=0.22, linewidth=1.0)
        ql_ax.plot(
            steps,
            smoothed_ql_values[:, ql_index],
            color=color,
            linewidth=2.0,
            label=f"q{order}",
        )

    ql_ax.set_ylabel("Mean q_l")
    ql_ax.set_ylim(bottom=0.0)
    ql_ax.grid(alpha=0.24, linewidth=0.5)
    ql_ax.legend(loc="upper right", ncol=max(1, min(3, len(ql_orders))), frameon=False)

    fraction_labels = ["bcc", "fcc", "hcp", "unknown"]
    for phase_index, phase_name in enumerate(fraction_labels):
        linestyle = "--" if phase_name == "unknown" else "-"
        fraction_ax.plot(
            steps,
            phase_fractions[:, phase_index],
            color=PHASE_COLORS[phase_name],
            linewidth=2.0 if phase_name != "unknown" else 1.5,
            linestyle=linestyle,
            label=f"{phase_name.upper()} fraction" if phase_name != "unknown" else "Unknown fraction",
        )

    for event in events:
        for axis in axes:
            axis.axvline(event.step, color="#111111", linestyle=":", linewidth=1.1, alpha=0.7)
        ql_ax.text(
            event.step,
            ql_ax.get_ylim()[1] * 0.98,
            f"{event.from_phase}->{event.to_phase}",
            rotation=90,
            va="top",
            ha="center",
            fontsize=8,
        )

    fraction_ax.set_xlabel("MD step")
    fraction_ax.set_ylabel("Atom fraction")
    fraction_ax.set_ylim(0.0, 1.05)
    fraction_ax.grid(alpha=0.24, linewidth=0.5)
    fraction_ax.legend(loc="upper right", ncol=2, frameon=False)

    fig.savefig(path, dpi=200)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    input_path = args.input.resolve()
    output_dir = (args.output_dir.resolve() if args.output_dir else default_output_dir(input_path))
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.stride < 1:
        raise ValueError("--stride must be at least 1")
    if args.max_frames is not None and args.max_frames < 1:
        raise ValueError("--max-frames must be positive when provided")
    if args.sigma <= 0.0:
        raise ValueError("--sigma must be positive")
    if args.ql_neighbor_count < 2:
        raise ValueError("--ql-neighbor-count must be at least 2")
    if any(order < 0 for order in args.ql_orders):
        raise ValueError("--ql-orders must be non-negative integers")

    all_frames = read_frames(input_path)
    frames = all_frames[:: args.stride]
    if args.max_frames is not None:
        frames = frames[: args.max_frames]
    if not frames:
        raise ValueError("No frames selected for analysis.")

    templates = build_phase_templates()
    analysis = analyze_frames(
        frames,
        templates=templates,
        sigma=float(args.sigma),
        atom_score_threshold=float(args.atom_score_threshold),
        atom_margin_threshold=float(args.atom_margin_threshold),
        ql_orders=tuple(args.ql_orders),
        ql_neighbor_count=int(args.ql_neighbor_count),
    )

    steps = analysis["steps"]
    mean_scores = analysis["mean_scores"]
    mean_ql_values = analysis["mean_ql_values"]
    phase_fractions = analysis["phase_fractions"]
    smoothed_scores = moving_average(mean_scores, window=int(args.smooth_window))
    smoothed_ql_values = moving_average(mean_ql_values, window=int(args.smooth_window))
    labels = dominant_labels(
        smoothed_scores,
        threshold=float(args.frame_score_threshold),
        margin_threshold=float(args.frame_margin_threshold),
    )
    events = detect_transitions(steps, labels, persistence=int(args.transition_persistence))

    transition_steps = {event.step for event in events}
    csv_path = output_dir / "phase_trace.csv"
    transitions_path = output_dir / "phase_transitions.csv"
    plot_path = output_dir / "phase_trace.png"

    write_phase_trace_csv(
        csv_path,
        steps=steps,
        mean_scores=mean_scores,
        smoothed_scores=smoothed_scores,
        mean_ql_values=mean_ql_values,
        smoothed_ql_values=smoothed_ql_values,
        ql_orders=tuple(args.ql_orders),
        phase_fractions=phase_fractions,
        labels=labels,
        transition_steps=transition_steps,
    )
    write_transition_csv(transitions_path, events)
    plot_phase_trace(
        plot_path,
        input_name=input_path.name,
        steps=steps,
        mean_scores=mean_scores,
        smoothed_scores=smoothed_scores,
        mean_ql_values=mean_ql_values,
        smoothed_ql_values=smoothed_ql_values,
        ql_orders=tuple(args.ql_orders),
        phase_fractions=phase_fractions,
        events=events,
    )

    dominant_counts = {phase: labels.count(phase) for phase in (*PHASES, "unknown")}
    print(f"Analyzed {len(frames)} frames from {input_path}")
    print(f"Outputs written to {output_dir}")
    print(
        "Dominant frames: "
        + ", ".join(f"{phase}={count}" for phase, count in dominant_counts.items() if count > 0)
    )
    if events:
        print("Detected transitions:")
        for event in events:
            print(f"  step {event.step}: {event.from_phase} -> {event.to_phase}")
    else:
        print("Detected transitions: none")


if __name__ == "__main__":
    main()
