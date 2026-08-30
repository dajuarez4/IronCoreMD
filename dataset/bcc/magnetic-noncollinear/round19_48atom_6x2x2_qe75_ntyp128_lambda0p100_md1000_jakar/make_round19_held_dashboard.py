#!/usr/bin/env python3
"""Build a Round-19 HELD trajectory, phonons, still image, and animated dashboard."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-round19-held-dashboard")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
import numpy as np
from scipy.optimize import linear_sum_assignment

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize


HERE = Path(__file__).resolve().parent
IRONCOREMD = HERE.parents[3]
WORKSPACE = IRONCOREMD.parent
HELD_ROOT = WORKSPACE / "HELD"
if str(HELD_ROOT) not in sys.path:
    sys.path.insert(0, str(HELD_ROOT))

from HELD import fit_case, plot_heatmap, plot_mean_dispersion  # noqa: E402


CASE = "a2p39_6x2x2_m2_48species_lambda0p100_md1000"
DEFAULT_INPUT = HERE / "cases" / CASE / "md_1000steps.in"
DEFAULT_OUTPUT = HERE / "md_1000steps.out"
DEFAULT_NPZ = HERE / "round19_held_trajectory.npz"
DEFAULT_COEFFICIENTS = HERE / "round19_held_coefficients.csv"
DEFAULT_STEPS_CSV = HERE / "round19_held_steps.csv"
DEFAULT_DISPERSION_DATA = HERE / "round19_held_mean_dispersion.dat"
DEFAULT_DISPERSION_PNG = HERE / "round19_held_mean_dispersion.png"
DEFAULT_HEATMAP_CACHE = HERE / "round19_held_phonon_frames.npz"
DEFAULT_HEATMAP_PNG = HERE / "round19_held_heatmap.png"
DEFAULT_SUMMARY = HERE / "round19_held_summary.json"
DEFAULT_PNG = HERE / "round19_held_md_dashboard_final.png"
DEFAULT_GIF = HERE / "round19_held_md_dashboard.gif"
NATOMS = 48
REPETITIONS = (6, 2, 2)
REQUESTED_STEPS = 1000
REQUESTED_MOMENT = 2.0
DT_AU = 20.67
AU_TIME_SECONDS = 2.4188843265864e-17
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"

COLORS = {
    "navy": "#0b2d4d", "orange": "#f47c20", "purple": "#7b4ab5",
    "green": "#2a9d68", "blue": "#2878b5", "red": "#c44e52",
    "gray": "#687386", "light": "#d8dee7",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stride", type=int, default=1,
                        help="Keep every Nth frame; the final frame is always included")
    parser.add_argument("--fps", type=int, default=5)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--points-per-segment", type=int, default=60)
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ)
    parser.add_argument("--gif", type=Path, default=DEFAULT_GIF)
    parser.add_argument("--png", type=Path, default=DEFAULT_PNG)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args()


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def scalars(pattern: str, text: str, flags: int = re.I) -> np.ndarray:
    return np.asarray([number(value) for value in re.findall(pattern, text, flags)], dtype=float)


def input_positions_cell_symbols(path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    lines = path.read_text().splitlines()
    p0 = lines.index("ATOMIC_POSITIONS crystal") + 1
    position_rows = [row.split() for row in lines[p0:p0 + NATOMS]]
    positions = np.asarray([[number(value) for value in row[1:4]] for row in position_rows])
    symbols = ["Fe"] * NATOMS
    c0 = lines.index("CELL_PARAMETERS angstrom") + 1
    cell = np.asarray([[number(value) for value in row.split()[:3]] for row in lines[c0:c0 + 3]])
    return positions, cell, symbols


def output_position_blocks(text: str) -> np.ndarray:
    lines = text.splitlines()
    blocks: list[list[list[float]]] = []
    for index, line in enumerate(lines):
        if not line.startswith("ATOMIC_POSITIONS (crystal)"):
            continue
        rows: list[list[float]] = []
        for candidate in lines[index + 1:index + 1 + NATOMS]:
            fields = candidate.split()
            if len(fields) < 4:
                break
            rows.append([number(value) for value in fields[1:4]])
        if len(rows) == NATOMS:
            blocks.append(rows)
    return np.asarray(blocks, dtype=float).reshape(-1, NATOMS, 3)


def force_blocks(text: str) -> np.ndarray:
    """Read only the total force block following each QE force header."""
    lines = text.splitlines()
    blocks: list[list[list[float]]] = []
    force_pattern = re.compile(
        rf"atom\s+\d+\s+type\s+\d+\s+force\s*=\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        re.I,
    )
    for index, line in enumerate(lines):
        if "Forces acting on atoms (cartesian axes, Ry/au)" not in line:
            continue
        rows: list[list[float]] = []
        cursor = index + 1
        while cursor < len(lines) and len(rows) < NATOMS:
            match = force_pattern.search(lines[cursor])
            if match:
                rows.append([number(value) for value in match.groups()])
            elif rows and lines[cursor].strip():
                break
            cursor += 1
        if len(rows) == NATOMS:
            blocks.append(rows)
    return np.asarray(blocks, dtype=float).reshape(-1, NATOMS, 3)


def vector_frames(text: str, label: str) -> np.ndarray:
    rows = re.findall(
        rf"^\s*{re.escape(label)}\s*:\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        text, re.I | re.M,
    )
    values = np.asarray([[number(value) for value in row] for row in rows], dtype=float)
    complete = len(values) // NATOMS
    return values[:complete * NATOMS].reshape(complete, NATOMS, 3)


def converged_cell_history(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    vectors: list[list[float]] = []
    absolute: list[float] = []
    iterations: list[int] = []
    current_vector = [np.nan, np.nan, np.nan]
    current_absolute = np.nan
    count = 0
    for line in text.splitlines():
        match = re.search(
            rf"total magnetization\s*=\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})", line
        )
        if match:
            current_vector = [number(value) for value in match.groups()]
        match = re.search(rf"absolute magnetization\s*=\s*({FLOAT})", line)
        if match:
            current_absolute = number(match.group(1))
        if "iteration #" in line:
            count += 1
        if "convergence has been achieved" in line:
            vectors.append(current_vector)
            absolute.append(current_absolute)
            iterations.append(count)
            count = 0
    return np.asarray(vectors), np.asarray(absolute), np.asarray(iterations)


def canonical_bcc_fractional(repetitions: tuple[int, int, int]) -> np.ndarray:
    nx, ny, nz = repetitions
    sites = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                sites.append([ix / nx, iy / ny, iz / nz])
                sites.append([(ix + 0.5) / nx, (iy + 0.5) / ny, (iz + 0.5) / nz])
    return np.asarray(sites, dtype=float)


def ideal_sites_in_atom_order(initial: np.ndarray, cell: np.ndarray) -> np.ndarray:
    canonical = canonical_bcc_fractional(REPETITIONS)
    delta = initial[:, None, :] - canonical[None, :, :]
    delta -= np.rint(delta)
    squared_distance = np.sum((delta @ cell) ** 2, axis=2)
    rows, columns = linear_sum_assignment(squared_distance)
    if not np.array_equal(rows, np.arange(NATOMS)):
        raise ValueError("Unexpected atom rows from ideal-site assignment")
    distances = np.sqrt(squared_distance[rows, columns])
    if float(np.max(distances)) > 1.0:
        raise ValueError(f"Ideal-site assignment is unsafe; max distance={np.max(distances):.3f} Å")
    return canonical[columns]


def selected_indices(count: int, stride: int) -> np.ndarray:
    if stride < 1:
        raise ValueError("--stride must be positive")
    indices = np.arange(0, count, stride, dtype=int)
    if indices[-1] != count - 1:
        indices = np.append(indices, count - 1)
    return indices


def prepare_trajectory(args: argparse.Namespace) -> dict[str, np.ndarray]:
    text = args.output.read_text(errors="replace")
    initial, cell, symbols = input_positions_cell_symbols(args.input)
    propagated = output_position_blocks(text)
    forces = force_blocks(text)
    total_magnetization, absolute_magnetization, scf_iterations = converged_cell_history(text)
    count = min(len(propagated), len(forces), len(total_magnetization))
    if count < 2:
        raise ValueError(f"Need at least two complete QE frames, found {count}")

    # QE evaluates force record 1 at the input positions. Each printed position
    # block is the configuration propagated after that force evaluation.
    force_aligned_positions = np.concatenate((initial[None], propagated[:-1]), axis=0)[:count]
    local_magnetization = vector_frames(text, "magnetization")
    constrained_moments = vector_frames(text, "constrained moment")
    if len(local_magnetization) >= count + 1:
        local_magnetization = local_magnetization[1:count + 1]
    else:
        local_magnetization = local_magnetization[:count]
    if len(constrained_moments) >= count + 1:
        constrained_moments = constrained_moments[1:count + 1]
    else:
        constrained_moments = constrained_moments[:count]

    temperature = scalars(
        rf"^\s*temperature\s*=\s*({FLOAT})\s*K", text, re.I | re.M
    )[:count]
    pressure = scalars(rf"(?m)(?:^|\s)P\s*=\s*({FLOAT})", text)[:count] / 10.0
    if min(len(local_magnetization), len(constrained_moments), len(temperature), len(pressure)) < count:
        raise ValueError("One or more QE observables do not cover every complete force frame")

    ideal = ideal_sites_in_atom_order(initial, cell)
    dt_ps = DT_AU * 2.0 * AU_TIME_SECONDS * 1.0e12
    all_iterations = np.arange(1, count + 1, dtype=int)
    all_time = np.arange(count, dtype=float) * dt_ps
    keep = selected_indices(count, args.stride)
    payload: dict[str, np.ndarray] = {
        "input_cell_parameters": np.asarray(cell, dtype=np.float64),
        "input_cell_unit": np.asarray("angstrom", dtype="U16"),
        "initial_positions_alat": np.asarray(ideal, dtype=np.float64),
        "initial_cell_alat": np.eye(3, dtype=np.float64),
        "positions": np.asarray(force_aligned_positions[keep], dtype=np.float64),
        "positions_unit": np.full(len(keep), "crystal", dtype="U16"),
        "forces_ry_au": np.asarray(forces[:count][keep], dtype=np.float64),
        "symbols": np.asarray(symbols, dtype="U2"),
        "species": np.asarray(symbols, dtype="U2"),
        "iteration": all_iterations[keep],
        "time_ps": all_time[keep],
        "temperature_K": temperature[keep],
        "pressure_GPa": pressure[keep],
        "mag_total_vector_Bohr": total_magnetization[:count][keep],
        "mag_total_Bohr": np.linalg.norm(total_magnetization[:count][keep], axis=1),
        "abs_mag_total_Bohr": absolute_magnetization[:count][keep],
        "local_magnetization_Bohr": local_magnetization[keep],
        "printed_constrained_moments_Bohr": constrained_moments[keep],
        "scf_iterations": scf_iterations[:count][keep],
        "source_complete_frames": np.asarray([count], dtype=np.int32),
        "source_requested_frames": np.asarray([REQUESTED_STEPS], dtype=np.int32),
        "source_job_done": np.asarray(["JOB DONE" in text], dtype=bool),
    }
    args.npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.npz, **payload)
    return payload


def run_held(args: argparse.Namespace) -> tuple[dict[str, object], dict[str, object]]:
    result, info = fit_case(
        phase="bcc", npz_path=args.npz, output_csv=DEFAULT_COEFFICIENTS,
        output_steps_csv=DEFAULT_STEPS_CSV, aggregate="mean", skip=0, every=1,
        verbose=True,
    )
    dispersion = plot_mean_dispersion(
        phase="bcc", npz_path=args.npz, held_csv=DEFAULT_COEFFICIENTS,
        output_data=DEFAULT_DISPERSION_DATA, output_plot=DEFAULT_DISPERSION_PNG,
        points_per_segment=args.points_per_segment,
    )
    heatmap = plot_heatmap(
        phase="bcc", npz_path=args.npz, held_csv=DEFAULT_COEFFICIENTS,
        cache_npz=DEFAULT_HEATMAP_CACHE, output_plot=DEFAULT_HEATMAP_PNG,
        points_per_segment=args.points_per_segment, y_bins=600,
        force_recompute=True, verbose=True,
    )
    frequencies = np.asarray(heatmap["step_frequencies_thz"], dtype=float)
    with np.load(args.npz, allow_pickle=False) as data:
        complete_frames = int(np.asarray(data["source_complete_frames"]).item())
        requested_frames = int(np.asarray(data["source_requested_frames"]).item())
        job_done = bool(np.asarray(data["source_job_done"]).item())
        constraint_mean = float(np.mean(np.linalg.norm(
            np.asarray(data["printed_constrained_moments_Bohr"], dtype=float), axis=2
        )))
    limitations = []
    if not job_done or complete_frames != requested_frames:
        limitations.append(
            f"The local source contains only {complete_frames}/{requested_frames} complete MD "
            "force frames and no valid complete-run marker."
        )
    if not np.isclose(constraint_mean, REQUESTED_MOMENT, rtol=0.01, atol=0.01):
        limitations.append(
            f"The run printed a {constraint_mean:.3f} Bohr-magneton constraint rather than "
            f"the requested {REQUESTED_MOMENT:.1f}."
        )
    limitations.extend([
        "The MD forces use conv_thr=4e-3 Ry and have substantial SCF-force corrections.",
        "The 6x2x2 cell is shorter than twice the default fifth-shell cutoff in y and z.",
        "The animated HELD branches are qualitative per-frame diagnostics, not publication-quality phonons.",
    ])
    summary: dict[str, object] = {
        "phase": "bcc", "method": "HELD per-frame diagnostic", "natoms": NATOMS,
        "supercell": list(REPETITIONS), "selected_frames": len(result.step_ids),
        "first_iteration": int(result.step_ids[0]), "last_iteration": int(result.step_ids[-1]),
        "stride": args.stride, "num_shells": int(info["num_shells"]),
        "shell_distances_ang": [float(value) for value in info["selected_shell_distances"]],
        "frequency_min_THz": float(np.min(frequencies)),
        "frequency_max_THz": float(np.max(frequencies)),
        "imaginary_values": int(np.count_nonzero(frequencies < -1.0e-8)),
        "source_complete_frames": complete_frames,
        "source_requested_frames": requested_frames,
        "source_job_done": job_done,
        "mean_printed_constraint_Bohr": constraint_mean,
        "limitations": limitations,
        "outputs": {
            "trajectory_npz": str(args.npz), "coefficients_csv": str(DEFAULT_COEFFICIENTS),
            "steps_csv": str(DEFAULT_STEPS_CSV), "dispersion_data": str(DEFAULT_DISPERSION_DATA),
            "dispersion_png": str(DEFAULT_DISPERSION_PNG), "phonon_cache": str(DEFAULT_HEATMAP_CACHE),
            "heatmap_png": str(DEFAULT_HEATMAP_PNG), "dashboard_png": str(args.png),
            "dashboard_gif": str(args.gif),
        },
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return dispersion, heatmap


def displacement_cartesian(positions: np.ndarray, reference: np.ndarray, cell: np.ndarray) -> np.ndarray:
    delta = positions - reference
    delta -= np.rint(delta)
    return delta @ cell


def draw_cell(axis, cell: np.ndarray) -> None:
    origin = np.zeros(3)
    a, b, c = cell
    edges = (
        (origin, a), (origin, b), (origin, c), (a, a + b), (a, a + c),
        (b, a + b), (b, b + c), (c, a + c), (c, b + c),
        (a + b, a + b + c), (a + c, a + b + c), (b + c, a + b + c),
    )
    for start, end in edges:
        axis.plot(*zip(start, end), color="#64748b", linewidth=1.05, alpha=0.75)


def style_history_axis(axis) -> None:
    axis.grid(color=COLORS["light"], linewidth=0.8, alpha=0.8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=10)


def render_dashboard(args: argparse.Namespace, heatmap: dict[str, object]) -> None:
    with np.load(args.npz, allow_pickle=False) as data:
        positions = np.asarray(data["positions"], dtype=float)
        forces = np.asarray(data["forces_ry_au"], dtype=float)
        iterations = np.asarray(data["iteration"], dtype=int)
        time_ps = np.asarray(data["time_ps"], dtype=float)
        temperature = np.asarray(data["temperature_K"], dtype=float)
        pressure = np.asarray(data["pressure_GPa"], dtype=float)
        magnetization = np.asarray(data["mag_total_vector_Bohr"], dtype=float)
        local_magnetization = np.asarray(data["local_magnetization_Bohr"], dtype=float)
        constrained = np.asarray(data["printed_constrained_moments_Bohr"], dtype=float)
        absolute_magnetization = np.asarray(data["abs_mag_total_Bohr"], dtype=float)
        ideal = np.asarray(data["initial_positions_alat"], dtype=float)
        cell = np.asarray(data["input_cell_parameters"], dtype=float)
        source_complete_frames = int(np.asarray(data["source_complete_frames"]).item())
        source_requested_frames = int(np.asarray(data["source_requested_frames"]).item())
        source_job_done = bool(np.asarray(data["source_job_done"]).item())
    phonon_path = np.asarray(heatmap["x_values"], dtype=float)
    phonons = np.asarray(heatmap["step_frequencies_thz"], dtype=float)
    if len(phonons) != len(positions):
        raise ValueError("HELD phonon frames do not match selected trajectory frames")

    ideal_cartesian = ideal @ cell
    displacement_ideal = np.asarray([displacement_cartesian(frame, ideal, cell) for frame in positions])
    displacement_start = np.asarray([displacement_cartesian(frame, positions[0], cell) for frame in positions])
    rms_ideal = np.sqrt(np.mean(np.sum(displacement_ideal**2, axis=2), axis=1))
    maximum_ideal = np.max(np.linalg.norm(displacement_ideal, axis=2), axis=1)
    rms_start = np.sqrt(np.mean(np.sum(displacement_start**2, axis=2), axis=1))
    magnetic_magnitude = np.linalg.norm(magnetization, axis=1)
    local_magnitude = np.linalg.norm(local_magnetization, axis=2)
    constrained_magnitude = np.mean(np.linalg.norm(constrained, axis=2), axis=1)
    maximum_local = max(float(np.max(local_magnitude)), 1.0e-12)
    force_rms = np.sqrt(np.mean(np.sum(forces**2, axis=2), axis=1))

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "axes.titleweight": "bold",
        "axes.titlecolor": COLORS["navy"], "axes.labelcolor": COLORS["navy"],
    })
    figure = plt.figure(figsize=(17, 10), facecolor="white")
    grid = figure.add_gridspec(
        4, 2, width_ratios=(1.55, 1.0), height_ratios=(0.9, 0.8, 0.8, 1.15),
        left=0.045, right=0.975, bottom=0.075, top=0.84, wspace=0.18, hspace=0.48,
    )
    structure_axis = figure.add_subplot(grid[:, 0], projection="3d")
    magnetic_axis = figure.add_subplot(grid[0, 1])
    temperature_axis = figure.add_subplot(grid[1, 1])
    displacement_axis = figure.add_subplot(grid[2, 1])
    phonon_axis = figure.add_subplot(grid[3, 1])

    component_lines = []
    for component, color, label in zip(
        magnetization.T, (COLORS["purple"], COLORS["green"], COLORS["blue"]),
        (r"$M_x$", r"$M_y$", r"$M_z$"),
    ):
        component_lines.append(magnetic_axis.plot(iterations, component, color=color, lw=2.3, label=label)[0])
    magnetic_axis.plot(iterations, magnetic_magnitude, color=COLORS["navy"], lw=1.8,
                       ls="--", label=r"$|\mathbf{M}|$")
    magnetic_axis.axhline(0.0, color="#334155", lw=0.9)
    magnetic_axis.set_title("Global Magnetization Components", loc="left", fontsize=14)
    magnetic_axis.set_ylabel(r"Magnetization ($\mu_B$/cell)", fontsize=11)
    magnetic_axis.legend(frameon=False, ncol=4, fontsize=10, loc="upper center")

    temperature_axis.plot(iterations, temperature, color=COLORS["red"], lw=2.5)
    temperature_axis.axhline(4000.0, color=COLORS["gray"], ls="--", alpha=0.65)
    temperature_axis.set_title("Ionic Temperature", loc="left", fontsize=14)
    temperature_axis.set_ylabel("Temperature (K)", fontsize=11)

    displacement_axis.plot(iterations, rms_ideal, color=COLORS["orange"], lw=2.5,
                           label="RMS from ideal BCC")
    displacement_axis.plot(iterations, maximum_ideal, color=COLORS["red"], lw=2.0,
                           label="Maximum from ideal BCC")
    displacement_axis.plot(iterations, rms_start, color=COLORS["blue"], lw=2.2,
                           label="RMS motion from first frame")
    displacement_axis.set_title("Displacement History", loc="left", fontsize=14)
    displacement_axis.set_xlabel("MD iteration", fontsize=11)
    displacement_axis.set_ylabel("Displacement (Å)", fontsize=11)
    displacement_axis.legend(frameon=False, fontsize=9.2, loc="upper center")

    mean_phonons = np.mean(phonons, axis=0)
    tick_positions = np.asarray(heatmap["tick_positions"], dtype=float)
    tick_labels = list(heatmap["tick_labels"])
    for position in tick_positions:
        phonon_axis.axvline(position, color="#a8a8a8", lw=0.7, alpha=0.65)
    phonon_axis.axhline(0.0, color="#334155", lw=0.8)
    for branch in range(mean_phonons.shape[1]):
        phonon_axis.plot(
            phonon_path, mean_phonons[:, branch], color=COLORS["navy"], lw=1.0,
            ls="--", alpha=0.30, label=f"{len(phonons)}-frame mean" if branch == 0 else None,
        )
    phonon_lines = [phonon_axis.plot([], [], color=COLORS["orange"], lw=2.0)[0]
                    for _ in range(phonons.shape[2])]
    lower = float(np.min(phonons)); upper = float(np.max(phonons))
    padding = 0.07 * max(upper - lower, 1.0)
    phonon_axis.set(xlim=(float(phonon_path[0]), float(phonon_path[-1])),
                    ylim=(lower - padding, upper + padding))
    phonon_axis.set_xticks(tick_positions); phonon_axis.set_xticklabels(tick_labels)
    phonon_axis.set_title("HELD Step Phonons — qualitative", loc="left", fontsize=14)
    phonon_axis.set_ylabel("Frequency (THz)", fontsize=11)
    phonon_axis.legend(frameon=False, fontsize=9, loc="upper left")

    for axis in (magnetic_axis, temperature_axis, displacement_axis, phonon_axis):
        style_history_axis(axis)
    for axis in (magnetic_axis, temperature_axis, displacement_axis):
        axis.set_xlim(iterations[0] - 0.5, iterations[-1] + 0.5)

    cursors = [axis.axvline(iterations[0], color=COLORS["orange"], lw=1.4)
               for axis in (magnetic_axis, temperature_axis, displacement_axis)]
    magnetic_markers = [magnetic_axis.plot([], [], "o", ms=6, color=line.get_color())[0]
                        for line in component_lines]
    temperature_marker = temperature_axis.plot([], [], "o", color=COLORS["red"], ms=7)[0]
    displacement_markers = [displacement_axis.plot([], [], "o", color=color, ms=6)[0]
                            for color in (COLORS["orange"], COLORS["red"], COLORS["blue"])]

    normalization = Normalize(vmin=-1.0, vmax=1.0)
    color_map = plt.get_cmap("coolwarm")
    color_bar_axis = figure.add_axes([0.11, 0.085, 0.34, 0.022])
    scalar_map = plt.cm.ScalarMappable(norm=normalization, cmap=color_map)
    color_bar = figure.colorbar(scalar_map, cax=color_bar_axis, orientation="horizontal")
    color_bar.set_label(r"Evolved local-spin direction, $m_z/|\mathbf{m}|$", fontsize=10)
    color_bar.ax.tick_params(labelsize=9)

    status = figure.text(
        0.075, 0.82, "", ha="left", va="top", fontsize=11.5, color=COLORS["navy"],
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "white",
              "edgecolor": COLORS["light"], "alpha": 0.94},
    )
    figure.suptitle("Round 19 — Noncollinear BCC Fe HELD Dashboard", x=0.045, y=0.975,
                    ha="left", fontsize=24, fontweight="bold", color=COLORS["navy"])
    constraint_average = float(np.mean(constrained_magnitude))
    complete_label = "complete" if source_job_done and source_complete_frames == source_requested_frames else "partial"
    constraint_label = (
        f"constraint={constraint_average:.3f} μB (requested {REQUESTED_MOMENT:.1f} μB)"
        if not np.isclose(constraint_average, REQUESTED_MOMENT, rtol=0.01, atol=0.01)
        else f"constraint={constraint_average:.3f} μB"
    )
    figure.text(
        0.047, 0.925,
        f"48 atoms (6×2×2) • {complete_label} {source_complete_frames}/{source_requested_frames}-step "
        f"source • {len(iterations)} displayed frames • HELD diagnostic phonons • {constraint_label}",
        fontsize=11.7, color=COLORS["red"] if complete_label == "partial" else COLORS["gray"],
        fontweight="bold",
    )

    lengths = np.linalg.norm(cell, axis=1)
    spin_scale = 0.55 * 2.39

    def update(frame_index: int):
        structure_axis.clear()
        current_fractional = np.mod(positions[frame_index], 1.0)
        current_cartesian = current_fractional @ cell
        draw_cell(structure_axis, cell)
        structure_axis.scatter(
            ideal_cartesian[:, 0], ideal_cartesian[:, 1], ideal_cartesian[:, 2],
            s=22, facecolor="none", edgecolor="#94a3b8", lw=0.7, alpha=0.55,
        )
        structure_axis.scatter(
            current_cartesian[:, 0], current_cartesian[:, 1], current_cartesian[:, 2],
            s=62, color=COLORS["orange"], edgecolor="white", lw=0.6, depthshade=True,
        )
        local = local_magnetization[frame_index]
        norms = local_magnitude[frame_index]
        direction_z = np.divide(local[:, 2], norms, out=np.zeros_like(norms), where=norms > 0)
        displayed = spin_scale * local / maximum_local
        structure_axis.quiver(
            current_cartesian[:, 0], current_cartesian[:, 1], current_cartesian[:, 2],
            displayed[:, 0], displayed[:, 1], displayed[:, 2],
            color=color_map(normalization(direction_z)), lw=1.35, arrow_length_ratio=0.30, alpha=0.95,
        )
        vector = magnetization[frame_index]
        vector_norm = magnetic_magnitude[frame_index]
        if vector_norm > 0:
            center = 0.5 * np.sum(cell, axis=0)
            shown_global = 0.9 * 2.39 * vector / vector_norm
            structure_axis.quiver(*center, *shown_global, color=COLORS["navy"], lw=4.0,
                                  arrow_length_ratio=0.22)
            structure_axis.text(*(center + shown_global), r" global $\mathbf{M}$",
                                color=COLORS["navy"], fontsize=10, fontweight="bold")
        # Set the non-cubic box aspect before the numerical data limits. Some
        # Matplotlib 3D backends otherwise replace the limits with normalized
        # aspect values on repeated animation updates.
        structure_axis.set_box_aspect((3.0, 1.0, 1.0))
        margin = 0.12
        structure_axis.set_xlim(-margin, lengths[0] + margin)
        structure_axis.set_ylim(-margin, lengths[1] + margin)
        structure_axis.set_zlim(-margin, lengths[2] + margin)
        structure_axis.view_init(elev=22, azim=-61)
        structure_axis.set_xlabel("x (Å)", labelpad=9)
        structure_axis.set_ylabel("y (Å)", labelpad=9)
        structure_axis.set_zlabel("z (Å)", labelpad=9)
        structure_axis.set_title(
            f"Atomic Motion and Spin Texture — MD Iteration {iterations[frame_index]}",
            loc="left", fontsize=15, pad=18,
        )

        iteration = iterations[frame_index]
        for cursor in cursors:
            cursor.set_xdata([iteration, iteration])
        for component_index, marker in enumerate(magnetic_markers):
            marker.set_data([iteration], [magnetization[frame_index, component_index]])
        temperature_marker.set_data([iteration], [temperature[frame_index]])
        for marker, value in zip(displacement_markers,
                                 (rms_ideal[frame_index], maximum_ideal[frame_index], rms_start[frame_index])):
            marker.set_data([iteration], [value])
        for branch, line in enumerate(phonon_lines):
            line.set_data(phonon_path, phonons[frame_index, :, branch])

        status.set_text(
            f"Frame {frame_index + 1}/{len(iterations)}   •   t = {time_ps[frame_index]:.3f} ps\n"
            f"T = {temperature[frame_index]:.0f} K   •   P = {pressure[frame_index]:.1f} GPa\n"
            f"M = ({vector[0]:+.3f}, {vector[1]:+.3f}, {vector[2]:+.3f}) μB\n"
            f"|M| = {vector_norm:.3f} μB   •   Mabs = {absolute_magnetization[frame_index]:.2f} μB\n"
            f"RMS displacement = {rms_ideal[frame_index]:.3f} Å   •   force RMS = {force_rms[frame_index]:.3f} Ry/Bohr\n"
            f"Local moment mean = {np.mean(norms):.3f} μB   •   printed constraint = {constrained_magnitude[frame_index]:.3f} μB"
        )
        return (*cursors, *magnetic_markers, temperature_marker, *displacement_markers,
                *phonon_lines, status)

    args.png.parent.mkdir(parents=True, exist_ok=True)
    update(len(iterations) - 1)
    figure.savefig(args.png, dpi=220, facecolor="white")
    animation = FuncAnimation(figure, update, frames=len(iterations),
                              interval=1000 / args.fps, blit=False, repeat=True)
    animation.save(
        args.gif,
        writer=PillowWriter(fps=args.fps, metadata={
            "title": "Round 19 noncollinear BCC Fe HELD dashboard",
            "artist": "IronCoreMD / HELD",
        }),
        dpi=args.dpi,
    )
    plt.close(figure)


def main() -> int:
    args = arguments()
    if not HELD_ROOT.is_dir():
        raise FileNotFoundError(f"HELD checkout not found: {HELD_ROOT}")
    payload = prepare_trajectory(args)
    _dispersion, heatmap = run_held(args)
    render_dashboard(args, heatmap)
    print(f"trajectory={args.npz.resolve()}")
    print(f"animated_dashboard={args.gif.resolve()}")
    print(f"static_dashboard={args.png.resolve()}")
    print(f"summary={args.summary.resolve()}")
    print(f"frames={len(payload['iteration'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
