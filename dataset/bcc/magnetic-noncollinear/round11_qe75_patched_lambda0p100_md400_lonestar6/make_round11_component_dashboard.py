#!/usr/bin/env python3
"""Render the Round-11 structure and component-magnetization dashboard GIF."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-round11-dashboard")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize


ROOT = Path(__file__).resolve().parent
CASE = "a2p39_m2_lambda0p100_md400"
DEFAULT_NPZ = ROOT / "round11_results_light.npz"
DEFAULT_INPUT = ROOT / "cases" / CASE / "md_400steps.in"
DEFAULT_OUTPUT = (
    ROOT / "round11_qe75_patched_lambda0p100_md400_lonestar6"
    / "cases" / CASE / "md_400steps.out"
)
DEFAULT_PNG = ROOT / "round11_component_magnetization_dashboard_final.png"
DEFAULT_GIF = ROOT / "round11_component_magnetization_dashboard.gif"
NATOMS = 8
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
# Move the visualization boundary away from ideal BCC sites. This is only a
# periodic change of origin; it does not alter distances or the trajectory.
VISUAL_ORIGIN_FRACTIONAL = np.asarray([0.125, 0.125, 0.250], dtype=float)

COLORS = {
    "navy": "#0b2d4d",
    "orange": "#f47c20",
    "purple": "#7b4ab5",
    "green": "#2a9d68",
    "blue": "#2878b5",
    "red": "#c44e52",
    "gray": "#687386",
    "light": "#d8dee7",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--png", type=Path, default=DEFAULT_PNG)
    parser.add_argument("--gif", type=Path, default=DEFAULT_GIF)
    parser.add_argument("--stride", type=int, default=4)
    parser.add_argument("--fps", type=int, default=8)
    parser.add_argument("--dpi", type=int, default=120)
    return parser.parse_args()


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def input_positions_and_cell(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text().splitlines()
    position_start = lines.index("ATOMIC_POSITIONS crystal") + 1
    positions = np.asarray(
        [[number(value) for value in line.split()[1:4]]
         for line in lines[position_start:position_start + NATOMS]],
        dtype=float,
    )
    cell_start = lines.index("CELL_PARAMETERS angstrom") + 1
    cell = np.asarray(
        [[number(value) for value in line.split()[:3]]
         for line in lines[cell_start:cell_start + 3]],
        dtype=float,
    )
    return positions, cell


def output_position_blocks(text: str) -> np.ndarray:
    lines = text.splitlines()
    blocks: list[list[list[float]]] = []
    for index, line in enumerate(lines):
        if not line.startswith("ATOMIC_POSITIONS (crystal)"):
            continue
        rows = []
        for candidate in lines[index + 1:index + 1 + NATOMS]:
            fields = candidate.split()
            rows.append([number(value) for value in fields[1:4]])
        blocks.append(rows)
    return np.asarray(blocks, dtype=float)


def converged_scf_history(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    magnetization: list[list[float]] = []
    accuracy: list[float] = []
    cycles: list[int] = []
    current_magnetization = [np.nan, np.nan, np.nan]
    current_accuracy = np.nan
    cycle_iterations = 0
    for line in text.splitlines():
        match = re.search(
            rf"total magnetization\s*=\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})", line
        )
        if match:
            current_magnetization = [number(value) for value in match.groups()]
        match = re.search(rf"estimated scf accuracy\s*<\s*({FLOAT})\s*Ry", line)
        if match:
            current_accuracy = number(match.group(1))
        if "iteration #" in line:
            cycle_iterations += 1
        if "convergence has been achieved" in line:
            magnetization.append(current_magnetization)
            accuracy.append(current_accuracy)
            cycles.append(cycle_iterations)
            cycle_iterations = 0
    return (
        np.asarray(magnetization, dtype=float),
        np.asarray(accuracy, dtype=float),
        np.asarray(cycles, dtype=int),
    )


def draw_cell(axis, cell: np.ndarray) -> None:
    origin = np.zeros(3)
    a_vector, b_vector, c_vector = cell
    edges = (
        (origin, a_vector), (origin, b_vector), (origin, c_vector),
        (a_vector, a_vector + b_vector), (a_vector, a_vector + c_vector),
        (b_vector, a_vector + b_vector), (b_vector, b_vector + c_vector),
        (c_vector, a_vector + c_vector), (c_vector, b_vector + c_vector),
        (a_vector + b_vector, a_vector + b_vector + c_vector),
        (a_vector + c_vector, a_vector + b_vector + c_vector),
        (b_vector + c_vector, a_vector + b_vector + c_vector),
    )
    for start, end in edges:
        axis.plot(*zip(start, end), color="#64748b", linewidth=1.15, alpha=0.8)


def ideal_bcc_2x2x1() -> np.ndarray:
    return np.asarray(
        [
            [0.00, 0.00, 0.00], [0.25, 0.25, 0.50],
            [0.00, 0.50, 0.00], [0.25, 0.75, 0.50],
            [0.50, 0.00, 0.00], [0.75, 0.25, 0.50],
            [0.50, 0.50, 0.00], [0.75, 0.75, 0.50],
        ],
        dtype=float,
    )


def cartesian_displacement(
    positions: np.ndarray, reference: np.ndarray, cell: np.ndarray
) -> np.ndarray:
    delta = positions - reference
    delta -= np.rint(delta)
    return delta @ cell


def style_axis(axis) -> None:
    axis.grid(color=COLORS["light"], linewidth=0.8, alpha=0.8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=9)


def main() -> int:
    args = arguments()
    if args.stride < 1:
        raise ValueError("--stride must be positive")
    text = args.output.read_text(errors="replace")
    initial_positions, cell = input_positions_and_cell(args.input)
    propagated_positions = output_position_blocks(text)
    magnetization, accuracy, cycles = converged_scf_history(text)

    with np.load(args.npz, allow_pickle=False) as data:
        md_index = int(np.flatnonzero(np.asarray(data["source_kind"]) == "md")[0])
        temperature_all = np.asarray(data["temperature_K"], dtype=float)
        temperature_offsets = np.asarray(data["temperature_K_offset"], dtype=int)
        temperature = temperature_all[
            temperature_offsets[md_index]:temperature_offsets[md_index + 1]
        ]
        pressure_all = np.asarray(data["pressure_kbar"], dtype=float)
        pressure_offsets = np.asarray(data["pressure_kbar_offset"], dtype=int)
        pressure = pressure_all[pressure_offsets[md_index]:pressure_offsets[md_index + 1]] / 10.0
        local_source = np.asarray(data["local_vector_source_index"], dtype=int)
        local_all = np.asarray(data["local_magnetization_vectors_Bohr"], dtype=float)
        local_magnetization = local_all[local_source == md_index]

    # Each converged force/magnetic record precedes its ionic propagation.
    # Pair force record 1 with the input coordinates, then records 2--400 with
    # the coordinates printed after dynamics steps 1--399.
    positions = np.concatenate((initial_positions[None], propagated_positions[:-1]), axis=0)
    local_magnetization = local_magnetization[1:]
    frame_count = len(magnetization)
    arrays = {
        "positions": positions,
        "local magnetization": local_magnetization,
        "temperature": temperature[:frame_count],
        "pressure": pressure[:frame_count],
        "accuracy": accuracy,
        "cycles": cycles,
    }
    for label, values in arrays.items():
        if len(values) != frame_count:
            raise ValueError(f"{label} has {len(values)} records; expected {frame_count}")
    if frame_count != 400:
        raise ValueError(f"Expected 400 converged Round-11 MD SCFs, found {frame_count}")

    configuration = np.arange(frame_count)
    time_ps = configuration * 20.67 * (2.0 * 2.4188843265864e-17) * 1.0e12
    magnetic_magnitude = np.linalg.norm(magnetization, axis=1)
    local_magnitude = np.linalg.norm(local_magnetization, axis=2)
    mean_local_magnitude = np.mean(local_magnitude, axis=1)
    local_direction = local_magnetization / local_magnitude[:, :, None]
    angular_deflection = np.degrees(
        np.arccos(
            np.clip(
                np.sum(local_direction * local_direction[0][None, :, :], axis=2),
                -1.0,
                1.0,
            )
        )
    )
    maximum_angular_deflection = np.max(angular_deflection, axis=1)
    ideal = ideal_bcc_2x2x1()
    displacement_ideal = np.asarray(
        [cartesian_displacement(frame, ideal, cell) for frame in positions]
    )
    displacement_start = np.asarray(
        [cartesian_displacement(frame, positions[0], cell) for frame in positions]
    )
    rms_ideal = np.sqrt(np.mean(np.sum(displacement_ideal**2, axis=2), axis=1))
    maximum_ideal = np.max(np.linalg.norm(displacement_ideal, axis=2), axis=1)
    rms_start = np.sqrt(np.mean(np.sum(displacement_start**2, axis=2), axis=1))

    frame_indices = np.arange(0, frame_count, args.stride, dtype=int)
    if frame_indices[-1] != frame_count - 1:
        frame_indices = np.append(frame_indices, frame_count - 1)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titleweight": "bold",
            "axes.titlecolor": COLORS["navy"],
            "axes.labelcolor": COLORS["navy"],
        }
    )
    figure = plt.figure(figsize=(17, 10), facecolor="white")
    grid = figure.add_gridspec(
        5, 2, width_ratios=(1.55, 1.0), height_ratios=(1.0, 0.85, 0.85, 0.85, 0.85),
        left=0.045, right=0.975, bottom=0.075, top=0.84,
        wspace=0.18, hspace=0.63,
    )
    structure_axis = figure.add_subplot(grid[:, 0], projection="3d")
    magnetic_axis = figure.add_subplot(grid[0, 1])
    temperature_axis = figure.add_subplot(grid[1, 1])
    local_axis = figure.add_subplot(grid[2, 1])
    scf_axis = figure.add_subplot(grid[3, 1])
    displacement_axis = figure.add_subplot(grid[4, 1])

    component_lines = []
    for component, color, label in zip(
        magnetization.T,
        (COLORS["purple"], COLORS["green"], COLORS["blue"]),
        (r"$M_x$", r"$M_y$", r"$M_z$"),
    ):
        component_lines.append(
            magnetic_axis.plot(configuration, component, color=color, linewidth=2.1, label=label)[0]
        )
    magnetic_axis.plot(
        configuration, magnetic_magnitude, color=COLORS["navy"], linewidth=1.6,
        linestyle="--", label=r"$|\mathbf{M}|$",
    )
    magnetic_axis.axhline(0.0, color="#334155", linewidth=0.8)
    magnetic_axis.set_title("Global Magnetization Components", loc="left", fontsize=13)
    magnetic_axis.set_ylabel(r"$\mu_B$/cell", fontsize=10)
    magnetic_axis.legend(frameon=False, ncol=4, fontsize=9, loc="upper center")

    temperature_axis.plot(configuration, temperature[:frame_count], color=COLORS["red"], linewidth=2.1)
    temperature_axis.axhline(4000.0, color=COLORS["red"], linestyle="--", alpha=0.35, label="target")
    temperature_axis.axhline(np.mean(temperature[:frame_count]), color=COLORS["navy"], linestyle=":", alpha=0.7, label="mean")
    temperature_axis.set_title("Ionic Temperature", loc="left", fontsize=13)
    temperature_axis.set_ylabel("K", fontsize=10)
    temperature_axis.legend(frameon=False, fontsize=8.5, ncol=2, loc="upper right")

    local_axis.plot(configuration, mean_local_magnitude, color=COLORS["orange"], linewidth=2.2, label="mean local moment")
    local_axis.axhline(2.0, color=COLORS["gray"], linestyle="--", alpha=0.7, label="2 μB target")
    local_axis.set_title("Local-Moment Retention", loc="left", fontsize=13)
    local_axis.set_ylabel(r"$\mu_B$/Fe", fontsize=10)
    local_axis.legend(frameon=False, fontsize=8.5, ncol=2, loc="upper right")
    angle_axis = local_axis.twinx()
    angle_axis.plot(
        configuration, maximum_angular_deflection, color=COLORS["purple"],
        linewidth=1.4, alpha=0.78,
    )
    angle_axis.set_ylabel("max direction change (°)", fontsize=8.5, color=COLORS["purple"])
    angle_axis.tick_params(axis="y", labelsize=8, colors=COLORS["purple"])

    scf_axis.plot(configuration, cycles, color=COLORS["blue"], linewidth=2.0, label="electronic iterations")
    scf_axis.set_title("SCF Convergence", loc="left", fontsize=13)
    scf_axis.set_ylabel("iterations", fontsize=10)
    scf_accuracy_axis = scf_axis.twinx()
    scf_accuracy_axis.plot(configuration, accuracy, color=COLORS["purple"], linewidth=1.4, alpha=0.75, label="final accuracy")
    scf_accuracy_axis.axhline(4.0e-3, color=COLORS["red"], linestyle="--", alpha=0.35)
    scf_accuracy_axis.set_yscale("log")
    scf_accuracy_axis.set_ylabel("accuracy (Ry)", fontsize=9, color=COLORS["purple"])
    scf_accuracy_axis.tick_params(axis="y", labelsize=8, colors=COLORS["purple"])

    displacement_axis.plot(configuration, rms_ideal, color=COLORS["orange"], linewidth=2.1, label="RMS from ideal BCC")
    displacement_axis.plot(configuration, maximum_ideal, color=COLORS["red"], linewidth=1.7, label="maximum from ideal")
    displacement_axis.plot(configuration, rms_start, color=COLORS["blue"], linewidth=1.9, label="RMS from initial")
    displacement_axis.set_title("Displacement History", loc="left", fontsize=13)
    displacement_axis.set_xlabel("MD configuration", fontsize=10)
    displacement_axis.set_ylabel("Å", fontsize=10)
    displacement_axis.legend(frameon=False, fontsize=8.2, ncol=3, loc="upper center")

    for axis in (magnetic_axis, temperature_axis, local_axis, scf_axis, displacement_axis):
        style_axis(axis)
        axis.set_xlim(-2, frame_count + 1)

    cursors = [
        axis.axvline(0, color=COLORS["orange"], linewidth=1.35)
        for axis in (magnetic_axis, temperature_axis, local_axis, scf_axis, displacement_axis)
    ]
    magnetic_markers = [
        magnetic_axis.plot([], [], marker="o", markersize=5.5, color=line.get_color())[0]
        for line in component_lines
    ]
    temperature_marker = temperature_axis.plot([], [], "o", color=COLORS["red"], markersize=6)[0]
    local_marker = local_axis.plot([], [], "o", color=COLORS["orange"], markersize=6)[0]
    scf_marker = scf_axis.plot([], [], "o", color=COLORS["blue"], markersize=5.5)[0]
    displacement_markers = [
        displacement_axis.plot([], [], "o", color=color, markersize=5.2)[0]
        for color in (COLORS["orange"], COLORS["red"], COLORS["blue"])
    ]

    normalization = Normalize(vmin=-1.0, vmax=1.0)
    color_map = plt.get_cmap("coolwarm")
    color_bar_axis = figure.add_axes([0.11, 0.085, 0.34, 0.022])
    scalar_map = plt.cm.ScalarMappable(norm=normalization, cmap=color_map)
    color_bar = figure.colorbar(scalar_map, cax=color_bar_axis, orientation="horizontal")
    color_bar.set_label(r"Local-spin direction, $m_z/|\mathbf{m}|$", fontsize=10)
    color_bar.ax.tick_params(labelsize=9)

    status = figure.text(
        0.075, 0.82, "", ha="left", va="top", fontsize=12.0, color=COLORS["navy"],
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "white", "edgecolor": COLORS["light"], "alpha": 0.94},
    )
    figure.suptitle(
        "Round 11 Noncollinear BCC Fe AIMD Dashboard",
        x=0.045, y=0.975, ha="left", fontsize=24, fontweight="bold", color=COLORS["navy"],
    )
    figure.text(
        0.047, 0.925,
        "8 unique atoms • periodic 2×2×2 display • true local-spin directions (not exaggerated)",
        fontsize=12.2, color=COLORS["gray"],
    )

    display_cell = cell.copy()
    display_cell[2] *= 2.0
    # Keep scalar limits separate from the array passed to Matplotlib.
    # Axes3D.set_box_aspect() normalizes mutable array-like input in place on
    # some Matplotlib versions; reusing that array made every GIF frame after
    # the standalone PNG use limits near 1 Å and clip most atoms.
    cell_lengths = tuple(float(value) for value in np.linalg.norm(display_cell, axis=1))
    ideal_inside = np.mod(ideal - VISUAL_ORIGIN_FRACTIONAL, 1.0)
    ideal_display = np.concatenate(
        (ideal_inside, ideal_inside + np.asarray([0.0, 0.0, 1.0])), axis=0
    )
    ideal_display_cartesian = ideal_display @ cell
    maximum_local = float(np.max(local_magnitude))

    def update(frame_index: int):
        structure_axis.clear()
        # Use a shifted periodic origin that lies between lattice sites. This
        # keeps every atom inside the box without cutting through boundary
        # atoms or changing the physical periodic structure.
        current_fractional = np.mod(
            positions[frame_index] - VISUAL_ORIGIN_FRACTIONAL, 1.0
        )
        display_fractional = np.concatenate(
            (current_fractional, current_fractional + np.asarray([0.0, 0.0, 1.0])),
            axis=0,
        )
        current_cartesian = display_fractional @ cell
        draw_cell(structure_axis, display_cell)
        structure_axis.scatter(
            ideal_display_cartesian[:, 0], ideal_display_cartesian[:, 1], ideal_display_cartesian[:, 2],
            s=34, facecolor="none", edgecolor="#94a3b8", linewidth=0.8, alpha=0.65,
        )
        current_local = local_magnetization[frame_index]
        current_norm = local_magnitude[frame_index]
        direction_z = np.divide(
            current_local[:, 2], current_norm,
            out=np.zeros_like(current_norm), where=current_norm > 0.0,
        )
        atom_colors = color_map(normalization(direction_z))
        display_colors = np.concatenate((atom_colors, atom_colors.copy()), axis=0)
        display_colors[NATOMS:, 3] = 0.42
        structure_axis.scatter(
            current_cartesian[:, 0], current_cartesian[:, 1], current_cartesian[:, 2],
            s=105, color=display_colors, edgecolor="white", linewidth=0.7, depthshade=True,
        )
        displayed_spins = 0.68 * np.concatenate((current_local, current_local), axis=0) / maximum_local
        structure_axis.quiver(
            current_cartesian[:, 0], current_cartesian[:, 1], current_cartesian[:, 2],
            displayed_spins[:, 0], displayed_spins[:, 1], displayed_spins[:, 2],
            color=display_colors, linewidth=1.8, arrow_length_ratio=0.28,
        )
        vector = magnetization[frame_index]
        vector_norm = magnetic_magnitude[frame_index]
        if vector_norm > 0.0:
            center = 0.5 * np.sum(display_cell, axis=0)
            displayed_global = 0.55 * np.min(cell_lengths) * vector / vector_norm
            structure_axis.quiver(
                *center, *displayed_global, color=COLORS["navy"], linewidth=4.0,
                arrow_length_ratio=0.22,
            )
            structure_axis.text(
                *(center + displayed_global), r" global $\mathbf{M}$",
                color=COLORS["navy"], fontsize=10, fontweight="bold",
            )
        margin = 0.12
        structure_axis.set_xlim(-margin, cell_lengths[0] + margin)
        structure_axis.set_ylim(-margin, cell_lengths[1] + margin)
        structure_axis.set_zlim(-margin, cell_lengths[2] + margin)
        structure_axis.set_box_aspect(tuple(cell_lengths))
        structure_axis.view_init(elev=23, azim=-58)
        structure_axis.set_xlabel("x (Å)", labelpad=9)
        structure_axis.set_ylabel("y (Å)", labelpad=9)
        structure_axis.set_zlabel("z (Å)", labelpad=9)
        structure_axis.set_title(
            f"Periodic BCC Structure and Spin Texture — Configuration {frame_index}",
            loc="left", fontsize=16, pad=18,
        )

        for cursor in cursors:
            cursor.set_xdata([frame_index, frame_index])
        for component_index, marker in enumerate(magnetic_markers):
            marker.set_data([frame_index], [magnetization[frame_index, component_index]])
        temperature_marker.set_data([frame_index], [temperature[frame_index]])
        local_marker.set_data([frame_index], [mean_local_magnitude[frame_index]])
        scf_marker.set_data([frame_index], [cycles[frame_index]])
        for marker, value in zip(
            displacement_markers,
            (rms_ideal[frame_index], maximum_ideal[frame_index], rms_start[frame_index]),
        ):
            marker.set_data([frame_index], [value])

        status.set_text(
            f"Configuration {frame_index}/399   •   t = {time_ps[frame_index]:.4f} ps\n"
            f"T = {temperature[frame_index]:.0f} K   •   P = {pressure[frame_index]:.1f} GPa\n"
            f"Mx = {vector[0]:+.3f}   My = {vector[1]:+.3f}   Mz = {vector[2]:+.3f} μB/cell\n"
            f"|M| = {vector_norm:.3f} μB/cell   •   mean local moment = {mean_local_magnitude[frame_index]:.3f} μB/Fe\n"
            f"Maximum local-spin rotation from configuration 0 = {maximum_angular_deflection[frame_index]:.2f}°\n"
            f"SCF: {cycles[frame_index]} iterations   •   accuracy = {accuracy[frame_index]:.2e} Ry\n"
            f"RMS displacement = {rms_ideal[frame_index]:.3f} Å"
        )
        return (*cursors, *magnetic_markers, temperature_marker, local_marker,
                scf_marker, *displacement_markers, status)

    args.png.parent.mkdir(parents=True, exist_ok=True)
    update(frame_count - 1)
    figure.savefig(args.png, dpi=200, facecolor="white")

    animation = FuncAnimation(
        figure, lambda sampled: update(int(frame_indices[sampled])),
        frames=len(frame_indices), interval=1000 / args.fps, blit=False, repeat=True,
    )
    animation.save(
        args.gif,
        writer=PillowWriter(
            fps=args.fps,
            metadata={"title": "Round 11 component magnetization dashboard", "artist": "IronCoreMD"},
        ),
        dpi=args.dpi,
    )
    plt.close(figure)
    print(f"static_dashboard={args.png}")
    print(f"animated_dashboard={args.gif}")
    print(f"trajectory_configurations={frame_count}")
    print(f"animated_frames={len(frame_indices)} stride={args.stride}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
