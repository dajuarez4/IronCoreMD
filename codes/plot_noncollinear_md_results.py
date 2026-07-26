#!/usr/bin/env python3
"""Plot noncollinear AIMD spin, displacement, and thermodynamic evolution."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BOHR_TO_ANG = 0.529177210903
COLORS = {
    "navy": "#0b2d4d",
    "orange": "#f47c20",
    "blue": "#2878b5",
    "green": "#2a9d68",
    "purple": "#7b4ab5",
    "red": "#c44e52",
    "gray": "#687386",
    "light": "#d8dee7",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("npz", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def metadata(data: np.lib.npyio.NpzFile) -> dict[str, object]:
    raw = data["metadata_json"].item()
    return json.loads(str(raw))


def bcc_ideal_fractional(natoms: int) -> np.ndarray:
    repetitions = round((natoms / 2.0) ** (1.0 / 3.0))
    if 2 * repetitions**3 != natoms:
        raise ValueError(f"Expected a conventional BCC supercell, received {natoms} atoms")
    corners = []
    centers = []
    for z_index in range(repetitions):
        for y_index in range(repetitions):
            for x_index in range(repetitions):
                corners.append(
                    [x_index / repetitions, y_index / repetitions, z_index / repetitions]
                )
                centers.append(
                    [
                        (x_index + 0.5) / repetitions,
                        (y_index + 0.5) / repetitions,
                        (z_index + 0.5) / repetitions,
                    ]
                )
    return np.asarray(corners + centers, dtype=float)


def cell_angstrom(data: np.lib.npyio.NpzFile, meta: dict[str, object]) -> np.ndarray:
    if "input_cell_parameters" in data.files:
        cell = np.asarray(data["input_cell_parameters"], dtype=float)
        unit = str(np.asarray(data["input_cell_unit"]).item()).lower()
        if np.isfinite(cell).all() and unit in {"angstrom", "ang"}:
            return cell
    return (
        np.asarray(data["initial_cell_alat"], dtype=float)
        * float(meta["alat_bohr"])
        * BOHR_TO_ANG
    )


def minimum_image_displacements(
    positions_frac: np.ndarray, reference_frac: np.ndarray, cell: np.ndarray
) -> np.ndarray:
    delta = np.asarray(positions_frac, dtype=float) - np.asarray(reference_frac, dtype=float)
    delta -= np.rint(delta)
    return delta @ cell


def equalize_3d(axis, cell: np.ndarray) -> None:
    extent = float(np.max(np.linalg.norm(cell, axis=1)))
    center = 0.5 * np.sum(cell, axis=0)
    half = 0.57 * extent
    axis.set_xlim(center[0] - half, center[0] + half)
    axis.set_ylim(center[1] - half, center[1] + half)
    axis.set_zlim(center[2] - half, center[2] + half)
    axis.set_box_aspect((1, 1, 1))


def draw_cell(axis, cell: np.ndarray) -> None:
    origin = np.zeros(3)
    a, b, c = cell
    edges = [
        (origin, a), (origin, b), (origin, c), (a, a + b), (a, a + c),
        (b, a + b), (b, b + c), (c, a + c), (c, b + c),
        (a + b, a + b + c), (a + c, a + b + c), (b + c, a + b + c),
    ]
    for start, end in edges:
        axis.plot(*zip(start, end), color="#7b8794", linewidth=0.8, alpha=0.65)


def main() -> int:
    args = parse_args()
    with np.load(args.npz, allow_pickle=False) as data:
        meta = metadata(data)
        positions = np.asarray(data["positions"], dtype=float)
        iterations = np.asarray(data["iteration"], dtype=int)
        valid = np.asarray(data["frame_valid"], dtype=bool)
        positions = positions[valid]
        iterations = iterations[valid]
        temperature = np.asarray(data["temperature_K"], dtype=float)[valid]
        pressure = np.asarray(data["pressure_GPa"], dtype=float)[valid]
        absolute_magnetization = np.asarray(data["abs_mag_total_Bohr"], dtype=float)[valid]
        magnetization = np.asarray(data["mag_total_vector_Bohr"], dtype=float)[valid]
        cell = cell_angstrom(data, meta)

    ideal = bcc_ideal_fractional(positions.shape[1])
    displacement_ideal = minimum_image_displacements(positions, ideal, cell)
    displacement_start = minimum_image_displacements(positions, positions[0], cell)
    rms_ideal = np.sqrt(np.mean(np.sum(displacement_ideal**2, axis=2), axis=1))
    rms_start = np.sqrt(np.mean(np.sum(displacement_start**2, axis=2), axis=1))
    maximum_ideal = np.max(np.linalg.norm(displacement_ideal, axis=2), axis=1)
    magnetization_norm = np.linalg.norm(magnetization, axis=1)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11.5,
            "axes.titleweight": "bold",
            "axes.titlecolor": COLORS["navy"],
            "axes.labelcolor": COLORS["navy"],
        }
    )
    figure = plt.figure(figsize=(15.2, 10.2), facecolor="white")
    grid = figure.add_gridspec(2, 2, hspace=0.32, wspace=0.25)
    spin_axis = figure.add_subplot(grid[0, 0])
    displacement_axis = figure.add_subplot(grid[0, 1])
    structure_axis = figure.add_subplot(grid[1, 0], projection="3d")
    thermo_axis = figure.add_subplot(grid[1, 1])

    for component, color, label in zip(
        magnetization.T,
        (COLORS["purple"], COLORS["green"], COLORS["blue"]),
        (r"$M_x$", r"$M_y$", r"$M_z$"),
    ):
        spin_axis.plot(iterations, component, color=color, linewidth=2.0, label=label)
    spin_axis.plot(
        iterations,
        magnetization_norm,
        color=COLORS["navy"],
        linewidth=2.3,
        linestyle="--",
        label=r"$|\mathbf{M}|$",
    )
    spin_axis.plot(
        iterations,
        absolute_magnetization,
        color=COLORS["orange"],
        linewidth=2.0,
        label=r"$M_\mathrm{abs}$",
    )
    spin_axis.axhline(0.0, color="#334155", linewidth=0.8)
    spin_axis.set_title("Converged Global Magnetization", loc="left")
    spin_axis.set_xlabel("MD iteration")
    spin_axis.set_ylabel(r"Magnetization ($\mu_B$/cell)")
    spin_axis.legend(frameon=False, ncol=3, fontsize=9)

    displacement_axis.plot(
        iterations, rms_ideal, color=COLORS["orange"], linewidth=2.2,
        label="RMS from ideal BCC",
    )
    displacement_axis.plot(
        iterations, maximum_ideal, color=COLORS["red"], linewidth=1.8,
        label="Maximum from ideal BCC",
    )
    displacement_axis.plot(
        iterations, rms_start, color=COLORS["blue"], linewidth=2.0,
        label="RMS change from first MD frame",
    )
    displacement_axis.set_title("Atomic Displacements", loc="left")
    displacement_axis.set_xlabel("MD iteration")
    displacement_axis.set_ylabel("Displacement (Å)")
    displacement_axis.legend(frameon=False, fontsize=9)

    final_frac = np.mod(positions[-1], 1.0)
    final_cart = final_frac @ cell
    ideal_cart = ideal @ cell
    final_displacement = minimum_image_displacements(final_frac, ideal, cell)
    draw_cell(structure_axis, cell)
    structure_axis.scatter(
        final_cart[:, 0], final_cart[:, 1], final_cart[:, 2],
        s=25, color=COLORS["orange"], edgecolor="white", linewidth=0.35,
        depthshade=True,
    )
    structure_axis.quiver(
        ideal_cart[:, 0], ideal_cart[:, 1], ideal_cart[:, 2],
        final_displacement[:, 0], final_displacement[:, 1], final_displacement[:, 2],
        color=COLORS["blue"], linewidth=0.75, arrow_length_ratio=0.20, alpha=0.72,
    )
    equalize_3d(structure_axis, cell)
    structure_axis.view_init(elev=20, azim=-58)
    structure_axis.set_title(
        f"Displacement Field at MD Iteration {iterations[-1]}", loc="left"
    )
    structure_axis.set_xlabel("x (Å)")
    structure_axis.set_ylabel("y (Å)")
    structure_axis.set_zlabel("z (Å)")

    pressure_axis = thermo_axis.twinx()
    thermo_axis.plot(
        iterations, temperature, color=COLORS["red"], linewidth=2.1,
        label="Temperature",
    )
    pressure_axis.plot(
        iterations, pressure, color=COLORS["blue"], linewidth=2.1,
        label="Pressure",
    )
    thermo_axis.axhline(np.mean(temperature), color=COLORS["red"], linestyle="--", alpha=0.45)
    pressure_axis.axhline(np.mean(pressure), color=COLORS["blue"], linestyle="--", alpha=0.45)
    thermo_axis.set_title("Thermodynamic Evolution", loc="left")
    thermo_axis.set_xlabel("MD iteration")
    thermo_axis.set_ylabel("Temperature (K)", color=COLORS["red"])
    pressure_axis.set_ylabel("Pressure (GPa)", color=COLORS["blue"])
    lines = thermo_axis.lines[:1] + pressure_axis.lines[:1]
    thermo_axis.legend(lines, [line.get_label() for line in lines], frameon=False, loc="best")

    for axis in (spin_axis, displacement_axis, thermo_axis):
        axis.grid(color=COLORS["light"], linewidth=0.7, alpha=0.75)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.suptitle(
        "Noncollinear BCC Fe AIMD: Spin and Lattice Evolution",
        x=0.06, y=0.985, ha="left", fontsize=21, fontweight="bold", color=COLORS["navy"],
    )
    figure.text(
        0.06, 0.947,
        (
            f"54 atoms  •  {len(iterations)} complete force frames  •  "
            f"T = {np.mean(temperature):.0f} ± {np.std(temperature):.0f} K  •  "
            f"P = {np.mean(pressure):.1f} ± {np.std(pressure):.1f} GPa"
        ),
        color=COLORS["gray"], fontsize=11.5,
    )
    figure.subplots_adjust(left=0.07, right=0.93, bottom=0.07, top=0.90)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(figure)
    print(f"figure={args.output}")
    print(f"complete_frames={len(iterations)}")
    print(f"mean_temperature_K={np.mean(temperature):.6f}")
    print(f"mean_pressure_GPa={np.mean(pressure):.6f}")
    print(f"final_spin_magnitude_Bohr={magnetization_norm[-1]:.6f}")
    print(f"final_rms_displacement_from_ideal_A={rms_ideal[-1]:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

