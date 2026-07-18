#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-ironcore")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


FLOAT_PATTERN = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"


def qe_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def indexed_values(text: str, keyword: str) -> dict[int, float]:
    matches = re.findall(
        rf"{keyword}\s*\(\s*(\d+)\s*\)\s*=\s*({FLOAT_PATTERN})",
        text,
        flags=re.IGNORECASE,
    )
    return {int(index): qe_float(value) for index, value in matches}


def parse_three_column_block(
    lines: list[str], header: str, count: int
) -> tuple[list[str], np.ndarray]:
    header_index = next(
        index
        for index, line in enumerate(lines)
        if line.strip().upper().startswith(header)
    )
    labels = []
    values = []
    for line in lines[header_index + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 4:
            break
        try:
            vector = [qe_float(parts[column]) for column in range(1, 4)]
        except ValueError:
            break
        labels.append(parts[0])
        values.append(vector)
        if len(values) == count:
            break
    if len(values) != count:
        raise ValueError(f"Expected {count} rows after {header}, found {len(values)}.")
    return labels, np.asarray(values, dtype=float)


def parse_cell(lines: list[str]) -> np.ndarray:
    header_index = next(
        index
        for index, line in enumerate(lines)
        if line.strip().upper().startswith("CELL_PARAMETERS")
    )
    cell = []
    for line in lines[header_index + 1 : header_index + 4]:
        values = re.findall(FLOAT_PATTERN, line)
        if len(values) < 3:
            raise ValueError("Incomplete CELL_PARAMETERS block.")
        cell.append([qe_float(value) for value in values[:3]])
    return np.asarray(cell, dtype=float)


def draw_cell(axis, cell: np.ndarray) -> None:
    corners = np.asarray(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ],
        dtype=float,
    ) @ cell
    edges = (
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 4),
        (1, 5),
        (2, 4),
        (2, 6),
        (3, 5),
        (3, 6),
        (4, 7),
        (5, 7),
        (6, 7),
    )
    for start, end in edges:
        axis.plot(
            *zip(corners[start], corners[end]),
            color="#5d6775",
            linewidth=0.9,
            alpha=0.58,
        )


def set_equal_3d_axes(axis, points: np.ndarray) -> None:
    minimum = np.min(points, axis=0)
    maximum = np.max(points, axis=0)
    center = 0.5 * (minimum + maximum)
    radius = 0.55 * float(np.max(maximum - minimum))
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_box_aspect((1, 1, 1))


def write_spin_csv(
    path: Path,
    labels: list[str],
    fractional_positions: np.ndarray,
    cartesian_positions: np.ndarray,
    magnitudes: np.ndarray,
    theta_deg: np.ndarray,
    phi_deg: np.ndarray,
    spins: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "atom",
                "label",
                "frac_x",
                "frac_y",
                "frac_z",
                "x_ang",
                "y_ang",
                "z_ang",
                "starting_magnetization",
                "theta_deg",
                "phi_deg",
                "spin_x",
                "spin_y",
                "spin_z",
            ]
        )
        for atom_index in range(len(labels)):
            writer.writerow(
                [
                    atom_index + 1,
                    labels[atom_index],
                    *fractional_positions[atom_index],
                    *cartesian_positions[atom_index],
                    magnitudes[atom_index],
                    theta_deg[atom_index],
                    phi_deg[atom_index],
                    *spins[atom_index],
                ]
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot atom-resolved noncollinear QE starting-spin directions."
    )
    parser.add_argument("input_file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--arrow-length", type=float, default=0.72)
    parser.add_argument("--elev", type=float, default=22.0)
    parser.add_argument("--azim", type=float, default=-55.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    text = args.input_file.read_text(errors="replace")
    lines = text.splitlines()
    magnitudes_by_type = indexed_values(text, "starting_magnetization")
    theta_by_type = indexed_values(text, "angle1")
    phi_by_type = indexed_values(text, "angle2")
    if not magnitudes_by_type:
        raise ValueError("No starting_magnetization values found.")
    type_count = len(magnitudes_by_type)
    if len(theta_by_type) != type_count or len(phi_by_type) != type_count:
        raise ValueError("Incomplete angle1/angle2 spin specification.")

    labels, fractional_positions = parse_three_column_block(
        lines, "ATOMIC_POSITIONS", type_count
    )
    cell = parse_cell(lines)
    cartesian_positions = fractional_positions @ cell

    type_indices = np.asarray(
        [int(re.search(r"(\d+)$", label).group(1)) for label in labels], dtype=int
    )
    magnitudes = np.asarray(
        [magnitudes_by_type[index] for index in type_indices], dtype=float
    )
    theta_deg = np.asarray([theta_by_type[index] for index in type_indices], dtype=float)
    phi_deg = np.asarray([phi_by_type[index] for index in type_indices], dtype=float)
    theta = np.deg2rad(theta_deg)
    phi = np.deg2rad(phi_deg)
    unit_spins = np.column_stack(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
    )
    spins = magnitudes[:, None] * unit_spins
    net_spin = np.sum(spins, axis=0)
    relative_net_spin = float(np.linalg.norm(net_spin) / np.sum(magnitudes))

    csv_path = args.csv or args.output.with_suffix(".csv")
    write_spin_csv(
        csv_path,
        labels,
        fractional_positions,
        cartesian_positions,
        magnitudes,
        theta_deg,
        phi_deg,
        spins,
    )

    figure = plt.figure(figsize=(15.2, 7.4), facecolor="white")
    real_axis = figure.add_subplot(1, 2, 1, projection="3d")
    sphere_axis = figure.add_subplot(1, 2, 2, projection="3d")
    color_map = plt.get_cmap("coolwarm")
    normalization = Normalize(vmin=-1.0, vmax=1.0)
    colors = color_map(normalization(unit_spins[:, 2]))

    draw_cell(real_axis, cell)
    real_axis.scatter(
        cartesian_positions[:, 0],
        cartesian_positions[:, 1],
        cartesian_positions[:, 2],
        s=12,
        color="#252b33",
        alpha=0.62,
        depthshade=True,
    )
    for position, direction, color in zip(
        cartesian_positions, unit_spins, colors
    ):
        real_axis.quiver(
            *position,
            *(args.arrow_length * direction),
            color=color,
            linewidth=1.15,
            arrow_length_ratio=0.30,
            normalize=False,
        )
    real_axis.set_title(
        "Spin Directions in the 4×4×4 BCC Supercell",
        fontsize=14,
        fontweight="bold",
        pad=14,
    )
    real_axis.set_xlabel("x (Å)")
    real_axis.set_ylabel("y (Å)")
    real_axis.set_zlabel("z (Å)")
    real_axis.view_init(elev=args.elev, azim=args.azim)
    set_equal_3d_axes(real_axis, np.vstack([np.zeros(3), cell]))

    polar = np.linspace(0.0, np.pi, 45)
    azimuth = np.linspace(0.0, 2.0 * np.pi, 80)
    polar_grid, azimuth_grid = np.meshgrid(polar, azimuth)
    sphere_x = np.sin(polar_grid) * np.cos(azimuth_grid)
    sphere_y = np.sin(polar_grid) * np.sin(azimuth_grid)
    sphere_z = np.cos(polar_grid)
    sphere_axis.plot_wireframe(
        sphere_x,
        sphere_y,
        sphere_z,
        rstride=5,
        cstride=8,
        color="#aeb7c4",
        linewidth=0.45,
        alpha=0.42,
    )
    sphere_axis.scatter(
        unit_spins[:, 0],
        unit_spins[:, 1],
        unit_spins[:, 2],
        c=unit_spins[:, 2],
        cmap=color_map,
        norm=normalization,
        s=32,
        edgecolor="white",
        linewidth=0.35,
        alpha=0.92,
        depthshade=False,
    )
    net_norm = float(np.linalg.norm(net_spin))
    if net_norm > 0.0:
        net_direction = net_spin / net_norm
        sphere_axis.quiver(
            0,
            0,
            0,
            *net_direction,
            color="#111827",
            linewidth=2.8,
            arrow_length_ratio=0.13,
            label="Net starting direction",
        )
    sphere_axis.set_title(
        "Orientation Distribution",
        fontsize=14,
        fontweight="bold",
        pad=14,
    )
    sphere_axis.set_xlabel(r"$S_x/|S|$")
    sphere_axis.set_ylabel(r"$S_y/|S|$")
    sphere_axis.set_zlabel(r"$S_z/|S|$")
    sphere_axis.set_xlim(-1.08, 1.08)
    sphere_axis.set_ylim(-1.08, 1.08)
    sphere_axis.set_zlim(-1.08, 1.08)
    sphere_axis.set_box_aspect((1, 1, 1))
    sphere_axis.view_init(elev=20, azim=-48)

    color_bar = figure.colorbar(
        plt.cm.ScalarMappable(norm=normalization, cmap=color_map),
        ax=[real_axis, sphere_axis],
        fraction=0.025,
        pad=0.02,
    )
    color_bar.set_label(r"Normalized $S_z$")
    figure.suptitle(
        "Noncollinear BCC Fe: Atom-Resolved Starting-Spin Texture",
        fontsize=19,
        fontweight="bold",
        color="#0b2d4d",
        y=0.98,
    )
    figure.text(
        0.5,
        0.025,
        (
            f"{len(labels)} Fe atoms  •  |Sᵢ| = {np.mean(magnitudes):.2f}  •  "
            f"|ΣSᵢ|/Σ|Sᵢ| = {relative_net_spin:.4f}  •  "
            f"source: {args.input_file.name}"
        ),
        ha="center",
        fontsize=10.5,
        color="#5b6472",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(figure)

    print(f"atoms={len(labels)}")
    print(f"net_spin=({net_spin[0]:.10f}, {net_spin[1]:.10f}, {net_spin[2]:.10f})")
    print(f"net_spin_magnitude={net_norm:.10f}")
    print(f"relative_net_spin={relative_net_spin:.10f}")
    print(f"figure={args.output}")
    print(f"csv={csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
