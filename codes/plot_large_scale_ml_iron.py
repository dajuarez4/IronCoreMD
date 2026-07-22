#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize the scale of a large BCC iron simulation enabled by ML."
    )
    parser.add_argument("--cells", type=int, default=18)
    parser.add_argument("--lattice-parameter", type=float, default=2.55)
    parser.add_argument("--output-dir", type=Path, default=Path("dataset/figures"))
    parser.add_argument("--poster-dir", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=200)
    return parser.parse_args()


def bcc_supercell(cells, lattice_parameter):
    grid = np.indices((cells, cells, cells), dtype=np.float64).reshape(3, -1).T
    corners = grid
    centers = grid + 0.5
    positions = np.vstack((corners, centers)) * lattice_parameter
    return positions


def displaced_positions(positions, lattice_parameter):
    phase = positions / lattice_parameter
    displacement = np.column_stack(
        (
            np.sin(0.63 * phase[:, 1] + 0.41 * phase[:, 2]),
            np.sin(0.57 * phase[:, 2] + 0.37 * phase[:, 0] + 1.2),
            np.sin(0.61 * phase[:, 0] + 0.43 * phase[:, 1] + 2.1),
        )
    )
    return positions + 0.075 * lattice_parameter * displacement


def project_isometric(positions):
    azimuth = np.deg2rad(43.0)
    elevation = np.deg2rad(24.0)
    rotation_z = np.array(
        [
            [np.cos(azimuth), -np.sin(azimuth), 0.0],
            [np.sin(azimuth), np.cos(azimuth), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    rotation_x = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(elevation), -np.sin(elevation)],
            [0.0, np.sin(elevation), np.cos(elevation)],
        ]
    )
    rotated = positions @ rotation_z.T @ rotation_x.T
    depth = rotated[:, 2]
    perspective = 1.0 + 0.0028 * (depth - depth.min())
    return rotated[:, 0] * perspective, rotated[:, 1] * perspective, depth


def plot_large_scale(cells, lattice_parameter, output_base, dpi):
    positions = bcc_supercell(cells, lattice_parameter)
    positions = displaced_positions(positions, lattice_parameter)
    projected_x, projected_y, depth = project_isometric(positions)
    order = np.argsort(depth)
    projected_x = projected_x[order]
    projected_y = projected_y[order]
    depth = depth[order]

    colormap = LinearSegmentedColormap.from_list(
        "iron_depth",
        ["#17375E", "#087E8B", "#56CFE1", "#FFD166", "#F05D5E"],
    )
    normalized_depth = (depth - depth.min()) / (np.ptp(depth) + 1.0e-12)
    atom_sizes = 3.0 + 9.0 * normalized_depth**1.7

    figure = plt.figure(figsize=(19.2, 10.8), facecolor="#07111E")
    axis = figure.add_axes((0.015, 0.015, 0.97, 0.97), facecolor="#07111E")
    axis.scatter(
        projected_x,
        projected_y,
        c=normalized_depth,
        cmap=colormap,
        s=atom_sizes,
        alpha=0.88,
        linewidths=0,
        rasterized=True,
    )

    center = np.array([projected_x.mean(), projected_y.mean()])
    radius = 0.105 * max(np.ptp(projected_x), np.ptp(projected_y))
    distances = np.sqrt((projected_x - center[0]) ** 2 + (projected_y - center[1]) ** 2)
    highlighted = distances < radius
    axis.scatter(
        projected_x[highlighted],
        projected_y[highlighted],
        s=atom_sizes[highlighted] * 2.1,
        facecolors="none",
        edgecolors="#FFF2B2",
        linewidths=0.35,
        alpha=0.70,
    )

    axis.set_aspect("equal")
    axis.set_axis_off()
    padding_x = np.ptp(projected_x) * 0.05
    padding_y = np.ptp(projected_y) * 0.07
    axis.set_xlim(projected_x.min() - padding_x, projected_x.max() + padding_x)
    axis.set_ylim(projected_y.min() - padding_y, projected_y.max() + padding_y)

    number_atoms = len(positions)

    for suffix in ("png", "pdf", "svg"):
        figure.savefig(
            output_base.with_suffix(f".{suffix}"),
            dpi=dpi if suffix == "png" else None,
            facecolor=figure.get_facecolor(),
        )
    plt.close(figure)
    return number_atoms


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_base = args.output_dir / "iron_ml_large_scale_simulation_space"
    number_atoms = plot_large_scale(
        args.cells,
        args.lattice_parameter,
        output_base,
        args.dpi,
    )

    if args.poster_dir is not None:
        args.poster_dir.mkdir(parents=True, exist_ok=True)
        for suffix in ("png", "pdf"):
            source = output_base.with_suffix(f".{suffix}")
            shutil.copy2(source, args.poster_dir / source.name)

    print(f"Visualized {number_atoms:,} BCC Fe atoms")
    print(f"Figure: {output_base.with_suffix('.png')}")


if __name__ == "__main__":
    main()
