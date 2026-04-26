#!/usr/bin/env python3
"""
Create a No-Vito-style GIF from compressed QE AIMD NPZ data.

The script expects the NPZ format produced by the QE parser used in this
workspace. It renders atoms in 3D, draws the simulation cell, and overlays a
small info box with timestep, time, temperature, and pressure for each frame.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter

BOHR_TO_ANG = 0.529177210903


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render a GIF from QE AIMD NPZ data.")
    parser.add_argument("npz", type=Path, help="Path to the input NPZ archive.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output GIF path. Default: <npz_stem>_qe_md.gif beside the input file.",
    )
    parser.add_argument("--start", type=int, default=0, help="First frame index to render.")
    parser.add_argument("--stop", type=int, default=0, help="Stop before this frame index. 0 means all.")
    parser.add_argument("--every", type=int, default=1, help="Render every Nth frame.")
    parser.add_argument("--fps", type=int, default=12, help="GIF frames per second.")
    parser.add_argument("--dpi", type=int, default=140, help="Output DPI.")
    parser.add_argument("--figsize", type=float, nargs=2, default=(7.8, 7.2), metavar=("W", "H"))
    parser.add_argument("--elev", type=float, default=17.0, help="3D camera elevation.")
    parser.add_argument("--azim", type=float, default=-62.0, help="3D camera azimuth.")
    parser.add_argument("--atom-size", type=float, default=36.0, help="Atom marker size.")
    parser.add_argument("--title", type=str, default=None, help="Optional plot title.")
    return parser.parse_args()


def load_metadata(raw) -> dict:
    value = raw.item() if hasattr(raw, "item") else raw
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return json.loads(str(value))


def safe_text(value: float, fmt: str, suffix: str = "") -> str:
    if value is None or not np.isfinite(value):
        return f"n/a{suffix}"
    return f"{format(float(value), fmt)}{suffix}"


def color_map(symbols: np.ndarray) -> dict[str, str]:
    base = {
        "Fe": "#e76f51",
        "Ni": "#f4a261",
        "Ti": "#4d96ff",
        "C": "#9aa0a6",
        "O": "#8ecae6",
    }
    fallback = plt.get_cmap("tab10")
    unique = sorted({str(sym) for sym in symbols.tolist()})
    cmap = {}
    for i, sym in enumerate(unique):
        cmap[sym] = base.get(sym, fallback(i % 10))
    return cmap


def convert_cell_to_angstrom(cell: np.ndarray, unit: str, alat_bohr: float | None) -> np.ndarray:
    norm = (unit or "").strip().lower()
    if norm in {"angstrom", "ang"}:
        return np.asarray(cell, dtype=float)
    if norm in {"bohr", "a.u.", "au"}:
        return np.asarray(cell, dtype=float) * BOHR_TO_ANG
    if norm == "alat":
        if alat_bohr is None:
            raise ValueError("Found CELL_PARAMETERS in alat but alat_bohr is unavailable.")
        return np.asarray(cell, dtype=float) * alat_bohr * BOHR_TO_ANG
    if norm in {"none", ""}:
        return np.asarray(cell, dtype=float)
    raise ValueError(f"Unsupported cell unit: {unit!r}")


def fixed_cell_angstrom(data, meta: dict) -> np.ndarray:
    alat_bohr = meta.get("alat_bohr")
    if "cell_parameters" in data and "cell_parameters_unit" in data:
        cells = np.asarray(data["cell_parameters"], dtype=float)
        units = np.asarray(data["cell_parameters_unit"])
        for idx in range(len(cells)):
            cell = cells[idx]
            if np.isfinite(cell).all():
                return convert_cell_to_angstrom(cell, str(units[idx]), alat_bohr)

    if "input_cell_parameters" in data and "input_cell_unit" in data:
        return convert_cell_to_angstrom(
            np.asarray(data["input_cell_parameters"], dtype=float),
            str(data["input_cell_unit"]),
            alat_bohr,
        )

    if "initial_cell_alat" in data and alat_bohr is not None:
        return np.asarray(data["initial_cell_alat"], dtype=float) * float(alat_bohr) * BOHR_TO_ANG

    raise ValueError("Could not infer a valid simulation cell from the NPZ file.")


def positions_to_cartesian(frame_positions: np.ndarray, unit: str, cell_ang: np.ndarray, alat_bohr: float | None) -> np.ndarray:
    norm = (unit or "").strip().lower()
    pos = np.asarray(frame_positions, dtype=float)

    if norm == "crystal":
        frac = pos - np.floor(pos)
        return frac @ cell_ang
    if norm in {"angstrom", "ang"}:
        return pos
    if norm in {"bohr", "a.u.", "au"}:
        return pos * BOHR_TO_ANG
    if norm == "alat":
        if alat_bohr is None:
            raise ValueError("Found positions in alat but alat_bohr is unavailable.")
        return pos * float(alat_bohr) * BOHR_TO_ANG
    raise ValueError(f"Unsupported positions unit: {unit!r}")


def cell_vertices(cell_ang: np.ndarray) -> np.ndarray:
    a, b, c = np.asarray(cell_ang, dtype=float)
    origin = np.zeros(3)
    return np.array(
        [
            origin,
            a,
            b,
            c,
            a + b,
            a + c,
            b + c,
            a + b + c,
        ],
        dtype=float,
    )


def plot_cell_edges(ax, cell_ang: np.ndarray) -> None:
    a, b, c = np.asarray(cell_ang, dtype=float)
    origin = np.zeros(3)
    edges = [
        (origin, a),
        (origin, b),
        (origin, c),
        (a, a + b),
        (a, a + c),
        (b, a + b),
        (b, b + c),
        (c, a + c),
        (c, b + c),
        (a + b, a + b + c),
        (a + c, a + b + c),
        (b + c, a + b + c),
    ]
    for p0, p1 in edges:
        ax.plot(
            [p0[0], p1[0]],
            [p0[1], p1[1]],
            [p0[2], p1[2]],
            color="white",
            linewidth=0.8,
            alpha=0.45,
        )


def set_axes_limits(ax, cell_ang: np.ndarray) -> None:
    verts = cell_vertices(cell_ang)
    mins = verts.min(axis=0)
    maxs = verts.max(axis=0)
    span = maxs - mins
    max_span = float(np.max(span))
    if max_span <= 0:
        max_span = 1.0
    centers = 0.5 * (mins + maxs)
    half = 0.58 * max_span

    ax.set_xlim(centers[0] - half, centers[0] + half)
    ax.set_ylim(centers[1] - half, centers[1] + half)
    ax.set_zlim(centers[2] - half, centers[2] + half)
    ax.set_box_aspect((1.0, 1.0, 1.0))


def draw_frame(
    ax,
    cart_positions: np.ndarray,
    cell_ang: np.ndarray,
    symbols: np.ndarray,
    symbol_colors: dict[str, str],
    atom_size: float,
    title: str,
    info_lines: list[str],
    elev: float,
    azim: float,
) -> None:
    ax.clear()
    ax.set_facecolor("black")
    set_axes_limits(ax, cell_ang)
    plot_cell_edges(ax, cell_ang)

    for sym in sorted({str(s) for s in symbols.tolist()}):
        mask = symbols == sym
        ax.scatter(
            cart_positions[mask, 0],
            cart_positions[mask, 1],
            cart_positions[mask, 2],
            s=atom_size,
            alpha=0.95,
            color=symbol_colors[sym],
            edgecolors="white",
            linewidths=0.25,
            depthshade=True,
        )

    ax.set_title(title, color="white", pad=14, fontsize=13)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    ax.text2D(
        0.03,
        0.97,
        "\n".join(info_lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        color="white",
        fontsize=10.5,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "black",
            "edgecolor": "white",
            "alpha": 0.72,
        },
    )


def main() -> None:
    args = parse_args()
    if args.every <= 0:
        raise ValueError("--every must be a positive integer.")

    npz_path = args.npz.resolve()
    if not npz_path.exists():
        raise FileNotFoundError(npz_path)

    out_path = args.output.resolve() if args.output else npz_path.with_name(f"{npz_path.stem}_qe_md.gif")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    data = np.load(npz_path, allow_pickle=False)
    meta = load_metadata(data["metadata_json"])
    symbols = np.asarray(data["species"])
    positions = np.asarray(data["positions"], dtype=float)
    units = np.asarray(data["positions_unit"])
    iterations = np.asarray(data["iteration"], dtype=int) if "iteration" in data else np.arange(1, len(positions) + 1)
    time_ps = np.asarray(data["time_ps"], dtype=float) if "time_ps" in data else np.full(len(positions), np.nan)
    temperature = np.asarray(data["temperature_K"], dtype=float) if "temperature_K" in data else np.full(len(positions), np.nan)
    pressure = np.asarray(data["pressure_GPa"], dtype=float) if "pressure_GPa" in data else np.full(len(positions), np.nan)

    stop = args.stop if args.stop > 0 else len(positions)
    frame_indices = np.arange(args.start, min(stop, len(positions)), args.every, dtype=int)
    if len(frame_indices) == 0:
        raise ValueError("No frames selected. Check --start/--stop/--every.")

    cell_ang = fixed_cell_angstrom(data, meta)
    alat_bohr = meta.get("alat_bohr")
    title = args.title or f"QE AIMD: {npz_path.stem}"
    symbol_colors = color_map(symbols)

    fig = plt.figure(figsize=tuple(args.figsize))
    fig.patch.set_facecolor("black")
    ax = fig.add_subplot(111, projection="3d")

    writer = PillowWriter(fps=args.fps)
    with writer.saving(fig, str(out_path), dpi=args.dpi):
        for iframe in frame_indices:
            cart = positions_to_cartesian(positions[iframe], str(units[iframe]), cell_ang, alat_bohr)
            info_lines = [
                f"Timestep: {int(iterations[iframe])}",
                f"Time: {safe_text(time_ps[iframe], '.3f', ' ps')}",
                f"T: {safe_text(temperature[iframe], '.1f', ' K')}",
                f"P: {safe_text(pressure[iframe], '.2f', ' GPa')}",
            ]
            draw_frame(
                ax=ax,
                cart_positions=cart,
                cell_ang=cell_ang,
                symbols=symbols,
                symbol_colors=symbol_colors,
                atom_size=args.atom_size,
                title=title,
                info_lines=info_lines,
                elev=args.elev,
                azim=args.azim,
            )
            writer.grab_frame()

    plt.close(fig)
    print(f"Wrote {out_path}")
    print(f"Frames: {len(frame_indices)}")


if __name__ == "__main__":
    main()
