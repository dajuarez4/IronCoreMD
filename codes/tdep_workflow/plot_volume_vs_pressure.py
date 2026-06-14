#!/usr/bin/env python3
"""Plot phase-resolved volume versus pressure using TDEP free energies and QE MD pressures."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from tdep_common import (
    default_dataset_dir,
    free_energy_csv_name,
    pressure_csv_name,
    pressure_eos_plot_name,
    pressure_plot_name,
    read_csv_rows,
    resolve_path,
    source_npz_for_folder,
)
from tdep_phases import PHASE_SPECS, get_phase_spec

EV_A3_TO_GPA = 160.21766208


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot phase-resolved volume vs pressure from TDEP free energies.")
    parser.add_argument("--phase", choices=sorted(PHASE_SPECS), default="bcc", help="Crystal phase. Default: bcc.")
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=None,
        help="Directory containing NPZ/TDEP data. Default: <repo>/dataset/<phase>.",
    )
    parser.add_argument("--free-energy-csv", type=Path, default=None, help="CSV written by plot_free_energy_vs_volume.py.")
    parser.add_argument("--output", type=Path, default=None, help="Output PNG path.")
    parser.add_argument("--csv", type=Path, default=None, help="Output CSV path.")
    parser.add_argument("--eos-only-output", type=Path, default=None, help="EOS-only output PNG path.")
    parser.add_argument("--temperature-label", default="5000", help="Temperature label, e.g. 4500, 5000, 5500.")
    parser.add_argument(
        "--eos-label",
        default=None,
        help="Legend label for the EOS-only plot. Default: Iron <phase> at <temperature>K.",
    )
    return parser.parse_args()


def read_md_pressure(npz_path: Path) -> tuple[float, float]:
    data = np.load(npz_path, allow_pickle=False)
    pressure = np.asarray(data["pressure_GPa"], dtype=float)
    finite = pressure[np.isfinite(pressure)]
    if finite.size == 0:
        raise ValueError(f"No finite pressure values in {npz_path}")
    return float(finite.mean()), float(finite.std())


def birch_murnaghan_energy(volume: np.ndarray, e0: float, v0: float, b0: float, b0_prime: float) -> np.ndarray:
    eta = (v0 / volume) ** (2.0 / 3.0)
    strain = eta - 1.0
    return e0 + (9.0 * v0 * b0 / 16.0) * (b0_prime * strain**3 + strain**2 * (6.0 - 4.0 * eta))


def birch_murnaghan_pressure(volume: np.ndarray, v0: float, b0: float, b0_prime: float) -> np.ndarray:
    eta = (v0 / volume) ** (2.0 / 3.0)
    return 1.5 * b0 * (eta**3.5 - eta**2.5) * (1.0 + 0.75 * (b0_prime - 4.0) * (eta - 1.0))


def fit_birch_murnaghan(volumes: np.ndarray, free_energy: np.ndarray) -> tuple[float, float, float, float]:
    coeff = np.polyfit(volumes, free_energy, deg=2)
    if coeff[0] > 0.0:
        vertex = -coeff[1] / (2.0 * coeff[0])
        v0_guess = max(volumes.min() * 0.9, min(volumes.max() * 2.0, vertex))
        b0_guess = max(1.0e-3, v0_guess * 2.0 * coeff[0])
        e0_guess = float(np.polyval(coeff, v0_guess))
    else:
        v0_guess = float(volumes[np.argmin(free_energy)])
        b0_guess = 0.5
        e0_guess = float(np.min(free_energy))

    params, _ = curve_fit(
        birch_murnaghan_energy,
        volumes,
        free_energy,
        p0=(e0_guess, v0_guess, b0_guess, 4.0),
        bounds=(
            [free_energy.min() - 10.0, volumes.min() * 0.8, 1.0e-6, -5.0],
            [free_energy.max() + 10.0, volumes.max() * 5.0, 20.0, 20.0],
        ),
        maxfev=50000,
    )
    return tuple(float(value) for value in params)


def fit_birch_murnaghan_pressure(
    volumes: np.ndarray,
    pressures_gpa: np.ndarray,
    v0_guess: float,
    b0_gpa_guess: float,
    b0_prime_guess: float,
) -> tuple[float, float, float]:
    params, _ = curve_fit(
        birch_murnaghan_pressure,
        volumes,
        pressures_gpa,
        p0=(v0_guess, b0_gpa_guess, b0_prime_guess),
        bounds=(
            [volumes.min() * 0.8, 1.0e-6, -5.0],
            [volumes.max() * 5.0, 5000.0, 20.0],
        ),
        maxfev=50000,
    )
    return tuple(float(value) for value in params)


def collect_rows(dataset_dir: Path, free_energy_csv: Path) -> tuple[list[dict[str, float | str]], np.ndarray, np.ndarray]:
    rows_in = read_csv_rows(free_energy_csv)
    rows_in.sort(key=lambda row: float(row["volume_per_atom_A3"]))

    volumes = np.array([float(row["volume_per_atom_A3"]) for row in rows_in], dtype=float)
    free_energy = np.array([float(row["F_total_eV_atom"]) for row in rows_in], dtype=float)
    total_volumes = np.array([float(row["total_volume_A3"]) for row in rows_in], dtype=float)

    _, v0, b0, b0_prime = fit_birch_murnaghan(volumes, free_energy)
    pressure_fit = birch_murnaghan_pressure(volumes, v0, b0, b0_prime) * EV_A3_TO_GPA

    md_mean_pressure = []
    md_std_pressure = []
    for row_in in rows_in:
        mean_pressure, std_pressure = read_md_pressure(source_npz_for_folder(dataset_dir / row_in["folder"]))
        md_mean_pressure.append(mean_pressure)
        md_std_pressure.append(std_pressure)
    md_mean_pressure = np.array(md_mean_pressure, dtype=float)
    md_std_pressure = np.array(md_std_pressure, dtype=float)

    md_v0, md_b0_gpa, md_b0_prime = fit_birch_murnaghan_pressure(
        volumes,
        md_mean_pressure,
        v0_guess=v0,
        b0_gpa_guess=b0 * EV_A3_TO_GPA,
        b0_prime_guess=b0_prime,
    )
    md_pressure_fit = birch_murnaghan_pressure(volumes, md_v0, md_b0_gpa, md_b0_prime)

    dense_volumes = np.linspace(volumes.min(), volumes.max(), 500)
    dense_pressure = birch_murnaghan_pressure(dense_volumes, v0, b0, b0_prime) * EV_A3_TO_GPA
    dense_md_pressure = birch_murnaghan_pressure(dense_volumes, md_v0, md_b0_gpa, md_b0_prime)
    volume_scale = float(np.mean(total_volumes / volumes))
    dense_total_volumes = dense_volumes * volume_scale

    rows_out: list[dict[str, float | str]] = []
    for row_in, total_volume, pressure, mean_pressure, std_pressure, md_fit_pressure in zip(
        rows_in,
        total_volumes,
        pressure_fit,
        md_mean_pressure,
        md_std_pressure,
        md_pressure_fit,
    ):
        rows_out.append(
            {
                "folder": row_in["folder"],
                "lattice_a_A": float(row_in["lattice_a_A"]),
                "volume_per_atom_A3": float(row_in["volume_per_atom_A3"]),
                "total_volume_A3": float(total_volume),
                "F_total_eV_atom": float(row_in["F_total_eV_atom"]),
                "pressure_from_fit_GPa": float(pressure),
                "mean_md_pressure_GPa": mean_pressure,
                "std_md_pressure_GPa": std_pressure,
                "pressure_from_md_bm_fit_GPa": float(md_fit_pressure),
            }
        )

    return (
        rows_out,
        np.column_stack([dense_pressure, dense_total_volumes]),
        np.column_stack([dense_md_pressure, dense_total_volumes]),
    )


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    fieldnames = [
        "folder",
        "lattice_a_A",
        "volume_per_atom_A3",
        "total_volume_A3",
        "F_total_eV_atom",
        "pressure_from_fit_GPa",
        "mean_md_pressure_GPa",
        "std_md_pressure_GPa",
        "pressure_from_md_bm_fit_GPa",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot(
    path: Path,
    rows: list[dict[str, float | str]],
    dense_curve: np.ndarray,
    dense_md_curve: np.ndarray,
    temperature_label: str,
    phase: str,
) -> None:
    spec = get_phase_spec(phase)
    fit_pressure = np.array([float(row["pressure_from_fit_GPa"]) for row in rows], dtype=float)
    md_pressure = np.array([float(row["mean_md_pressure_GPa"]) for row in rows], dtype=float)
    md_std = np.array([float(row["std_md_pressure_GPa"]) for row in rows], dtype=float)
    volume = np.array([float(row["total_volume_A3"]) for row in rows], dtype=float)

    order_fit = np.argsort(fit_pressure)
    order_md = np.argsort(md_pressure)
    order_curve = np.argsort(dense_curve[:, 0])
    order_md_curve = np.argsort(dense_md_curve[:, 0])

    fig, ax = plt.subplots(figsize=(7.0, 6.2), constrained_layout=True)
    ax.plot(
        dense_curve[order_curve, 0],
        dense_curve[order_curve, 1],
        "--",
        color="black",
        linewidth=1.5,
        label=r"$P = (-\partial F / \partial V)$",
    )
    ax.plot(
        dense_md_curve[order_md_curve, 0],
        dense_md_curve[order_md_curve, 1],
        "-",
        color="#1f5f8a",
        linewidth=1.7,
        label="BM fit to AIMD mean pressure",
    )
    ax.plot(fit_pressure[order_fit], volume[order_fit], "o", color="#c24b2a", markersize=6.5)
    ax.errorbar(
        md_pressure[order_md],
        volume[order_md],
        xerr=md_std[order_md],
        fmt="s",
        ms=7.0,
        mfc="#1f5f8a",
        mec="#1f5f8a",
        ecolor="#1f5f8a",
        elinewidth=1.0,
        capsize=2.5,
        label="_nolegend_",
    )
    ax.set_xlabel("Pressure (GPa)")
    ax.set_ylabel(r"Volume ($\AA^3$)")
    ax.set_title(f"{spec.title} Fe Volume vs Pressure at {temperature_label} K")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, loc="best")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_eos_only(path: Path, rows: list[dict[str, float | str]], dense_md_curve: np.ndarray, eos_label: str) -> None:
    md_pressure = np.array([float(row["mean_md_pressure_GPa"]) for row in rows], dtype=float)
    md_std = np.array([float(row["std_md_pressure_GPa"]) for row in rows], dtype=float)
    volume = np.array([float(row["total_volume_A3"]) for row in rows], dtype=float)
    marker_size = 7.2

    order_md = np.argsort(md_pressure)
    order_curve = np.argsort(dense_md_curve[:, 0])

    fig, ax = plt.subplots(figsize=(7.0, 6.2), constrained_layout=True)
    all_pressures = np.concatenate([md_pressure, dense_md_curve[:, 0]])
    pressure_span = float(all_pressures.max() - all_pressures.min())
    volume_span = float(volume.max() - volume.min())
    xpad = 0.05 * pressure_span if pressure_span > 0.0 else 5.0
    ypad = 0.05 * volume_span if volume_span > 0.0 else 0.2
    ax.set_xlabel("Pressure (GPa)", fontsize=20)
    ax.set_ylabel(r"Volume ($\AA^3$)", fontsize=20)
    ax.set_xlim(all_pressures.min() - xpad, all_pressures.max() + xpad)
    ax.set_ylim(volume.min() - ypad, volume.max() + ypad)
    ax.tick_params(axis="both", labelsize=16)
    ax.grid(True, color="#d9c7f0", alpha=0.45)
    fig.canvas.draw()

    sorted_volume = volume[order_md]
    axes_bbox = ax.get_window_extent()
    y_range = float(ax.get_ylim()[1] - ax.get_ylim()[0])
    marker_height_pixels = marker_size * 0.98 * fig.dpi / 72.0
    box_height = marker_height_pixels * (y_range / axes_bbox.height)

    ax.plot(
        dense_md_curve[order_curve, 0],
        dense_md_curve[order_curve, 1],
        ":",
        color="black",
        linewidth=1.5,
        label="_nolegend_",
    )
    ax.barh(
        sorted_volume,
        2.0 * md_std[order_md],
        left=md_pressure[order_md] - md_std[order_md],
        height=box_height,
        color="#d95c5c",
        alpha=0.6,
        edgecolor="none",
        label="_nolegend_",
    )
    ax.plot(
        md_pressure[order_md],
        sorted_volume,
        "s",
        ms=marker_size,
        mfc="#1f5aa6",
        mec="#123b73",
        mew=0.9,
        label=eos_label,
    )
    ax.legend(frameon=False, loc="upper right", fontsize=14)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    spec = get_phase_spec(args.phase)
    dataset_dir = args.dataset_dir.resolve() if args.dataset_dir is not None else default_dataset_dir(args.phase).resolve()
    free_energy_csv = resolve_path(
        dataset_dir,
        args.free_energy_csv if args.free_energy_csv is not None else free_energy_csv_name(args.temperature_label),
    )
    output = resolve_path(
        dataset_dir,
        args.output if args.output is not None else pressure_plot_name(args.temperature_label, args.phase),
    )
    csv_path = resolve_path(
        dataset_dir,
        args.csv if args.csv is not None else pressure_csv_name(args.temperature_label, args.phase),
    )
    eos_output = resolve_path(
        dataset_dir,
        args.eos_only_output if args.eos_only_output is not None else pressure_eos_plot_name(args.temperature_label, args.phase),
    )
    eos_label = args.eos_label if args.eos_label else f"Iron {spec.long_name} at {args.temperature_label}K"

    rows, dense_curve, dense_md_curve = collect_rows(dataset_dir, free_energy_csv)
    write_csv(csv_path, rows)
    plot(output, rows, dense_curve, dense_md_curve, args.temperature_label, args.phase)
    plot_eos_only(eos_output, rows, dense_md_curve, eos_label)

    print(f"Wrote {csv_path}")
    print(f"Wrote {output}")
    print(f"Wrote {eos_output}")


if __name__ == "__main__":
    main()
