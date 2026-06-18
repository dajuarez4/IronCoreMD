#!/usr/bin/env python3
"""Overlay phase-resolved free-energy and pressure-volume curves for multiple temperatures."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from tdep_common import (
    DEFAULT_COMPARISON_TEMPERATURES,
    comparison_output_name,
    default_dataset_dir,
    discover_temperature_series,
    free_energy_csv_name,
    pressure_csv_name,
    read_csv_rows,
    resolve_path,
)
from tdep_phases import PHASE_SPECS, get_phase_spec

SERIES_STYLES = {
    "4000": {"point_color": "#8f1d3f", "line_color": "#5d0f28", "box_color": "#e9bcc9", "marker": "P"},
    "4500": {"point_color": "#c65f1a", "line_color": "#8f3d05", "box_color": "#f0c9b2", "marker": "s"},
    "5000": {"point_color": "#1f5aa6", "line_color": "#123b73", "box_color": "#bfd1ec", "marker": "o"},
    "5500": {"point_color": "#2f7d4a", "line_color": "#1c5732", "box_color": "#bfdcbc", "marker": "^"},
    "6000": {"point_color": "#7a3fc0", "line_color": "#4f2285", "box_color": "#d7c0ef", "marker": "D"},
}
FALLBACK_STYLES = [
    {"point_color": "#b04b30", "line_color": "#7d2f1c", "box_color": "#e9c0b6", "marker": "v"},
    {"point_color": "#4e8a2d", "line_color": "#2f5c18", "box_color": "#cee4be", "marker": "P"},
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Overlay phase-resolved free-energy and pressure-volume curves.")
    parser.add_argument("--phase", choices=sorted(PHASE_SPECS), default="bcc", help="Crystal phase. Default: bcc.")
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=None,
        help="Directory containing thermodynamic CSVs. Default: <repo>/dataset/<phase>.",
    )
    parser.add_argument(
        "--temperatures",
        nargs="+",
        default=None,
        help="Temperature labels to overlay, e.g. 4500 5000 5500. Default: all available CSV pairs.",
    )
    parser.add_argument("--free-energy-output", type=Path, default=None)
    parser.add_argument("--pressure-output", type=Path, default=None)
    return parser.parse_args()


def normalize_temperatures(labels: list[str]) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    for label in labels:
        clean = str(label).removesuffix("K")
        if clean not in seen:
            seen.add(clean)
            ordered.append(clean)
    return ordered


def style_for_temperature(temperature_label: str, fallback_index: int) -> dict[str, str]:
    if temperature_label in SERIES_STYLES:
        return SERIES_STYLES[temperature_label]
    if fallback_index < len(FALLBACK_STYLES):
        return FALLBACK_STYLES[fallback_index]
    return FALLBACK_STYLES[-1]


def build_series(dataset_dir: Path, temperatures: list[str], phase: str) -> list[dict[str, object]]:
    series = []
    for fallback_index, temperature in enumerate(temperatures):
        style = style_for_temperature(temperature, fallback_index)
        series.append(
            {
                "temperature_label": temperature,
                "label": f"{temperature} K",
                "free_energy_rows": read_csv_rows(resolve_path(dataset_dir, free_energy_csv_name(temperature))),
                "pressure_rows": read_csv_rows(resolve_path(dataset_dir, pressure_csv_name(temperature, phase))),
                **style,
            }
        )
    return series


def birch_murnaghan_energy(volume: np.ndarray, e0: float, v0: float, b0: float, b0_prime: float) -> np.ndarray:
    eta = (v0 / volume) ** (2.0 / 3.0)
    strain = eta - 1.0
    return e0 + (9.0 * v0 * b0 / 16.0) * (b0_prime * strain**3 + strain**2 * (6.0 - 4.0 * eta))


def birch_murnaghan_pressure(volume: np.ndarray, v0: float, b0: float, b0_prime: float) -> np.ndarray:
    eta = (v0 / volume) ** (2.0 / 3.0)
    return 1.5 * b0 * (eta**3.5 - eta**2.5) * (1.0 + 0.75 * (b0_prime - 4.0) * (eta - 1.0))


def fit_birch_murnaghan_energy(volumes: np.ndarray, free_energies: np.ndarray) -> tuple[float, float, float, float]:
    coeff = np.polyfit(volumes, free_energies, deg=2)
    if coeff[0] > 0.0:
        vertex = -coeff[1] / (2.0 * coeff[0])
        v0_guess = max(volumes.min() * 0.9, min(volumes.max() * 2.0, vertex))
        b0_guess = max(1.0e-3, v0_guess * 2.0 * coeff[0])
        e0_guess = float(np.polyval(coeff, v0_guess))
    else:
        v0_guess = float(volumes[np.argmin(free_energies)])
        b0_guess = 0.5
        e0_guess = float(np.min(free_energies))

    params, _ = curve_fit(
        birch_murnaghan_energy,
        volumes,
        free_energies,
        p0=(e0_guess, v0_guess, b0_guess, 4.0),
        bounds=(
            [free_energies.min() - 10.0, volumes.min() * 0.8, 1.0e-6, -5.0],
            [free_energies.max() + 10.0, volumes.max() * 5.0, 20.0, 20.0],
        ),
        maxfev=50000,
    )
    return tuple(float(value) for value in params)


def fit_birch_murnaghan_pressure(volumes: np.ndarray, pressures_gpa: np.ndarray) -> tuple[float, float, float]:
    coeff = np.polyfit(volumes, pressures_gpa, deg=2)
    if coeff[0] != 0.0:
        vertex = -coeff[1] / (2.0 * coeff[0])
        v0_guess = max(volumes.min() * 0.9, min(volumes.max() * 2.0, vertex))
    else:
        v0_guess = float(volumes[np.argmin(pressures_gpa)])
    params, _ = curve_fit(
        birch_murnaghan_pressure,
        volumes,
        pressures_gpa,
        p0=(v0_guess, 100.0, 4.0),
        bounds=(
            [volumes.min() * 0.8, 1.0e-6, -5.0],
            [volumes.max() * 5.0, 5000.0, 20.0],
        ),
        maxfev=50000,
    )
    return tuple(float(value) for value in params)


def plot_free_energy(path: Path, series: list[dict[str, object]], phase: str) -> None:
    spec = get_phase_spec(phase)
    plt.rcParams.update({"font.family": "serif", "font.size": 12.5})
    fig, ax = plt.subplots(figsize=(7.3, 6.1), constrained_layout=True)

    all_volumes = []
    all_energies = []
    for item in series:
        rows = item["free_energy_rows"]
        volumes = np.array([float(row["total_volume_A3"]) for row in rows], dtype=float)
        free_energies = np.array([float(row["F_total_eV_atom"]) for row in rows], dtype=float)
        order = np.argsort(volumes)
        e0_fit, v0_fit, b0_fit, b0_prime_fit = fit_birch_murnaghan_energy(volumes, free_energies)
        fit_volumes = np.linspace(volumes.min(), volumes.max(), 400)
        fit_free_energies = birch_murnaghan_energy(fit_volumes, e0_fit, v0_fit, b0_fit, b0_prime_fit)

        ax.plot(fit_volumes, fit_free_energies, "-", color=item["line_color"], linewidth=1.9, label="_nolegend_")
        ax.scatter(
            volumes[order],
            free_energies[order],
            s=60,
            color=item["point_color"],
            edgecolors=item["line_color"],
            linewidths=0.9,
            marker=item["marker"],
            zorder=3,
            label=item["label"],
        )
        all_volumes.append(volumes)
        all_energies.append(free_energies)

    volumes_all = np.concatenate(all_volumes)
    energies_all = np.concatenate(all_energies)
    xspan = float(volumes_all.max() - volumes_all.min())
    yspan = float(energies_all.max() - energies_all.min())
    xpad = 0.03 * xspan if xspan > 0.0 else 0.2
    ypad = 0.05 * yspan if yspan > 0.0 else 0.02

    ax.set_title(f"{spec.title} Fe Free Helmholtz Energy vs Volume")
    ax.set_xlabel(r"Conventional-cell volume ($\AA^3$)", fontsize=18)
    ax.set_ylabel("Free Helmholtz energy (eV/atom)", fontsize=18)
    ax.set_xlim(volumes_all.min() - xpad, volumes_all.max() + xpad)
    ax.set_ylim(energies_all.min() - ypad, energies_all.max() + ypad)
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(True, color="#d9c7f0", alpha=0.35)
    ax.legend(frameon=False, loc="upper right", fontsize=13, ncol=min(3, len(series)), columnspacing=0.9)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def add_pressure_series(
    ax: plt.Axes,
    rows: list[dict[str, str]],
    *,
    label: str,
    point_face: str,
    point_edge: str,
    box_color: str,
    line_color: str,
    marker: str,
    marker_size: float,
    box_height: float,
) -> None:
    md_pressure = np.array([float(row["mean_md_pressure_GPa"]) for row in rows], dtype=float)
    md_std = np.array([float(row["std_md_pressure_GPa"]) for row in rows], dtype=float)
    volume = np.array([float(row["total_volume_A3"]) for row in rows], dtype=float)

    order = np.argsort(md_pressure)
    sorted_pressure = md_pressure[order]
    sorted_std = md_std[order]
    sorted_volume = volume[order]

    v0_fit, b0_fit, b0_prime_fit = fit_birch_murnaghan_pressure(volume, md_pressure)
    fit_volumes = np.linspace(volume.min(), volume.max(), 400)
    fit_pressures = birch_murnaghan_pressure(fit_volumes, v0_fit, b0_fit, b0_prime_fit)
    fit_order = np.argsort(fit_pressures)

    ax.plot(fit_pressures[fit_order], fit_volumes[fit_order], ":", color=line_color, linewidth=1.8, label="_nolegend_", zorder=1)
    ax.barh(
        sorted_volume,
        2.0 * sorted_std,
        left=sorted_pressure - sorted_std,
        height=box_height,
        color=box_color,
        alpha=0.58,
        edgecolor="none",
        label="_nolegend_",
        zorder=2,
    )
    ax.plot(
        sorted_pressure,
        sorted_volume,
        marker,
        ms=marker_size,
        mfc=point_face,
        mec=point_edge,
        mew=0.9,
        label=label,
        zorder=3,
    )


def plot_pressure(path: Path, series: list[dict[str, object]], phase: str) -> None:
    spec = get_phase_spec(phase)
    fig, ax = plt.subplots(figsize=(7.0, 6.2), constrained_layout=True)
    marker_size = 7.0

    pressure_arrays = []
    volume_arrays = []
    for item in series:
        rows = item["pressure_rows"]
        pressure_arrays.append(np.array([float(row["mean_md_pressure_GPa"]) for row in rows], dtype=float))
        volume_arrays.append(np.array([float(row["total_volume_A3"]) for row in rows], dtype=float))
    pressures_all = np.concatenate(pressure_arrays)
    volumes_all = np.concatenate(volume_arrays)
    xspan = float(pressures_all.max() - pressures_all.min())
    yspan = float(volumes_all.max() - volumes_all.min())
    xpad = 0.05 * xspan if xspan > 0.0 else 5.0
    ypad = 0.05 * yspan if yspan > 0.0 else 0.2

    ax.set_xlabel("Pressure (GPa)", fontsize=20)
    ax.set_ylabel(r"Volume ($\AA^3$)", fontsize=20)
    ax.set_xlim(pressures_all.min() - xpad, pressures_all.max() + xpad)
    ax.set_ylim(volumes_all.min() - ypad, volumes_all.max() + ypad)
    ax.tick_params(axis="both", labelsize=16)
    ax.grid(True, color="#d9c7f0", alpha=0.45)
    fig.canvas.draw()

    axes_bbox = ax.get_window_extent()
    y_range = float(ax.get_ylim()[1] - ax.get_ylim()[0])
    marker_height_pixels = marker_size * 0.98 * fig.dpi / 72.0
    box_height = marker_height_pixels * (y_range / axes_bbox.height)

    for item in series:
        add_pressure_series(
            ax,
            item["pressure_rows"],
            label=f"Iron {spec.long_name} at {item['temperature_label']}K",
            point_face=item["point_color"],
            point_edge=item["line_color"],
            box_color=item["box_color"],
            line_color=item["line_color"],
            marker=item["marker"],
            marker_size=marker_size,
            box_height=box_height,
        )
    ax.legend(frameon=False, loc="best", fontsize=13)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    dataset_dir = args.dataset_dir.resolve() if args.dataset_dir is not None else default_dataset_dir(args.phase).resolve()
    temperatures = (
        normalize_temperatures(args.temperatures)
        if args.temperatures is not None
        else discover_temperature_series(dataset_dir, args.phase)
    )
    if not temperatures:
        temperatures = list(DEFAULT_COMPARISON_TEMPERATURES)
    free_energy_output = resolve_path(
        dataset_dir,
        args.free_energy_output
        if args.free_energy_output is not None
        else comparison_output_name("free_energy_vs_volume", temperatures, ".png"),
    )
    pressure_output = resolve_path(
        dataset_dir,
        args.pressure_output
        if args.pressure_output is not None
        else comparison_output_name("volume_vs_pressure", temperatures, f"_{args.phase}.png"),
    )

    series = build_series(dataset_dir, temperatures, args.phase)
    plot_free_energy(free_energy_output, series, args.phase)
    plot_pressure(pressure_output, series, args.phase)
    print(f"Wrote {free_energy_output}")
    print(f"Wrote {pressure_output}")


if __name__ == "__main__":
    main()
