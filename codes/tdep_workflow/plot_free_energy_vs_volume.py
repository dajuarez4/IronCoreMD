#!/usr/bin/env python3
"""Plot BCC TDEP free energy versus volume and lattice parameter."""

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
    classify_free_energy,
    default_dataset_dir,
    default_tdep_folders,
    free_energy_csv_name,
    free_energy_lattice_plot_name,
    free_energy_plot_name,
    prefer_unique_lattice_points,
    read_free_energy,
    read_u0_second_order,
    read_uc_volume_per_atom,
    relative_free_energy_lattice_plot_name,
    relative_free_energy_plot_name,
    resolve_path,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot BCC TDEP free energy vs volume and lattice parameter.")
    parser.add_argument("folders", nargs="*", type=Path, help="TDEP folders. Default: all folders at temperature-label.")
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=default_dataset_dir(),
        help="Directory containing the BCC TDEP folders. Default: <repo>/dataset/bcc.",
    )
    parser.add_argument("--output", type=Path, default=None, help="Absolute or dataset-relative output PNG path.")
    parser.add_argument("--relative-output", type=Path, default=None, help="Relative free-energy output PNG path.")
    parser.add_argument("--csv", type=Path, default=None, help="Output CSV path.")
    parser.add_argument("--lattice-output", type=Path, default=None, help="Lattice-parameter free-energy PNG path.")
    parser.add_argument(
        "--relative-lattice-output",
        type=Path,
        default=None,
        help="Relative free-energy vs lattice PNG path.",
    )
    parser.add_argument("--temperature-label", default="5000", help="Temperature label, e.g. 4500, 5000, 5500.")
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    fieldnames = [
        "folder",
        "lattice_a_A",
        "volume_per_atom_A3",
        "total_volume_A3",
        "T_K",
        "F_vib_eV_atom",
        "U0_2nd_eV_atom",
        "F_total_eV_atom",
        "F_relative_meV_atom",
        "S_eV_K_atom",
        "Cv_eV_K_atom",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def birch_murnaghan_energy(volume: np.ndarray, e0: float, v0: float, b0: float, b0_prime: float) -> np.ndarray:
    eta = (v0 / volume) ** (2.0 / 3.0)
    strain = eta - 1.0
    return e0 + (9.0 * v0 * b0 / 16.0) * (b0_prime * strain**3 + strain**2 * (6.0 - 4.0 * eta))


def fit_birch_murnaghan(volumes: np.ndarray, free_energies: np.ndarray) -> tuple[float, float, float, float]:
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


def collect_rows(folders: list[Path]) -> tuple[list[dict[str, float | str]], list[tuple[str, str]]]:
    rows: list[dict[str, float | str]] = []
    skipped: list[tuple[str, str]] = []
    for folder in folders:
        lattice_a, volume_per_atom, total_volume = read_uc_volume_per_atom(folder / "infile.ucposcar")
        temperature, f_vib, entropy, cv = read_free_energy(folder / "outfile.free_energy")
        u0 = read_u0_second_order(folder / "outfile.U0")
        status = classify_free_energy(temperature, f_vib, entropy, cv, u0)
        if status != "ok":
            skipped.append((folder.name, status))
            continue
        rows.append(
            {
                "folder": folder.name,
                "lattice_a_A": lattice_a,
                "volume_per_atom_A3": volume_per_atom,
                "total_volume_A3": total_volume,
                "T_K": temperature,
                "F_vib_eV_atom": f_vib,
                "U0_2nd_eV_atom": u0,
                "F_total_eV_atom": u0 + f_vib,
                "S_eV_K_atom": entropy,
                "Cv_eV_K_atom": cv,
            }
        )
    if not rows:
        raise ValueError("No valid TDEP folders were available for plotting.")
    rows.sort(key=lambda row: float(row["volume_per_atom_A3"]))
    fmin = min(float(row["F_total_eV_atom"]) for row in rows)
    for row in rows:
        row["F_relative_meV_atom"] = (float(row["F_total_eV_atom"]) - fmin) * 1000.0
    return rows, skipped


def save_single_panel(
    path: Path,
    x_points: np.ndarray,
    y_points: np.ndarray,
    x_fit: np.ndarray,
    y_fit: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
    point_color: str,
    fit_color: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 5.2), constrained_layout=True)
    ax.scatter(x_points, y_points, color=point_color, s=75, zorder=3)
    ax.plot(x_fit, y_fit, "-", color=fit_color, linewidth=1.8, label="Birch-Murnaghan fit")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_vs_volume(path: Path, relative_path: Path, rows: list[dict[str, float | str]], temperature_label: str) -> None:
    volumes = np.array([float(row["total_volume_A3"]) for row in rows], dtype=float)
    free_energies = np.array([float(row["F_total_eV_atom"]) for row in rows], dtype=float)
    e0_fit, v0_fit, b0_fit, b0_prime_fit = fit_birch_murnaghan(volumes, free_energies)
    fit_volumes = np.linspace(volumes.min(), volumes.max(), 400)
    fit_free_energies = birch_murnaghan_energy(fit_volumes, e0_fit, v0_fit, b0_fit, b0_prime_fit)
    relative = (free_energies - e0_fit) * 1000.0
    fit_relative = (fit_free_energies - e0_fit) * 1000.0

    save_single_panel(
        path=path,
        x_points=volumes,
        y_points=free_energies,
        x_fit=fit_volumes,
        y_fit=fit_free_energies,
        xlabel=r"Volume ($\AA^3$)",
        ylabel="Free Helmholtz energy (eV)",
        title=f"BCC Fe Free Energy vs Volume at {temperature_label} K",
        point_color="#1f6f5b",
        fit_color="#c24b2a",
    )
    save_single_panel(
        path=relative_path,
        x_points=volumes,
        y_points=relative,
        x_fit=fit_volumes,
        y_fit=fit_relative,
        xlabel=r"Volume ($\AA^3$)",
        ylabel="Relative free energy (meV/atom)",
        title=f"BCC Fe Relative Free Energy vs Volume at {temperature_label} K",
        point_color="#1f6f5b",
        fit_color="#c24b2a",
    )


def plot_vs_lattice(path: Path, relative_path: Path, rows: list[dict[str, float | str]], temperature_label: str) -> None:
    lattice_a = np.array([float(row["lattice_a_A"]) for row in rows], dtype=float)
    free_energies = np.array([float(row["F_total_eV_atom"]) for row in rows], dtype=float)
    volumes = np.array([float(row["volume_per_atom_A3"]) for row in rows], dtype=float)
    e0_fit, v0_fit, b0_fit, b0_prime_fit = fit_birch_murnaghan(volumes, free_energies)
    fit_volumes = np.linspace(volumes.min(), volumes.max(), 400)
    fit_lattice_a = (2.0 * fit_volumes) ** (1.0 / 3.0)
    fit_free_energies = birch_murnaghan_energy(fit_volumes, e0_fit, v0_fit, b0_fit, b0_prime_fit)
    relative = (free_energies - e0_fit) * 1000.0
    fit_relative = (fit_free_energies - e0_fit) * 1000.0

    save_single_panel(
        path=path,
        x_points=lattice_a,
        y_points=free_energies,
        x_fit=fit_lattice_a,
        y_fit=fit_free_energies,
        xlabel="Lattice parameter, a (A)",
        ylabel="Free Helmholtz energy (eV)",
        title=f"BCC Fe Free Energy vs Lattice Parameter at {temperature_label} K",
        point_color="#305c9b",
        fit_color="#b6572d",
    )
    save_single_panel(
        path=relative_path,
        x_points=lattice_a,
        y_points=relative,
        x_fit=fit_lattice_a,
        y_fit=fit_relative,
        xlabel="Lattice parameter, a (A)",
        ylabel="Relative free energy (meV/atom)",
        title=f"BCC Fe Relative Free Energy vs Lattice Parameter at {temperature_label} K",
        point_color="#305c9b",
        fit_color="#b6572d",
    )


def main() -> None:
    args = parse_args()
    dataset_dir = args.dataset_dir.resolve()
    folders = [path.resolve() if path.is_absolute() else (dataset_dir / path).resolve() for path in args.folders]
    if not folders:
        folders = default_tdep_folders(dataset_dir, args.temperature_label)
    folders = prefer_unique_lattice_points(folders)
    if not folders:
        raise FileNotFoundError(f"No TDEP folders found for {args.temperature_label} K in {dataset_dir}")

    output = resolve_path(dataset_dir, args.output if args.output is not None else free_energy_plot_name(args.temperature_label))
    relative_output = resolve_path(
        dataset_dir,
        args.relative_output if args.relative_output is not None else relative_free_energy_plot_name(args.temperature_label),
    )
    csv_path = resolve_path(dataset_dir, args.csv if args.csv is not None else free_energy_csv_name(args.temperature_label))
    lattice_output = resolve_path(
        dataset_dir,
        args.lattice_output if args.lattice_output is not None else free_energy_lattice_plot_name(args.temperature_label),
    )
    relative_lattice_output = resolve_path(
        dataset_dir,
        args.relative_lattice_output
        if args.relative_lattice_output is not None
        else relative_free_energy_lattice_plot_name(args.temperature_label),
    )

    rows, skipped = collect_rows(folders)
    write_csv(csv_path, rows)
    plot_vs_volume(output, relative_output, rows, args.temperature_label)
    plot_vs_lattice(lattice_output, relative_lattice_output, rows, args.temperature_label)

    for folder_name, status in skipped:
        print(f"Skipped {folder_name}: {status}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {output}")
    print(f"Wrote {relative_output}")
    print(f"Wrote {lattice_output}")
    print(f"Wrote {relative_lattice_output}")


if __name__ == "__main__":
    main()
