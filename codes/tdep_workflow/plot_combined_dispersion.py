#!/usr/bin/env python3
"""Overlay phase-resolved TDEP phonon dispersions and DOS in one figure."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec

from tdep_common import (
    classify_free_energy,
    default_dataset_dir,
    default_tdep_folders,
    dispersion_plot_name,
    folder_case_label,
    folder_sort_key,
    lattice_parameter_from_folder,
    prefer_unique_lattice_points,
    read_free_energy,
    read_u0_second_order,
    resolve_path,
)
from tdep_phases import PHASE_SPECS, get_phase_spec


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Overlay TDEP phonon dispersions from several folders.")
    parser.add_argument("folders", nargs="*", type=Path, help="TDEP folders. Default: all valid folders at temperature-label.")
    parser.add_argument("--phase", choices=sorted(PHASE_SPECS), default="bcc", help="Crystal phase. Default: bcc.")
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=None,
        help="Directory containing the TDEP folders. Default: <repo>/dataset/<phase>.",
    )
    parser.add_argument("--output", type=Path, default=None, help="Output PNG path.")
    parser.add_argument("--temperature-label", default="5000", help="Temperature label, e.g. 4500, 5000, 5500.")
    return parser.parse_args()


def choose_reference_folder(folders: list[Path], temperature_label: str, phase: str) -> Path:
    del temperature_label
    ordered = sorted(folders, key=lambda folder: folder_sort_key(folder, phase))
    return ordered[len(ordered) // 2]


def read_dispersion(path: Path) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected dispersion format in {path}")
    return data


def read_total_dos(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected DOS format in {path}")
    return data[:, 0], data[:, 1]


def read_xticks(path: Path) -> tuple[list[str], np.ndarray]:
    labels: list[str] = []
    positions: list[float] = []
    pattern = re.compile(r'\(\s*"([^"]+)"\s+([0-9Ee+\-.]+)\s*\)')
    for line in path.read_text().splitlines():
        match = pattern.search(line)
        if match:
            labels.append(match.group(1))
            positions.append(float(match.group(2)))
    if not labels:
        raise ValueError(f"Could not parse xtics from {path}")
    return labels, np.array(positions, dtype=float)


def rescale_path(x: np.ndarray, source_ticks: np.ndarray, reference_ticks: np.ndarray) -> np.ndarray:
    if len(source_ticks) != len(reference_ticks):
        raise ValueError("Tick count mismatch while aligning path coordinates")
    return np.interp(x, source_ticks, reference_ticks)


def folder_has_valid_free_energy(folder: Path) -> bool:
    temperature, f_vib, entropy, cv = read_free_energy(folder / "outfile.free_energy")
    u0 = read_u0_second_order(folder / "outfile.U0")
    return classify_free_energy(temperature, f_vib, entropy, cv, u0) == "ok"


def plot(output: Path, folders: list[Path], temperature_label: str, phase: str) -> None:
    spec = get_phase_spec(phase)
    reference_folder = choose_reference_folder(folders, temperature_label, phase)
    ref_labels, ref_ticks = read_xticks(reference_folder / "outfile.dispersion_relations.gnuplot")

    plt.rcParams.update({"font.family": "serif", "font.size": 12.5})
    if spec.key == "hcp":
        fig = plt.figure(figsize=(10.6, 6.5), constrained_layout=True)
        grid = GridSpec(1, 2, figure=fig, width_ratios=[4.9, 1.35], wspace=0.06)
        ax = fig.add_subplot(grid[0, 0])
        ax_dos = fig.add_subplot(grid[0, 1], sharey=ax)
        ax_cbar = None
    else:
        fig = plt.figure(figsize=(10.9, 6.5), constrained_layout=True)
        grid = GridSpec(1, 3, figure=fig, width_ratios=[4.9, 1.35, 0.14], wspace=0.08)
        ax = fig.add_subplot(grid[0, 0])
        ax_dos = fig.add_subplot(grid[0, 1], sharey=ax)
        ax_cbar = fig.add_subplot(grid[0, 2])

    cmap = plt.get_cmap("viridis")
    if spec.key == "hcp":
        folder_colors = {
            folder: cmap(index / max(1, len(folders) - 1))
            for index, folder in enumerate(sorted(folders, key=lambda folder: folder_sort_key(folder, phase)))
        }
        norm = None
        lattice_parameters = None
    else:
        lattice_parameters = np.array([lattice_parameter_from_folder(folder, phase=phase) for folder in folders], dtype=float)
        norm = Normalize(vmin=float(np.min(lattice_parameters)), vmax=float(np.max(lattice_parameters)))
        folder_colors = {folder: cmap(norm(lattice_parameter_from_folder(folder, phase=phase))) for folder in folders}
    min_frequency = 0.0
    max_frequency = 0.0
    max_dos = 0.0

    for folder in folders:
        data = read_dispersion(folder / "outfile.dispersion_relations")
        _, folder_ticks = read_xticks(folder / "outfile.dispersion_relations.gnuplot")
        x_aligned = rescale_path(data[:, 0], folder_ticks, ref_ticks)
        color = folder_colors[folder]
        min_frequency = min(min_frequency, float(np.min(data[:, 1:])))
        max_frequency = max(max_frequency, float(np.max(data[:, 1:])))

        for band_index in range(1, data.shape[1]):
            ax.plot(
                x_aligned,
                data[:, band_index],
                color=color,
                linewidth=1.45 if spec.key == "hcp" else 1.15,
                alpha=0.94 if spec.key == "hcp" else 0.9,
                label=folder_case_label(folder, phase) if spec.key == "hcp" and band_index == 1 else None,
            )

        dos_frequency, total_dos = read_total_dos(folder / "outfile.phonon_dos")
        max_frequency = max(max_frequency, float(np.max(dos_frequency)))
        max_dos = max(max_dos, float(np.max(total_dos)))
        ax_dos.plot(total_dos, dos_frequency, color=color, linewidth=1.55, alpha=0.95)

    for xpos in ref_ticks:
        ax.axvline(xpos, color="#888888", linewidth=0.8, alpha=0.45, zorder=0)

    ax.axhline(0.0, color="black", linewidth=0.85, alpha=0.8)
    ax_dos.axhline(0.0, color="black", linewidth=0.85, alpha=0.8)
    ax.set_xlim(ref_ticks[0], ref_ticks[-1])
    ax.set_xticks(ref_ticks)
    ax.set_xticklabels(ref_labels)
    ymin = min_frequency * 1.08 if min_frequency < 0.0 else 0.0
    ax.set_ylim(ymin, max_frequency * 1.03)
    ax.set_ylabel("Frequency (THz)")
    ax.set_title(f"{spec.title} Fe Phonon Dispersion and DOS at {temperature_label} K")
    ax.grid(axis="y", alpha=0.25)

    ax_dos.set_xlim(0.0, max_dos * 1.08)
    ax_dos.set_xlabel("DOS")
    ax_dos.set_title("Total DOS")
    ax_dos.grid(axis="y", alpha=0.2)
    ax_dos.tick_params(axis="y", which="both", left=False, labelleft=False)

    if spec.key == "hcp":
        ax.legend(frameon=False, loc="upper right")
    else:
        colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=ax_cbar)
        colorbar.set_label("a (Å)")
        colorbar.set_ticks(np.linspace(lattice_parameters[0], lattice_parameters[-1], 6))

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    dataset_dir = args.dataset_dir.resolve() if args.dataset_dir is not None else default_dataset_dir(args.phase).resolve()
    folders = [path.resolve() if path.is_absolute() else (dataset_dir / path).resolve() for path in args.folders]
    if not folders:
        folders = default_tdep_folders(dataset_dir, args.temperature_label, phase=args.phase)
    folders = prefer_unique_lattice_points(folders, phase=args.phase)
    if not folders:
        raise FileNotFoundError(f"No TDEP folders found for {args.temperature_label} K in {dataset_dir}")

    output = resolve_path(
        dataset_dir,
        args.output if args.output is not None else dispersion_plot_name(args.temperature_label, phase=args.phase),
    )
    plot(output, folders, args.temperature_label, args.phase)

    unstable_folders = [folder.name for folder in folders if not folder_has_valid_free_energy(folder)]
    for folder_name in unstable_folders:
        print(f"Included {folder_name}: invalid free energy / imaginary modes retained in phonon overlay")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
