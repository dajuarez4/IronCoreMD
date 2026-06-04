#!/usr/bin/env python3
"""Run the harmonic BCC TDEP workflow from QE NPZ archives."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from npz_to_tdep_bcc import sanitize_stem, write_tdep_folder
from tdep_common import (
    default_dataset_dir,
    discover_bcc_npz_files,
    find_tdep_root,
    free_energy_csv_name,
    free_energy_lattice_plot_name,
    free_energy_plot_name,
    normalize_temperature_label,
    pressure_csv_name,
    pressure_eos_plot_name,
    pressure_plot_name,
    relative_free_energy_lattice_plot_name,
    relative_free_energy_plot_name,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the harmonic BCC TDEP workflow from QE NPZ files.")
    parser.add_argument(
        "targets",
        nargs="*",
        help="NPZ file names, NPZ stems, or tdep_* folder names. Default: all QE NPZ files at temperature-label.",
    )
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=default_dataset_dir(),
        help="Directory containing QE NPZ archives and TDEP folders. Default: <repo>/dataset/bcc.",
    )
    parser.add_argument("--temperature-label", default="5000", help="Temperature label, e.g. 4500, 5000, 5500.")
    parser.add_argument(
        "--temperature-K",
        type=float,
        default=None,
        help="Temperature written into infile.meta/stat. Default: float(temperature-label).",
    )
    parser.add_argument("--skip", type=int, default=0, help="Skip this many initial frames when rebuilding TDEP inputs.")
    parser.add_argument("--every", type=int, default=1, help="Use every Nth frame after --skip.")
    parser.add_argument("--max-frames", type=int, default=0, help="Maximum frames after filtering. 0 means all.")
    parser.add_argument("--keep-invalid", action="store_true", help="Keep frames containing NaN/Inf values.")
    parser.add_argument("--supercell", type=int, nargs=3, default=(4, 4, 4), metavar=("NX", "NY", "NZ"))
    parser.add_argument("--rc2", type=float, default=5.0, help="Force-constant cutoff passed to extract_forceconstants.")
    parser.add_argument("--dos-qgrid", type=int, nargs=3, default=(32, 32, 32), metavar=("QX", "QY", "QZ"))
    parser.add_argument(
        "--tdep-root",
        type=Path,
        default=None,
        help="Path to the TDEP build/src directory. Default: auto-detect repo-local or sibling checkout.",
    )
    parser.add_argument("--no-convert", action="store_true", help="Skip the NPZ -> TDEP input regeneration step.")
    parser.add_argument("--no-tdep", action="store_true", help="Skip the TDEP binary execution step.")
    parser.add_argument("--no-plots", action="store_true", help="Skip the plot regeneration step.")
    parser.add_argument(
        "--comparison-temperatures",
        nargs="+",
        default=None,
        help="Temperature labels passed to the comparison plot script. Default: auto-discover available CSV pairs.",
    )
    parser.add_argument(
        "--no-comparison-plots",
        action="store_true",
        help="Refresh only the single-temperature plots, not the multi-temperature comparison.",
    )
    return parser.parse_args()


def resolve_targets(dataset_dir: Path, targets: list[str], temperature_label: str) -> list[Path]:
    if not targets:
        return discover_bcc_npz_files(dataset_dir, temperature_label)

    resolved: list[Path] = []
    for target in targets:
        path = Path(target)
        if path.is_absolute():
            candidate = path
        else:
            name = path.name
            if name.startswith("tdep_"):
                name = name.removeprefix("tdep_")
            if not name.endswith(".npz"):
                name = f"{name}.npz"
            candidate = dataset_dir / name
        resolved.append(candidate.resolve())
    return resolved


def tdep_folder_for_npz(dataset_dir: Path, npz_path: Path) -> Path:
    return dataset_dir / f"tdep_{sanitize_stem(npz_path.stem)}"


def ensure_forceconstant_link(folder: Path) -> None:
    link = folder / "infile.forceconstant"
    target = folder / "outfile.forceconstant"
    if link.exists() or link.is_symlink():
        link.unlink()
    link.symlink_to(target.resolve())


def run_logged(command: list[str], cwd: Path, log_path: Path) -> None:
    with log_path.open("w") as handle:
        subprocess.run(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, check=True)


def refresh_single_temperature_plots(dataset_dir: Path, temperature_label: str) -> None:
    script_dir = Path(__file__).resolve().parent
    commands = [
        [
            sys.executable,
            str(script_dir / "plot_free_energy_vs_volume.py"),
            "--dataset-dir",
            str(dataset_dir),
            "--temperature-label",
            temperature_label,
            "--output",
            str(free_energy_plot_name(temperature_label)),
            "--relative-output",
            str(relative_free_energy_plot_name(temperature_label)),
            "--csv",
            str(free_energy_csv_name(temperature_label)),
            "--lattice-output",
            str(free_energy_lattice_plot_name(temperature_label)),
            "--relative-lattice-output",
            str(relative_free_energy_lattice_plot_name(temperature_label)),
        ],
        [
            sys.executable,
            str(script_dir / "plot_volume_vs_pressure.py"),
            "--dataset-dir",
            str(dataset_dir),
            "--temperature-label",
            temperature_label,
            "--free-energy-csv",
            str(free_energy_csv_name(temperature_label)),
            "--output",
            str(pressure_plot_name(temperature_label)),
            "--csv",
            str(pressure_csv_name(temperature_label)),
            "--eos-only-output",
            str(pressure_eos_plot_name(temperature_label)),
        ],
        [
            sys.executable,
            str(script_dir / "plot_combined_dispersion.py"),
            "--dataset-dir",
            str(dataset_dir),
            "--temperature-label",
            temperature_label,
        ],
    ]
    for command in commands:
        subprocess.run(command, check=True)


def refresh_comparison_plot(dataset_dir: Path, temperatures: list[str]) -> None:
    script_dir = Path(__file__).resolve().parent
    command = [
        sys.executable,
        str(script_dir / "plot_temperature_comparison.py"),
        "--dataset-dir",
        str(dataset_dir),
    ]
    if temperatures:
        command.extend(["--temperatures", *temperatures])
    subprocess.run(command, check=True)


def main() -> None:
    args = parse_args()
    dataset_dir = args.dataset_dir.resolve()
    temperature_label = normalize_temperature_label(args.temperature_label)
    comparison_temperatures = (
        [normalize_temperature_label(item) for item in args.comparison_temperatures]
        if args.comparison_temperatures is not None
        else []
    )
    temperature_override = args.temperature_K if args.temperature_K is not None else float(temperature_label)
    npz_files = resolve_targets(dataset_dir, args.targets, temperature_label)
    tdep_folders = [tdep_folder_for_npz(dataset_dir, npz_path) for npz_path in npz_files]

    if not args.no_convert:
        for npz_path, folder in zip(npz_files, tdep_folders):
            write_tdep_folder(
                npz_path=npz_path,
                outdir=folder,
                supercell=tuple(args.supercell),
                temperature_override=temperature_override,
                skip=args.skip,
                every=args.every,
                max_frames=args.max_frames,
                keep_invalid=args.keep_invalid,
            )
            print(f"[convert] {npz_path.name} -> {folder}")

    if not args.no_tdep:
        tdep_root = find_tdep_root(args.tdep_root)
        extract_bin = tdep_root / "extract_forceconstants" / "extract_forceconstants"
        dispersion_bin = tdep_root / "phonon_dispersion_relations" / "phonon_dispersion_relations"
        for folder in tdep_folders:
            print(f"[tdep] {folder.name}")
            run_logged([str(extract_bin), "-rc2", str(args.rc2)], folder, folder / "extract_forceconstants.log")
            ensure_forceconstant_link(folder)
            run_logged([str(dispersion_bin)], folder, folder / "phonon_dispersion_relations.log")
            run_logged(
                [
                    str(dispersion_bin),
                    "--dos",
                    "--qpoint_grid",
                    *(str(value) for value in args.dos_qgrid),
                    "--temperature",
                    temperature_label,
                ],
                folder,
                folder / f"free_energy_{temperature_label}K.log",
            )

    if not args.no_plots:
        refresh_single_temperature_plots(dataset_dir, temperature_label)
        if not args.no_comparison_plots:
            refresh_comparison_plot(dataset_dir, comparison_temperatures)


if __name__ == "__main__":
    main()
