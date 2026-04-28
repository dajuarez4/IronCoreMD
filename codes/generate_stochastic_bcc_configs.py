#!/usr/bin/env python3
"""Generate stochastic BCC configurations from an existing TDEP folder.

This uses TDEP's canonical_configuration utility directly on the fitted
second-order force constants. For the current BCC workflow this is cleaner
than the downloaded phonopy scripts because the data is already in TDEP
format (`infile.ucposcar`, `infile.ssposcar`, `infile.forceconstant`).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


OUTPUT_FORMATS = {
    "vasp": 1,
    "abinit": 2,
    "aims": 4,
    "siesta": 5,
    "qe": 6,
    "parsec": 7,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate stochastic canonical configurations from a BCC TDEP folder."
    )
    parser.add_argument(
        "folder",
        nargs="?",
        type=Path,
        default=Path("tdep_2.40_5000K"),
        help="BCC TDEP folder. Default: tdep_2.40_5000K in the parent directory.",
    )
    parser.add_argument(
        "--temperature",
        "-t",
        type=float,
        default=None,
        help="Sampling temperature in K. Default: read from infile.meta.",
    )
    parser.add_argument(
        "--nconf",
        "-n",
        type=int,
        default=20,
        help="Number of configurations to generate. Default: 20.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where generated files are written. Default: <folder>/stochastic_qe_<T>K.",
    )
    parser.add_argument(
        "--format",
        choices=sorted(OUTPUT_FORMATS),
        default="qe",
        help="Output configuration format. Default: qe.",
    )
    parser.add_argument(
        "--quantum",
        action="store_true",
        help="Use Bose-Einstein statistics instead of classical Maxwell-Boltzmann.",
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=None,
        help="Path to TDEP canonical_configuration binary. Default: repo-local build.",
    )
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def resolve_tdep_folder(folder: Path) -> Path:
    return folder.resolve() if folder.is_absolute() else (Path(__file__).resolve().parents[1] / folder).resolve()


def read_temperature_from_meta(path: Path) -> float:
    lines = path.read_text().splitlines()
    if len(lines) < 4:
        raise ValueError(f"Expected at least 4 lines in {path}")
    return float(lines[3].split()[0])


def require_inputs(folder: Path) -> None:
    needed = ["infile.ucposcar", "infile.ssposcar", "infile.forceconstant"]
    missing = [name for name in needed if not (folder / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required input files in {folder}: {', '.join(missing)}")


def recreate_directory(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def stage_inputs(src: Path, dst: Path) -> None:
    for name in ("infile.ucposcar", "infile.ssposcar", "infile.forceconstant"):
        (dst / name).symlink_to((src / name).resolve())


def generated_files(path: Path) -> list[Path]:
    prefixes = ("contcar_conf", "qe_conf", "aims_conf", "parsec_conf", "siesta_conf")
    return sorted(p for p in path.iterdir() if p.is_file() and p.name.startswith(prefixes))


def main() -> None:
    args = parse_args()
    folder = resolve_tdep_folder(args.folder)
    require_inputs(folder)

    temperature = args.temperature
    if temperature is None:
        temperature = read_temperature_from_meta(folder / "infile.meta")

    output_dir = args.output_dir
    if output_dir is None:
        output_dir = folder / f"stochastic_{args.format}_{int(round(temperature))}K"
    elif not output_dir.is_absolute():
        output_dir = (folder / output_dir).resolve()
    else:
        output_dir = output_dir.resolve()

    binary = args.binary
    if binary is None:
        binary = repo_root() / "tdep/build/src/canonical_configuration/canonical_configuration"
    binary = binary.resolve()
    if not binary.exists():
        raise FileNotFoundError(f"Could not find canonical_configuration binary at {binary}")

    recreate_directory(output_dir)
    stage_inputs(folder, output_dir)

    command = [
        str(binary),
        "-n",
        str(args.nconf),
        "-t",
        str(temperature),
        "-of",
        str(OUTPUT_FORMATS[args.format]),
    ]
    if args.quantum:
        command.append("--quantum")

    subprocess.run(command, cwd=output_dir, check=True)

    produced = generated_files(output_dir)
    print(f"Source folder: {folder}")
    print(f"Output folder: {output_dir}")
    print(f"Temperature (K): {temperature}")
    print(f"Configurations requested: {args.nconf}")
    print(f"Configurations written: {len(produced)}")
    if produced:
        print(f"First file: {produced[0].name}")
        print(f"Last file: {produced[-1].name}")


if __name__ == "__main__":
    main()
