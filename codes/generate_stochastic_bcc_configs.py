#!/usr/bin/env python3
"""Generate stochastic BCC configurations from an existing TDEP folder.

This uses TDEP's canonical_configuration utility directly on the fitted
second-order force constants. For the current BCC workflow this is cleaner
than the downloaded phonopy scripts because the data is already in TDEP
format (`infile.ucposcar`, `infile.ssposcar`, `infile.forceconstant`).
"""

from __future__ import annotations

import argparse
import math
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
        "--mindist",
        type=float,
        default=None,
        help=(
            "Minimum allowed pair distance in units of the nearest-neighbor "
            "distance. Passed to TDEP --mindist."
        ),
    )
    parser.add_argument(
        "--binary",
        type=Path,
        default=None,
        help="Path to TDEP canonical_configuration binary. Default: repo-local build.",
    )
    parser.add_argument(
        "--supercell",
        nargs=3,
        type=int,
        metavar=("NX", "NY", "NZ"),
        default=None,
        help=(
            "Optional BCC conventional-cell repetitions. When provided, write a new "
            "infile.ssposcar instead of reusing the source supercell."
        ),
    )
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_tdep_folder(folder: Path) -> Path:
    if folder.is_absolute():
        return folder.resolve()
    cwd_candidate = (Path.cwd() / folder).resolve()
    if cwd_candidate.exists():
        return cwd_candidate
    return (Path(__file__).resolve().parents[1] / folder).resolve()


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


def read_bcc_lattice_parameter(ucposcar: Path) -> float:
    lines = ucposcar.read_text().splitlines()
    if len(lines) < 5:
        raise ValueError(f"Incomplete unit-cell POSCAR: {ucposcar}")
    scale = float(lines[1].split()[0])
    first_vector = [float(value) for value in lines[2].split()[:3]]
    primitive_length = scale * math.sqrt(sum(value * value for value in first_vector))
    return 2.0 * primitive_length / math.sqrt(3.0)


def write_bcc_conventional_supercell(
    path: Path,
    lattice_parameter_ang: float,
    repetitions: tuple[int, int, int],
) -> None:
    nx, ny, nz = repetitions
    if min(repetitions) <= 0:
        raise ValueError(f"Supercell repetitions must be positive, got {repetitions}.")
    positions: list[tuple[float, float, float]] = []
    for offset in (0.0, 0.5):
        for z_index in range(nz):
            for y_index in range(ny):
                for x_index in range(nx):
                    positions.append(
                        (
                            (x_index + offset) / nx,
                            (y_index + offset) / ny,
                            (z_index + offset) / nz,
                        )
                    )
    lines = [
        f"Fe bcc ideal {nx}x{ny}x{nz} conventional supercell",
        "1.0",
        f"{nx * lattice_parameter_ang: .16f}  0.0000000000000000  0.0000000000000000",
        f"0.0000000000000000  {ny * lattice_parameter_ang: .16f}  0.0000000000000000",
        f"0.0000000000000000  0.0000000000000000  {nz * lattice_parameter_ang: .16f}",
        "Fe",
        str(len(positions)),
        "Direct",
    ]
    lines.extend(f"{x: .16f} {y: .16f} {z: .16f}" for x, y, z in positions)
    path.write_text("\n".join(lines) + "\n")


def stage_inputs(
    src: Path,
    dst: Path,
    supercell: tuple[int, int, int] | None = None,
) -> None:
    for name in ("infile.ucposcar", "infile.forceconstant"):
        (dst / name).symlink_to((src / name).resolve())
    if supercell is None:
        (dst / "infile.ssposcar").symlink_to((src / "infile.ssposcar").resolve())
    else:
        lattice_parameter = read_bcc_lattice_parameter(src / "infile.ucposcar")
        write_bcc_conventional_supercell(
            dst / "infile.ssposcar",
            lattice_parameter_ang=lattice_parameter,
            repetitions=supercell,
        )


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

    supercell = None if args.supercell is None else tuple(args.supercell)
    recreate_directory(output_dir)
    stage_inputs(folder, output_dir, supercell=supercell)

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
    if args.mindist is not None:
        command.extend(["--mindist", str(args.mindist)])

    subprocess.run(command, cwd=output_dir, check=True)

    produced = generated_files(output_dir)
    print(f"Source folder: {folder}")
    print(f"Output folder: {output_dir}")
    print(f"Temperature (K): {temperature}")
    print(f"Configurations requested: {args.nconf}")
    print(f"Configurations written: {len(produced)}")
    if args.mindist is not None:
        print(f"Minimum-distance ratio: {args.mindist}")
    if supercell is not None:
        print(f"BCC conventional supercell: {supercell[0]} x {supercell[1]} x {supercell[2]}")
        print(f"Atoms per configuration: {2 * math.prod(supercell)}")
    if produced:
        print(f"First file: {produced[0].name}")
        print(f"Last file: {produced[-1].name}")


if __name__ == "__main__":
    main()
