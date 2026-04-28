#!/usr/bin/env python3
"""Wrap stochastic QE position files into full BCC QE MD inputs."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


CONFIG_RE = re.compile(r"qe_conf(\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build full QE MD inputs from stochastic qe_conf files for BCC Fe."
    )
    parser.add_argument(
        "config_dir",
        nargs="?",
        type=Path,
        default=Path("tdep_2.40_5000K/stochastic_qe_5000K"),
        help="Directory containing qe_conf#### files. Default: tdep_2.40_5000K/stochastic_qe_5000K",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("md_inputs"),
        help="Directory for generated QE input files. Default: md_inputs inside config_dir.",
    )
    parser.add_argument("--temperature", "-t", type=float, default=5000.0, help="Target MD temperature in K.")
    parser.add_argument("--nstep", type=int, default=400, help="Number of MD steps. Default: 400.")
    parser.add_argument("--dt", type=float, default=20.670, help="QE dt in atomic units. Default: 20.670.")
    parser.add_argument(
        "--pseudo-dir",
        type=Path,
        default=None,
        help="Pseudo directory path. Default: repository root containing the Fe UPF file.",
    )
    parser.add_argument("--prefix-base", type=str, default="Fe_bcc_4x4x4_no_spin_stochastic_5000K")
    parser.add_argument("--ecutwfc", type=float, default=70.0)
    parser.add_argument("--ecutrho", type=float, default=560.0)
    parser.add_argument("--degauss", type=float, default=0.03)
    parser.add_argument(
        "--kgrid",
        nargs=3,
        type=int,
        default=(4, 4, 4),
        metavar=("KX", "KY", "KZ"),
        help="Automatic k-point grid. Default: 4 4 4.",
    )
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def resolve_dir(base: Path, value: Path) -> Path:
    return value.resolve() if value.is_absolute() else (base / value).resolve()


def read_qe_conf(path: Path) -> tuple[int, int, str, list[str], list[str]]:
    lines = path.read_text().splitlines()

    nat = None
    ntyp = None
    cell_header_index = None
    positions_header_index = None
    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("nat"):
            nat = int(stripped.split("=")[1].strip().rstrip(","))
        elif stripped.startswith("ntyp"):
            ntyp = int(stripped.split("=")[1].strip().rstrip(","))
        elif stripped.startswith("CELL_PARAMETERS"):
            cell_header_index = index
        elif stripped.startswith("ATOMIC_POSITIONS"):
            positions_header_index = index

    if nat is None or ntyp is None or cell_header_index is None or positions_header_index is None:
        raise ValueError(f"Could not parse {path}")

    cell_header = lines[cell_header_index].strip()
    cell_lines = lines[cell_header_index + 1 : cell_header_index + 4]
    position_header = lines[positions_header_index].strip()
    position_lines = []
    for line in lines[positions_header_index + 1 :]:
        if not line.strip():
            break
        position_lines.append(line.rstrip())

    if len(cell_lines) != 3:
        raise ValueError(f"Expected 3 CELL_PARAMETERS lines in {path}")
    if len(position_lines) != nat:
        raise ValueError(f"Expected {nat} atomic positions in {path}, found {len(position_lines)}")

    return nat, ntyp, cell_header, cell_lines, [position_header, *position_lines]


def conf_id(path: Path) -> str:
    match = CONFIG_RE.fullmatch(path.name)
    if not match:
        raise ValueError(f"Unexpected stochastic config name: {path.name}")
    return match.group(1)


def build_input_text(
    nat: int,
    ntyp: int,
    cell_header: str,
    cell_lines: list[str],
    position_block: list[str],
    prefix: str,
    pseudo_dir: str,
    temperature: float,
    nstep: int,
    dt: float,
    ecutwfc: float,
    ecutrho: float,
    degauss: float,
    kgrid: tuple[int, int, int],
) -> str:
    kx, ky, kz = kgrid
    pieces = [
        "&CONTROL",
        "    calculation   = 'md'",
        "    restart_mode  = 'from_scratch'",
        f"    prefix        = '{prefix}'",
        f"    pseudo_dir    = '{pseudo_dir}'",
        "    outdir        = './tmp/'",
        f"    dt            = {dt:.4f}d0",
        f"    nstep         = {nstep}",
        "    tstress       = .true.",
        "    tprnfor       = .true.",
        "/",
        "",
        "&SYSTEM",
        "    ibrav = 0",
        f"    nat   = {nat}",
        f"    ntyp  = {ntyp}",
        f"    ecutwfc     = {ecutwfc:.1f}",
        f"    ecutrho     = {ecutrho:.1f}",
        "    occupations = 'smearing'",
        "    smearing    = 'mv'",
        f"    degauss     = {degauss:.2f}",
        "    nosym       = .true.",
        "/",
        "",
        "&ELECTRONS",
        "    conv_thr         = 1.0d-4",
        "    mixing_beta      = 0.01d0",
        "    electron_maxstep = 400",
        "    diagonalization  = 'cg'",
        "    mixing_mode      = 'local-TF'",
        "    mixing_ndim      = 8",
        "/",
        "",
        "&IONS",
        "    pot_extrapolation = 'second_order'",
        "    wfc_extrapolation = 'second_order'",
        "    ion_temperature   = 'svr'",
        f"    tempw             = {temperature:.1f}",
        "    nraise            = 20",
        "/",
        "",
        "ATOMIC_SPECIES",
        "Fe  55.845  Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "",
    ]
    pieces.extend(position_block)
    pieces.extend(
        [
            "",
            cell_header,
            *cell_lines,
            "",
            "K_POINTS automatic",
            f"{kx} {ky} {kz} 0 0 0",
            "",
        ]
    )
    return "\n".join(pieces)


def main() -> None:
    args = parse_args()
    data_dir = Path(__file__).resolve().parents[1]
    config_dir = resolve_dir(data_dir, args.config_dir)
    output_dir = config_dir / args.output_dir if not args.output_dir.is_absolute() else args.output_dir.resolve()
    pseudo_dir = args.pseudo_dir.resolve() if args.pseudo_dir is not None else repo_root()

    configs = sorted(path for path in config_dir.iterdir() if path.is_file() and CONFIG_RE.fullmatch(path.name))
    if not configs:
        raise FileNotFoundError(f"No qe_conf#### files found in {config_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    for config in configs:
        cfg_id = conf_id(config)
        nat, ntyp, cell_header, cell_lines, position_block = read_qe_conf(config)
        prefix = f"{args.prefix_base}_conf{cfg_id}"
        text = build_input_text(
            nat=nat,
            ntyp=ntyp,
            cell_header=cell_header,
            cell_lines=cell_lines,
            position_block=position_block,
            prefix=prefix,
            pseudo_dir=str(pseudo_dir),
            temperature=args.temperature,
            nstep=args.nstep,
            dt=args.dt,
            ecutwfc=args.ecutwfc,
            ecutrho=args.ecutrho,
            degauss=args.degauss,
            kgrid=tuple(args.kgrid),
        )
        output_path = output_dir / f"{prefix}.in"
        output_path.write_text(text)

    print(f"Source config dir: {config_dir}")
    print(f"Output QE input dir: {output_dir}")
    print(f"QE inputs written: {len(configs)}")
    print(f"First file: {args.prefix_base}_conf{conf_id(configs[0])}.in")
    print(f"Last file: {args.prefix_base}_conf{conf_id(configs[-1])}.in")


if __name__ == "__main__":
    main()
