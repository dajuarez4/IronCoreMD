#!/usr/bin/env python3
"""Summarize TDEP Helmholtz free energies for one phase."""

from __future__ import annotations

import argparse
from pathlib import Path

from tdep_common import classify_free_energy, default_dataset_dir, read_free_energy, read_u0_second_order
from tdep_phases import PHASE_SPECS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize TDEP free-energy outputs for one phase.")
    parser.add_argument(
        "folders",
        nargs="*",
        type=Path,
        help="TDEP folders to summarize. Default: all tdep_* folders inside dataset-dir.",
    )
    parser.add_argument("--phase", choices=sorted(PHASE_SPECS), default="bcc", help="Crystal phase. Default: bcc.")
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=None,
        help="Directory containing the TDEP folders. Default: <repo>/dataset/<phase>.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    dataset_dir = args.dataset_dir.resolve() if args.dataset_dir is not None else default_dataset_dir(args.phase).resolve()
    folders = [path.resolve() if path.is_absolute() else (dataset_dir / path).resolve() for path in args.folders]
    if not folders:
        folders = sorted(path.resolve() for path in dataset_dir.glob("tdep_*") if path.is_dir())

    print("folder,T_K,F_vib_eV_atom,U0_2nd_eV_atom,F_total_eV_atom,S_eV_K_atom,Cv_eV_K_atom,status")
    for folder in folders:
        temperature, f_vib, entropy, cv = read_free_energy(folder / "outfile.free_energy")
        u0 = read_u0_second_order(folder / "outfile.U0")
        total = u0 + f_vib
        status = classify_free_energy(temperature, f_vib, entropy, cv, u0)
        print(
            f"{folder.name},{temperature:.5f},{f_vib:.12f},{u0:.12f},"
            f"{total:.12f},{entropy:.12e},{cv:.12e},{status}"
        )


if __name__ == "__main__":
    main()
