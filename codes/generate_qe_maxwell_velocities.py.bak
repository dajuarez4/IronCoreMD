#!/usr/bin/env python3
"""Generate Maxwell-Boltzmann ionic velocities for Quantum ESPRESSO inputs.

The script reads ATOMIC_SPECIES and ATOMIC_POSITIONS from a QE pw.x input,
samples independent Cartesian velocities from the Maxwell-Boltzmann
distribution at a target temperature, removes center-of-mass drift, optionally
rescales the result to hit the target temperature exactly, and writes:

1. An ``ATOMIC_VELOCITIES { a.u }`` block.
2. Optionally, a patched QE input with that block inserted or replaced.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

KB_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
AU_VELOCITY_SI = 2.18769126364e6

CARD_PREFIXES = (
    "&CONTROL",
    "&SYSTEM",
    "&ELECTRONS",
    "&IONS",
    "&CELL",
    "ATOMIC_SPECIES",
    "ATOMIC_POSITIONS",
    "K_POINTS",
    "CELL_PARAMETERS",
    "CONSTRAINTS",
    "OCCUPATIONS",
    "ATOMIC_VELOCITIES",
    "ATOMIC_FORCES",
    "SOLVENTS",
    "HUBBARD",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Maxwell-Boltzmann ATOMIC_VELOCITIES for a QE pw.x input."
    )
    parser.add_argument("input_qe", type=Path, help="Input QE file to read.")
    parser.add_argument(
        "--temperature",
        type=float,
        default=4500.0,
        help="Target temperature in K. Default: 4500.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=4500,
        help="Random seed used for reproducible velocities. Default: 4500.",
    )
    parser.add_argument(
        "--output-block",
        type=Path,
        default=Path("atomic_velocities_4500K.txt"),
        help="Output file for the ATOMIC_VELOCITIES card.",
    )
    parser.add_argument(
        "--patched-input",
        type=Path,
        default=None,
        help="Optional QE input written with ATOMIC_VELOCITIES inserted.",
    )
    parser.add_argument(
        "--keep-com-drift",
        action="store_true",
        help="Do not remove center-of-mass drift.",
    )
    parser.add_argument(
        "--no-rescale-exact",
        action="store_true",
        help="Do not rescale velocities to match the target temperature exactly after COM removal.",
    )
    return parser.parse_args()


def is_card_start(line: str) -> bool:
    stripped = line.strip()
    if not stripped:
        return False
    upper = stripped.upper()
    return any(upper.startswith(prefix) for prefix in CARD_PREFIXES)


def find_card_start(lines: list[str], card_name: str) -> int:
    target = card_name.upper()
    for index, line in enumerate(lines):
        if line.strip().upper().startswith(target):
            return index
    raise ValueError(f"Could not find {card_name} in QE input.")


def card_end_index(lines: list[str], start_index: int) -> int:
    index = start_index + 1
    while index < len(lines):
        stripped = lines[index].strip()
        if stripped and is_card_start(lines[index]):
            break
        index += 1
    return index


def parse_atomic_species(lines: list[str]) -> dict[str, float]:
    start = find_card_start(lines, "ATOMIC_SPECIES")
    end = card_end_index(lines, start)
    masses: dict[str, float] = {}
    for line in lines[start + 1 : end]:
        stripped = line.strip()
        if not stripped or stripped.startswith(("!", "#")):
            continue
        fields = stripped.split()
        if len(fields) < 3:
            continue
        masses[fields[0]] = float(fields[1])
    if not masses:
        raise ValueError("No ATOMIC_SPECIES entries found.")
    return masses


def parse_atomic_positions(lines: list[str]) -> list[str]:
    start = find_card_start(lines, "ATOMIC_POSITIONS")
    end = card_end_index(lines, start)
    labels: list[str] = []
    for line in lines[start + 1 : end]:
        stripped = line.strip()
        if not stripped or stripped.startswith(("!", "#")):
            continue
        fields = stripped.split()
        if len(fields) < 4:
            continue
        labels.append(fields[0])
    if not labels:
        raise ValueError("No ATOMIC_POSITIONS entries found.")
    return labels


def read_qe_species_sequence(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text().splitlines()
    mass_map = parse_atomic_species(lines)
    labels = parse_atomic_positions(lines)
    masses = np.array([mass_map[label] for label in labels], dtype=float)
    return labels, masses


def sigma_au(mass_amu: np.ndarray, temperature_k: float) -> np.ndarray:
    sigma_si = np.sqrt(KB_SI * temperature_k / (mass_amu * AMU_SI))
    return sigma_si / AU_VELOCITY_SI


def sample_maxwell_velocities_au(
    masses_amu: np.ndarray,
    temperature_k: float,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    sigma = sigma_au(masses_amu, temperature_k)
    return rng.normal(loc=0.0, scale=sigma[:, None], size=(len(masses_amu), 3))


def remove_center_of_mass_drift_au(velocities_au: np.ndarray, masses_amu: np.ndarray) -> np.ndarray:
    weighted_mean = np.average(velocities_au, axis=0, weights=masses_amu)
    return velocities_au - weighted_mean


def temperature_from_velocities_au(
    velocities_au: np.ndarray,
    masses_amu: np.ndarray,
    remove_com: bool = True,
) -> float:
    kinetic_joule = 0.5 * np.sum((masses_amu[:, None] * AMU_SI) * (velocities_au * AU_VELOCITY_SI) ** 2)
    dof = 3 * len(masses_amu) - (3 if remove_com else 0)
    if dof <= 0:
        raise ValueError("Non-positive number of degrees of freedom.")
    return 2.0 * kinetic_joule / (dof * KB_SI)


def rescale_to_temperature_au(
    velocities_au: np.ndarray,
    masses_amu: np.ndarray,
    target_temperature_k: float,
    remove_com: bool = True,
) -> np.ndarray:
    current_temperature = temperature_from_velocities_au(velocities_au, masses_amu, remove_com=remove_com)
    if current_temperature <= 0.0:
        raise ValueError("Current velocity set has non-positive temperature.")
    return velocities_au * math.sqrt(target_temperature_k / current_temperature)


def format_atomic_velocities_card(labels: list[str], velocities_au: np.ndarray) -> str:
    lines = ["ATOMIC_VELOCITIES { a.u }"]
    for label, velocity in zip(labels, velocities_au):
        lines.append(f"{label:4s} {velocity[0]: .16e} {velocity[1]: .16e} {velocity[2]: .16e}")
    return "\n".join(lines) + "\n"


def replace_or_append_atomic_velocities(qe_text: str, velocity_block: str) -> str:
    lines = qe_text.splitlines()
    try:
        start = find_card_start(lines, "ATOMIC_VELOCITIES")
    except ValueError:
        start = -1

    if start >= 0:
        end = card_end_index(lines, start)
        new_lines = lines[:start] + velocity_block.rstrip("\n").splitlines() + lines[end:]
        return "\n".join(new_lines) + "\n"

    positions_start = find_card_start(lines, "ATOMIC_POSITIONS")
    insert_at = card_end_index(lines, positions_start)
    for index in range(insert_at, len(lines)):
        if lines[index].strip().upper().startswith(("K_POINTS", "CELL_PARAMETERS", "CONSTRAINTS", "OCCUPATIONS")):
            insert_at = index
            break
    new_lines = lines[:insert_at] + [""] + velocity_block.rstrip("\n").splitlines() + [""] + lines[insert_at:]
    return "\n".join(new_lines) + "\n"


def main() -> None:
    args = parse_args()
    labels, masses_amu = read_qe_species_sequence(args.input_qe)
    velocities = sample_maxwell_velocities_au(masses_amu, args.temperature, args.seed)

    if not args.keep_com_drift:
        velocities = remove_center_of_mass_drift_au(velocities, masses_amu)

    if not args.no_rescale_exact:
        velocities = rescale_to_temperature_au(
            velocities,
            masses_amu,
            args.temperature,
            remove_com=not args.keep_com_drift,
        )

    measured_temperature = temperature_from_velocities_au(
        velocities,
        masses_amu,
        remove_com=not args.keep_com_drift,
    )
    velocity_block = format_atomic_velocities_card(labels, velocities)

    output_block = args.output_block if args.output_block.is_absolute() else args.input_qe.parent / args.output_block
    output_block.write_text(velocity_block)

    print(f"Input QE file: {args.input_qe}")
    print(f"Atoms: {len(labels)}")
    print(f"Target temperature: {args.temperature:.2f} K")
    print(f"Measured temperature after processing: {measured_temperature:.6f} K")
    print(f"Removed COM drift: {not args.keep_com_drift}")
    print(f"Rescaled to exact target temperature: {not args.no_rescale_exact}")
    print(f"Wrote velocity card: {output_block}")

    if args.patched_input is not None:
        patched_path = args.patched_input if args.patched_input.is_absolute() else args.input_qe.parent / args.patched_input
        patched_path.write_text(replace_or_append_atomic_velocities(args.input_qe.read_text(), velocity_block))
        print(f"Wrote patched QE input: {patched_path}")


if __name__ == "__main__":
    main()
