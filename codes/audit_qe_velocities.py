#!/usr/bin/env python3
"""Audit a Quantum ESPRESSO ATOMIC_VELOCITIES card under both atomic-unit conventions."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
import sys

import numpy as np

KB_J_PER_K = 1.380649e-23
AMU_KG = 1.66053906660e-27
BOHR_M = 5.29177210903e-11
HARTREE_ATOMIC_TIME_S = 2.4188843265864e-17
RYDBERG_ATOMIC_TIME_S = 4.8377686531728e-17
HARTREE_VELOCITY_UNIT_M_PER_S = BOHR_M / HARTREE_ATOMIC_TIME_S
QE_VELOCITY_UNIT_M_PER_S = BOHR_M / RYDBERG_ATOMIC_TIME_S
CARD_NAMES = (
    "&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL", "ATOMIC_SPECIES",
    "ATOMIC_POSITIONS", "ATOMIC_VELOCITIES", "K_POINTS", "CELL_PARAMETERS",
    "CONSTRAINTS", "OCCUPATIONS", "ATOMIC_FORCES", "SOLVENTS", "HUBBARD",
)


def card_rows(lines: list[str], name: str) -> list[list[str]]:
    start = next(index for index, line in enumerate(lines) if line.strip().upper().startswith(name))
    rows: list[list[str]] = []
    for line in lines[start + 1 :]:
        stripped = line.strip()
        if stripped and any(stripped.upper().startswith(card) for card in CARD_NAMES):
            break
        if stripped and not stripped.startswith(("!", "#")):
            rows.append(stripped.replace("D", "E").replace("d", "e").split())
    return rows


def parse_qe_input(path: Path) -> tuple[np.ndarray, np.ndarray, float | None]:
    text = path.read_text()
    lines = text.splitlines()
    masses = {row[0]: float(row[1]) for row in card_rows(lines, "ATOMIC_SPECIES") if len(row) >= 3}
    velocity_rows = [row for row in card_rows(lines, "ATOMIC_VELOCITIES") if len(row) >= 4]
    labels = [row[0] for row in velocity_rows]
    velocities = np.asarray([[float(value) for value in row[1:4]] for row in velocity_rows])
    atom_masses = np.asarray([masses[label] for label in labels])
    match = re.search(r"\btempw\s*=\s*([+\-0-9.eEdD]+)", text, flags=re.IGNORECASE)
    tempw = float(match.group(1).replace("D", "E").replace("d", "e")) if match else None
    return atom_masses, velocities, tempw


def kinetic_energy_j(masses_amu: np.ndarray, velocities_numeric: np.ndarray, unit_m_per_s: float) -> float:
    velocities_si = velocities_numeric * unit_m_per_s
    return float(0.5 * np.sum(masses_amu[:, None] * AMU_KG * velocities_si**2))


def temperature_k(kinetic_energy_joule: float, degrees_of_freedom: int) -> float:
    return 2.0 * kinetic_energy_joule / (degrees_of_freedom * KB_J_PER_K)


def audit(path: Path) -> dict[str, float | int | None | np.ndarray]:
    masses, velocities, tempw = parse_qe_input(path)
    com = np.average(velocities, axis=0, weights=masses)
    centered = velocities - com
    result: dict[str, float | int | None | np.ndarray] = {
        "atoms": len(masses), "tempw_k": tempw, "com_numeric": com,
    }
    for name, unit in (("hartree", HARTREE_VELOCITY_UNIT_M_PER_S), ("qe_rydberg", QE_VELOCITY_UNIT_M_PER_S)):
        kinetic_raw = kinetic_energy_j(masses, velocities, unit)
        kinetic_com = kinetic_energy_j(masses, centered, unit)
        result[f"{name}_kinetic_raw_j"] = kinetic_raw
        result[f"{name}_kinetic_com_j"] = kinetic_com
        result[f"{name}_temperature_3n_k"] = temperature_k(kinetic_raw, 3 * len(masses))
        result[f"{name}_temperature_3n_minus_3_k"] = temperature_k(kinetic_com, 3 * len(masses) - 3)
    qe_temperature = float(result["qe_rydberg_temperature_3n_minus_3_k"])
    result["scale_to_tempw"] = None if tempw is None else (tempw / qe_temperature) ** 0.5
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("qe_input", type=Path)
    args = parser.parse_args()
    result = audit(args.qe_input)
    print(f"Input: {args.qe_input}")
    print(f"Atoms: {result['atoms']}")
    print(f"tempw: {result['tempw_k']} K")
    print(f"COM velocity (numeric card units): {np.asarray(result['com_numeric'])}")
    for name in ("hartree", "qe_rydberg"):
        print(f"\n{name} interpretation:")
        print(f"  K before COM removal: {result[f'{name}_kinetic_raw_j']:.16e} J")
        print(f"  K after COM removal:  {result[f'{name}_kinetic_com_j']:.16e} J")
        print(f"  T using 3N:           {result[f'{name}_temperature_3n_k']:.9f} K")
        print(f"  T using 3N-3:         {result[f'{name}_temperature_3n_minus_3_k']:.9f} K")
    print(f"\nScale factor needed for tempw under QE units: {result['scale_to_tempw']}")


if __name__ == "__main__":
    main()
