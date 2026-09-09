#!/usr/bin/env python3
"""Generate a 32-atom Jakar pilot with the best Round-21 SCF settings."""

from __future__ import annotations

import csv
import importlib.util
import json
import shutil
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent
ROUND18 = PARENT / "round18_32atom_4x2x2_qe75_ntyp128_lambda0p100_md1000_lonestar6_vmsmall"
ROUND21 = PARENT / "round21_dlm01_scf_stabilization_jakar"
SOURCE = ROUND21 / "cases" / "conv1e-4_beta1e-2" / "md_100steps.in"
PSEUDO = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
CASE = "conv1e-4_beta1e-2_32atom"
PREFIX = "Fe_bcc32_r22_conv1e4_beta1e2"
LABEL_SEED = 220032
NATOMS = 32


def load_round18_helpers():
    path = ROUND18 / "generate_round18_lonestar6.py"
    spec = importlib.util.spec_from_file_location("round18_generator", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def input_block(lines: list[str], header: str, count: int) -> list[str]:
    start = lines.index(header) + 1
    return lines[start : start + count]


def round21_state() -> tuple[np.ndarray, np.ndarray]:
    lines = SOURCE.read_text().splitlines()
    positions = np.asarray([
        [float(value) for value in row.split()[1:4]]
        for row in input_block(lines, "ATOMIC_POSITIONS crystal", 16)
    ])
    velocities = np.asarray([
        [float(value) for value in row.split()[1:4]]
        for row in input_block(lines, "ATOMIC_VELOCITIES { a.u }", 16)
    ])
    left = positions.copy()
    left[:, 0] *= 0.5
    right = positions.copy()
    right[:, 0] = 0.5 * (right[:, 0] + 1.0)
    right_velocities = np.roll(velocities, 5, axis=0)[:, [2, 0, 1]]
    return (
        np.concatenate((left, right), axis=0) % 1.0,
        np.concatenate((velocities, right_velocities), axis=0),
    )


def main() -> None:
    helper = load_round18_helpers()
    for name in ("cases", "logs", "pseudo"):
        (ROOT / name).mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROUND21 / "pseudo" / PSEUDO, ROOT / "pseudo" / PSEUDO)

    positions, raw_velocities = round21_state()
    velocities, velocity_summary = helper.rescale_velocities(raw_velocities)
    velocity_summary["source"] = (
        "Round-21 16-atom difficult-seed velocities; the second tile uses a "
        "deterministic row and Cartesian-axis permutation; center-of-mass "
        "motion removed and the full 32-atom set rescaled to 4000 K"
    )
    rng = np.random.default_rng(LABEL_SEED)
    labels = np.asarray(helper.LABELS, dtype=object)
    rng.shuffle(labels)
    labels = labels.tolist()
    vectors = np.asarray([helper.direction(*angle) for angle in helper.ANGLES])

    summaries = {
        "structure_summary.json": {
            "source": str(SOURCE.relative_to(PARENT)),
            "construction": "Round-21 2x2x2 state tiled once along x",
            "supercell": "4x2x2 conventional BCC",
            "cell_A": [9.56, 4.78, 4.78],
            "natoms": NATOMS,
        },
        "velocity_summary.json": velocity_summary,
        "magnetic_seed_summary.json": {
            "unique_species": 32,
            "unique_directions": 32,
            "one_species_per_atom": True,
            "target_local_moment_Bohr": 2.0,
            "direction_vector_sum": vectors.sum(axis=0).tolist(),
            "target_moment_vector_sum_Bohr": (2.0 * vectors).sum(axis=0).tolist(),
            "position_label_seed": LABEL_SEED,
            "position_labels": labels,
        },
    }
    for filename, data in summaries.items():
        (ROOT / filename).write_text(json.dumps(data, indent=2) + "\n")

    folder = ROOT / "cases" / CASE
    folder.mkdir(exist_ok=True)
    if list(folder.glob("*.out")) or (folder / "qe_tmp").exists():
        raise RuntimeError(f"Refusing to overwrite results in {folder}")

    stages = (
        ("stage1_4e-3", 0.005, 4.0e-3, 500),
        ("stage2_1e-3", 0.020, 1.0e-3, 500),
        ("stage3_1e-4", 0.100, 1.0e-4, 1000),
    )
    for index, (stem, lambda_value, threshold, maxstep) in enumerate(stages):
        lines = helper.control("scf", PREFIX) + ["/"]
        lines += helper.system(lambda_value)
        lines += helper.electrons(index > 0, threshold, maxstep)
        lines += helper.cards(positions, labels, None)
        (folder / f"{stem}.in").write_text("\n".join(lines))

    lines = helper.control("md", PREFIX) + ["   dt=20.6700d0", "   nstep=100", "/"]
    lines += helper.system(0.100)
    lines += helper.electrons(True, 1.0e-4, 1000)
    lines += [
        "&IONS", "   ion_dynamics='verlet'", "   pot_extrapolation='atomic'",
        "   wfc_extrapolation='none'", "   ion_temperature='svr'",
        "   ion_velocities='from_input'", "   tempw=4000.0", "   nraise=5", "/",
    ]
    lines += helper.cards(positions, labels, velocities)
    (folder / "md_100steps.in").write_text("\n".join(lines))

    metadata = {
        "case": CASE, "source_round": 21, "supercell": "4x2x2",
        "natoms": 32, "ntyp": 32, "md_steps": 100,
        "md_conv_thr_Ry": 1.0e-4, "mixing_mode": "local-TF",
        "mixing_beta": 0.010, "mixing_ndim": 20,
        "electron_maxstep": 1000, "lambda": 0.100, "nraise": 5,
        "target_temperature_K": 4000.0, "mpi_tasks": 32,
        "qe_requirement": "QE 7.5 atomic-constraint fix with ntypx=128",
    }
    (folder / "case_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    with (ROOT / "case_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metadata), lineterminator="\n")
        writer.writeheader()
        writer.writerow(metadata)
    (ROOT / "cases.txt").write_text(CASE + "\n")
    print(f"Generated {CASE}: 32 atoms, 32 spin species, 100 MD steps")
    print(f"Velocity temperature: {velocity_summary['measured_temperature_K']:.6f} K")


if __name__ == "__main__":
    main()
