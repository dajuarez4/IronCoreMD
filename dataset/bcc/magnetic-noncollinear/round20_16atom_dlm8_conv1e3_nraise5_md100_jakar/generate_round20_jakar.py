#!/usr/bin/env python3
"""Generate the eight-replica Round-20 constrained-DLM Jakar workflow."""

from __future__ import annotations

import csv
import json
import math
import shutil
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
ROUND17 = ROOT.parent / "round17_16atom_2x2x2_qe75_ntyp128_lambda0p100_md1000_jakar"
ROUND17_INPUT = (
    ROUND17 / "cases" / "a2p39_2x2x2_m2_16species_lambda0p100_md1000"
    / "md_1000steps.in"
)
PSEUDO = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
LABELS = [f"Fe{i:02d}" for i in range(1, 17)]
NATOMS = 16
NREPLICAS = 8
TARGET_TEMPERATURE = 4000.0
SPIN_SEED_BASE = 200000
VELOCITY_SEED_BASE = 210000
LABEL_SEED_BASE = 220000

KB_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
BOHR_M = 5.29177210903e-11
RYDBERG_ATOMIC_TIME_S = 2.0 * 2.4188843265864e-17
QE_VELOCITY_UNIT_M_PER_S = BOHR_M / RYDBERG_ATOMIC_TIME_S

STAGES = (
    ("stage1_4e-3", 0.005, 4.0e-3, 500),
    ("stage2_1e-3", 0.020, 1.0e-3, 500),
    ("stage3_3e-4", 0.100, 3.0e-4, 600),
)


def input_block(lines: list[str], header: str, count: int) -> list[str]:
    start = lines.index(header) + 1
    return lines[start:start + count]


def initial_positions() -> np.ndarray:
    lines = ROUND17_INPUT.read_text().splitlines()
    return np.asarray([
        [float(value) for value in row.split()[1:4]]
        for row in input_block(lines, "ATOMIC_POSITIONS crystal", NATOMS)
    ])


def maxwell_velocities(seed: int) -> tuple[np.ndarray, dict[str, object]]:
    mass_amu = np.full(NATOMS, 55.845)
    sigma_si = np.sqrt(KB_SI * TARGET_TEMPERATURE / (mass_amu * AMU_SI))
    sigma_au = sigma_si / QE_VELOCITY_UNIT_M_PER_S
    values = np.random.default_rng(seed).normal(
        0.0, sigma_au[:, None], size=(NATOMS, 3)
    )
    values -= np.average(values, axis=0, weights=mass_amu)

    def temperature(velocities: np.ndarray) -> float:
        kinetic = 0.5 * np.sum(
            mass_amu[:, None] * AMU_SI
            * (velocities * QE_VELOCITY_UNIT_M_PER_S) ** 2
        )
        return 2.0 * kinetic / ((3 * NATOMS - 3) * KB_SI)

    before = temperature(values)
    values *= math.sqrt(TARGET_TEMPERATURE / before)
    return values, {
        "velocity_seed": seed,
        "target_temperature_K": TARGET_TEMPERATURE,
        "temperature_before_rescaling_K": before,
        "measured_temperature_K": temperature(values),
        "center_of_mass_velocity_au": np.average(
            values, axis=0, weights=mass_amu
        ).tolist(),
        "degrees_of_freedom": 3 * NATOMS - 3,
        "generator": "Maxwell-Boltzmann; COM removed; exactly rescaled",
    }


def balanced_angles(seed: int) -> tuple[list[tuple[float, float]], dict[str, object]]:
    """Generate eight random directions and their exact antipodes."""
    rng = np.random.default_rng(seed)
    directions: list[np.ndarray] = []
    for _ in range(NATOMS // 2):
        z = rng.uniform(-1.0, 1.0)
        phi = rng.uniform(0.0, 2.0 * np.pi)
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        vector = np.asarray((radius * math.cos(phi), radius * math.sin(phi), z))
        directions.extend((vector, -vector))
    order = rng.permutation(NATOMS)
    vectors = np.asarray(directions)[order]
    angles = [
        (
            math.degrees(math.acos(float(np.clip(vector[2], -1.0, 1.0)))),
            math.degrees(math.atan2(vector[1], vector[0])) % 360.0,
        )
        for vector in vectors
    ]
    cross_dot = vectors @ vectors.T
    np.fill_diagonal(cross_dot, -2.0)
    return angles, {
        "spin_seed": seed,
        "construction": "eight random unit vectors plus their antipodes, then shuffled",
        "unique_directions": len({(round(a, 10), round(b, 10)) for a, b in angles}),
        "direction_vector_sum": vectors.sum(axis=0).tolist(),
        "target_moment_vector_sum_Bohr": (2.0 * vectors).sum(axis=0).tolist(),
        "largest_dot_between_distinct_directions": float(cross_dot.max()),
        "angles_degrees": [[theta, phi] for theta, phi in angles],
    }


def shuffled_labels(seed: int) -> list[str]:
    labels = np.asarray(LABELS, dtype=object)
    np.random.default_rng(seed).shuffle(labels)
    return labels.tolist()


def system(lambda_value: float, angles: list[tuple[float, float]]) -> list[str]:
    lines = [
        "&SYSTEM", "   ibrav=0", "   nat=16", "   ntyp=16",
        "   ecutwfc=71.0", "   ecutrho=496.0",
        "   occupations='smearing'", "   smearing='fd'",
        "   degauss=0.025334", "   nosym=.true.",
        "   noncolin=.true.", "   lspinorb=.false.",
        "   constrained_magnetization='atomic'",
        f"   lambda={lambda_value:.8f}", "   report=1",
    ]
    for index, (theta, phi) in enumerate(angles, 1):
        lines.extend((
            f"   starting_magnetization({index})=0.12500000",
            f"   angle1({index})={theta:.10f}",
            f"   angle2({index})={phi:.10f}",
        ))
    return lines + ["/"]


def electrons(read_state: bool, threshold: float, maxstep: int) -> list[str]:
    return [
        "&ELECTRONS",
        f"   startingpot='{'file' if read_state else 'atomic'}'",
        f"   startingwfc='{'file' if read_state else 'atomic+random'}'",
        f"   conv_thr={threshold:.1e}", f"   electron_maxstep={maxstep}",
        "   scf_must_converge=.true.", "   mixing_mode='local-TF'",
        "   mixing_beta=0.010d0", "   mixing_ndim=20",
        "   diagonalization='cg'", "   diago_thr_init=1.0d-3",
        "   diago_cg_maxiter=200", "/",
    ]


def control(calculation: str, prefix: str) -> list[str]:
    return [
        "&CONTROL", f"   calculation='{calculation}'",
        "   restart_mode='from_scratch'", f"   prefix='{prefix}'",
        "   pseudo_dir='../../pseudo'", "   outdir='./qe_tmp'",
        "   disk_io='low'", "   verbosity='high'",
        "   tstress=.true.", "   tprnfor=.true.",
    ]


def cards(
    positions: np.ndarray, labels: list[str], velocities: np.ndarray | None
) -> list[str]:
    lines = ["ATOMIC_SPECIES"]
    lines.extend(f"{label} 55.845 {PSEUDO}" for label in LABELS)
    lines.append("ATOMIC_POSITIONS crystal")
    lines.extend(
        f"{label} {position[0]:.12f} {position[1]:.12f} {position[2]:.12f}"
        for label, position in zip(labels, positions)
    )
    if velocities is not None:
        lines.append("ATOMIC_VELOCITIES { a.u }")
        lines.extend(
            f"{label} {velocity[0]: .16e} {velocity[1]: .16e} {velocity[2]: .16e}"
            for label, velocity in zip(labels, velocities)
        )
    return lines + [
        "CELL_PARAMETERS angstrom",
        "4.780000000000 0.000000000000 0.000000000000",
        "0.000000000000 4.780000000000 0.000000000000",
        "0.000000000000 0.000000000000 4.780000000000",
        "K_POINTS automatic", "1 1 1 0 0 0", "",
    ]


def main() -> None:
    for name in ("cases", "logs", "pseudo"):
        (ROOT / name).mkdir(parents=True, exist_ok=True)
    source_pseudo = ROUND17 / "pseudo" / PSEUDO
    if not source_pseudo.is_file():
        raise FileNotFoundError(source_pseudo)
    shutil.copy2(source_pseudo, ROOT / "pseudo" / PSEUDO)

    positions = initial_positions()
    manifest: list[dict[str, object]] = []
    ensemble: list[dict[str, object]] = []
    case_names: list[str] = []

    for replica in range(NREPLICAS):
        spin_seed = SPIN_SEED_BASE + replica
        velocity_seed = VELOCITY_SEED_BASE + replica
        label_seed = LABEL_SEED_BASE + replica
        case = f"dlm{replica:02d}_spin{spin_seed}_vel{velocity_seed}"
        case_names.append(case)
        folder = ROOT / "cases" / case
        folder.mkdir(exist_ok=True)
        if list(folder.glob("*.out")) or (folder / "qe_tmp").exists():
            raise RuntimeError(f"Refusing to overwrite results in {folder}")

        angles, magnetic_summary = balanced_angles(spin_seed)
        velocities, velocity_summary = maxwell_velocities(velocity_seed)
        labels = shuffled_labels(label_seed)
        prefix = f"Fe_bcc16_r20_dlm{replica:02d}"

        for stage_index, (stem, lambda_value, threshold, maxstep) in enumerate(STAGES):
            lines = control("scf", prefix) + ["/"]
            lines += system(lambda_value, angles)
            lines += electrons(stage_index > 0, threshold, maxstep)
            lines += cards(positions, labels, None)
            (folder / f"{stem}.in").write_text("\n".join(lines))

        md_lines = control("md", prefix) + [
            "   dt=20.6700d0", "   nstep=100", "/",
        ]
        md_lines += system(0.100, angles)
        md_lines += electrons(True, 1.0e-3, 500)
        md_lines += [
            "&IONS", "   ion_dynamics='verlet'", "   pot_extrapolation='atomic'",
            "   wfc_extrapolation='none'", "   ion_temperature='svr'",
            "   ion_velocities='from_input'", "   tempw=4000.0",
            "   nraise=5", "/",
        ]
        md_lines += cards(positions, labels, velocities)
        (folder / "md_100steps.in").write_text("\n".join(md_lines))

        metadata = {
            "serial_index": replica,
            "case": case,
            "spin_seed": spin_seed,
            "velocity_seed": velocity_seed,
            "label_seed": label_seed,
            "target_temperature_K": TARGET_TEMPERATURE,
            "target_local_moment_Bohr": 2.0,
            "lambda_stage1": 0.005,
            "lambda_stage2": 0.020,
            "lambda_stage3_md": 0.100,
            "md_steps": 100,
            "md_conv_thr_Ry": 1.0e-3,
            "nraise": 5,
            "natoms": NATOMS,
            "ntyp": 16,
            "unique_spin_directions": 16,
            "one_species_per_atom": True,
            "supercell": "2x2x2",
            "mixing_mode": "local-TF",
            "mixing_beta": 0.010,
            "mixing_ndim": 20,
            "jakar_account": "jakar_general",
            "jakar_qos": "jakar_medium_general",
            "jakar_partition": "medium",
            "mpi_tasks": 16,
        }
        (folder / "case_metadata.json").write_text(
            json.dumps(metadata, indent=2) + "\n"
        )
        manifest.append(metadata)
        ensemble.append({
            "case": case,
            "position_labels": labels,
            "magnetic": magnetic_summary,
            "velocity": velocity_summary,
        })

    with (ROOT / "case_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    (ROOT / "cases.txt").write_text("\n".join(case_names) + "\n")
    (ROOT / "ensemble_summary.json").write_text(
        json.dumps(ensemble, indent=2) + "\n"
    )
    (ROOT / "structure_summary.json").write_text(json.dumps({
        "source": str(ROUND17_INPUT.relative_to(ROOT.parent.parent.parent.parent)),
        "cell_A": [4.78, 4.78, 4.78],
        "supercell": "2x2x2",
        "natoms": NATOMS,
        "same_initial_positions_for_all_replicas": True,
    }, indent=2) + "\n")
    print(f"Generated Round 20: {NREPLICAS} balanced-DLM replicas, MD100 each")
    print("Pilot case:", case_names[0])
    print("Ensemble cases:", ", ".join(case_names[1:]))


if __name__ == "__main__":
    main()
