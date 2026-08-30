#!/usr/bin/env python3
"""Generate the 16-species Round-17 Lonestar6 workflow from Round 15."""

from __future__ import annotations

import csv
import json
import math
import shutil
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
ROUND15 = ROOT.parent / "round15_qe75_patched_lambda0p100_md1000_lonestar6"
ROUND15_INPUT = (
    ROUND15 / "cases" / "a2p39_m2_lambda0p100_md1000" / "md_1000steps.in"
)
PSEUDO = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
LABELS = [f"Fe{i:02d}" for i in range(1, 17)]
CASE = "a2p39_2x2x2_m2_16species_lambda0p100_md1000"
NATOMS = 16
TARGET_TEMPERATURE = 4000.0
LABEL_SEED = 170016

KB_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
BOHR_M = 5.29177210903e-11
RYDBERG_ATOMIC_TIME_S = 2.0 * 2.4188843265864e-17
QE_VELOCITY_UNIT_M_PER_S = BOHR_M / RYDBERG_ATOMIC_TIME_S

# The first eight directions are exactly Round 15. The second eight are a
# common 47-degree azimuthal rotation. Each group contains four antipodal
# pairs and therefore has zero target-vector sum while all 16 directions are
# distinct.
ROUND15_ANGLES = (
    (72.8157590495, 72.7795960677),
    (73.5303944933, 92.3546343919),
    (33.6579186653, 331.8844103957),
    (95.5478199306, 0.9655109597),
    (107.1842409505, 252.7795960677),
    (106.4696055067, 272.3546343919),
    (146.3420813347, 151.8844103957),
    (84.4521800694, 180.9655109597),
)
ANGLES = ROUND15_ANGLES + tuple(
    (theta, (phi + 47.0) % 360.0) for theta, phi in ROUND15_ANGLES
)
STAGES = (
    ("stage1_4e-3", 0.005, 4.0e-3, 500),
    ("stage2_1e-3", 0.020, 1.0e-3, 500),
    ("stage3_3e-4", 0.100, 3.0e-4, 600),
)


def input_block(lines: list[str], header: str, count: int) -> list[str]:
    start = lines.index(header) + 1
    return lines[start : start + count]


def round15_state() -> tuple[np.ndarray, np.ndarray]:
    lines = ROUND15_INPUT.read_text().splitlines()
    positions = np.asarray(
        [[float(value) for value in row.split()[1:4]]
         for row in input_block(lines, "ATOMIC_POSITIONS crystal", 8)]
    )
    velocities = np.asarray(
        [[float(value) for value in row.split()[1:4]]
         for row in input_block(lines, "ATOMIC_VELOCITIES { a.u }", 8)]
    )
    first = positions.copy()
    first[:, 2] *= 0.5
    second = positions.copy()
    second[:, 2] = 0.5 * (second[:, 2] + 1.0)
    tiled_positions = np.concatenate((first, second), axis=0) % 1.0
    # Preserve the Round-15 speed distribution without making the upper tile
    # an exact dynamical clone of the lower one. A row permutation plus a
    # cyclic Cartesian-axis permutation is deterministic and norm preserving.
    upper_velocities = np.roll(velocities, 3, axis=0)[:, [1, 2, 0]]
    tiled_velocities = np.concatenate((velocities, upper_velocities), axis=0)
    return tiled_positions, tiled_velocities


def rescale_velocities(velocities: np.ndarray) -> tuple[np.ndarray, dict[str, object]]:
    mass_amu = np.full(NATOMS, 55.845)
    velocities = velocities.copy()
    velocities -= np.average(velocities, axis=0, weights=mass_amu)

    def temperature(values: np.ndarray) -> float:
        kinetic = 0.5 * np.sum(
            mass_amu[:, None] * AMU_SI
            * (values * QE_VELOCITY_UNIT_M_PER_S) ** 2
        )
        return 2.0 * kinetic / ((3 * NATOMS - 3) * KB_SI)

    before = temperature(velocities)
    velocities *= math.sqrt(TARGET_TEMPERATURE / before)
    after = temperature(velocities)
    com = np.average(velocities, axis=0, weights=mass_amu)
    return velocities, {
        "target_temperature_K": TARGET_TEMPERATURE,
        "temperature_before_rescaling_K": before,
        "measured_temperature_K": after,
        "center_of_mass_velocity_au": com.tolist(),
        "degrees_of_freedom": 3 * NATOMS - 3,
        "source": (
            "Round-15 eight-atom velocity set; upper tile uses a deterministic "
            "row and Cartesian-axis permutation; COM removed and exactly rescaled"
        ),
    }


def shuffled_labels() -> list[str]:
    rng = np.random.default_rng(LABEL_SEED)
    lower = np.asarray(LABELS[:8], dtype=object)
    upper = np.asarray(LABELS[8:], dtype=object)
    rng.shuffle(lower)
    rng.shuffle(upper)
    return np.concatenate((lower, upper)).tolist()


def direction(theta: float, phi: float) -> np.ndarray:
    theta_r, phi_r = np.radians((theta, phi))
    return np.asarray((
        np.sin(theta_r) * np.cos(phi_r),
        np.sin(theta_r) * np.sin(phi_r),
        np.cos(theta_r),
    ))


def system(lambda_value: float) -> list[str]:
    lines = [
        "&SYSTEM", "   ibrav=0", "   nat=16", "   ntyp=16",
        "   ecutwfc=71.0", "   ecutrho=496.0",
        "   occupations='smearing'", "   smearing='fd'",
        "   degauss=0.025334", "   nosym=.true.",
        "   noncolin=.true.", "   lspinorb=.false.",
        "   constrained_magnetization='atomic'",
        f"   lambda={lambda_value:.8f}", "   report=1",
    ]
    for index, (theta, phi) in enumerate(ANGLES, 1):
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


def control(calculation: str, prefix: str) -> list[str]:
    return [
        "&CONTROL", f"   calculation='{calculation}'",
        "   restart_mode='from_scratch'", f"   prefix='{prefix}'",
        "   pseudo_dir='../../pseudo'", "   outdir='./qe_tmp'",
        "   disk_io='low'", "   verbosity='high'",
        "   tstress=.true.", "   tprnfor=.true.",
    ]


def main() -> None:
    for name in ("cases", "logs", "pseudo"):
        (ROOT / name).mkdir(parents=True, exist_ok=True)
    source_pseudo = ROUND15 / "pseudo" / PSEUDO
    if not source_pseudo.is_file():
        raise FileNotFoundError(source_pseudo)
    shutil.copy2(source_pseudo, ROOT / "pseudo" / PSEUDO)

    positions, raw_velocities = round15_state()
    velocities, velocity_summary = rescale_velocities(raw_velocities)
    labels = shuffled_labels()
    vectors = np.asarray([direction(*angle) for angle in ANGLES])
    magnetic_summary = {
        "unique_species": 16,
        "unique_directions": 16,
        "one_species_per_atom": True,
        "target_local_moment_Bohr": 2.0,
        "first_eight_source": "Round 15",
        "second_eight_transform": "common +47 degree azimuthal rotation",
        "direction_vector_sum": vectors.sum(axis=0).tolist(),
        "target_moment_vector_sum_Bohr": (2.0 * vectors).sum(axis=0).tolist(),
        "position_label_seed": LABEL_SEED,
        "position_labels": labels,
    }
    structure_summary = {
        "source": str(ROUND15_INPUT.relative_to(ROOT.parent.parent.parent.parent)),
        "construction": "Round-15 2x2x1 structure tiled once along z",
        "cell_A": [4.78, 4.78, 4.78],
        "natoms": NATOMS,
    }
    (ROOT / "velocity_summary.json").write_text(json.dumps(velocity_summary, indent=2) + "\n")
    (ROOT / "magnetic_seed_summary.json").write_text(json.dumps(magnetic_summary, indent=2) + "\n")
    (ROOT / "structure_summary.json").write_text(json.dumps(structure_summary, indent=2) + "\n")
    (ROOT / "atomic_velocities_4000K.txt").write_text(
        "ATOMIC_VELOCITIES { a.u }\n" + "\n".join(
            f"{label} {velocity[0]: .16e} {velocity[1]: .16e} {velocity[2]: .16e}"
            for label, velocity in zip(labels, velocities)
        ) + "\n"
    )

    folder = ROOT / "cases" / CASE
    folder.mkdir(exist_ok=True)
    if list(folder.glob("*.out")) or (folder / "qe_tmp").exists():
        raise RuntimeError(f"Refusing to overwrite results in {folder}")
    prefix = "Fe_bcc16_r17_a2p39_m2_l0100_md1000"
    for stage_index, (stem, lambda_value, threshold, maxstep) in enumerate(STAGES):
        lines = control("scf", prefix) + ["/"]
        lines += system(lambda_value)
        lines += electrons(stage_index > 0, threshold, maxstep)
        lines += cards(positions, labels, None)
        (folder / f"{stem}.in").write_text("\n".join(lines))

    md_lines = control("md", prefix) + ["   dt=20.6700d0", "   nstep=1000", "/"]
    md_lines += system(0.100)
    md_lines += electrons(True, 4.0e-3, 500)
    md_lines += [
        "&IONS", "   ion_dynamics='verlet'", "   pot_extrapolation='atomic'",
        "   wfc_extrapolation='none'", "   ion_temperature='svr'",
        "   ion_velocities='from_input'", "   tempw=4000.0", "   nraise=20", "/",
    ]
    md_lines += cards(positions, labels, velocities)
    (folder / "md_1000steps.in").write_text("\n".join(md_lines))

    metadata = {
        "serial_index": 0,
        "case": CASE,
        "lattice_a_A": 2.39,
        "target_temperature_K": TARGET_TEMPERATURE,
        "target_local_moment_Bohr": 2.0,
        "starting_magnetization": 0.125,
        "lambda_stage1": 0.005,
        "lambda_stage2": 0.020,
        "lambda_stage3_md": 0.100,
        "md_steps": 1000,
        "qe_requirement": "QE 7.5 private patch mcons=zv*starting_magnetization and ntypx=128",
        "supercell": "2x2x2",
        "natoms": NATOMS,
        "ntyp": 16,
        "unique_spin_directions": 16,
        "one_species_per_atom": True,
        "nraise": 20,
        "mixing_mode": "local-TF",
        "mixing_beta": 0.010,
        "mixing_ndim": 20,
        "lonestar6_partition": "vm-small",
        "mpi_tasks": 16,
    }
    (folder / "case_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    with (ROOT / "case_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metadata), lineterminator="\n")
        writer.writeheader()
        writer.writerow(metadata)
    (ROOT / "cases.txt").write_text(CASE + "\n")
    print("Generated Round 17: 16 atoms, 16 species, 16 directions, MD1000")
    print("Velocity temperature:", velocity_summary["measured_temperature_K"])
    print("Target moment vector sum:", magnetic_summary["target_moment_vector_sum_Bohr"])


if __name__ == "__main__":
    main()
