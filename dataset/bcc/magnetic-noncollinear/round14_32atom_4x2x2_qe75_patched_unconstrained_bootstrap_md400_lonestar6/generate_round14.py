#!/usr/bin/env python3
"""Generate the unconstrained-bootstrap 32-atom noncollinear Round-14 workflow."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
ROUND12 = ROOT.parent / "round12_32atom_4x2x2_qe75_patched_lambda0p100_md400_lonestar6"
THERMAL_SOURCE = ROOT.parent.parent / "non-mag" / "2.39_4000K.npz"
PSEUDO = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
LABELS = [f"Fe{i:02d}" for i in range(1, 9)]
NATOMS = 32
LATTICE_A = 2.39
TARGET_TEMPERATURE = 4000.0
THERMAL_FRAME_INDEX = 27
DISPLACEMENT_SEED = 320889
VELOCITY_SEED = 320400
LABEL_SEED = 130032
KB_SI = 1.380649e-23
AMU_SI = 1.66053906660e-27
BOHR_M = 5.29177210903e-11
RYDBERG_ATOMIC_TIME_S = 2.0 * 2.4188843265864e-17
QE_VELOCITY_UNIT_M_PER_S = BOHR_M / RYDBERG_ATOMIC_TIME_S

# Four antipodal pairs. Repeating every direction four times gives an exactly
# zero vector sum while staying below QE's ten-species limit.
ANGLES = (
    (72.8157590495, 72.7795960677),
    (73.5303944933, 92.3546343919),
    (33.6579186653, 331.8844103957),
    (95.5478199306, 0.9655109597),
    (107.1842409505, 252.7795960677),
    (106.4696055067, 272.3546343919),
    (146.3420813347, 151.8844103957),
    (84.4521800694, 180.9655109597),
)
STAGES = (
    # name, lambda (Ry or None), conv_thr, maxstep, degauss, mixing_beta
    ("bootstrap_unconstrained_2e-2", None, 2.0e-2, 500, 0.080000, 0.050),
    ("bootstrap_l0p001_1e-2", 0.001, 1.0e-2, 500, 0.050000, 0.010),
    ("stage1_l0p005_4e-3", 0.005, 4.0e-3, 500, 0.050000, 0.010),
    ("stage2_l0p020_1e-3", 0.020, 1.0e-3, 500, 0.025334, 0.010),
    ("stage3_l0p100_3e-4", 0.100, 3.0e-4, 600, 0.025334, 0.010),
)
CASE = "a2p39_4x2x2_m2_unconstrained_bootstrap_lambda0p100_md400"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ideal_positions_and_sublattices() -> tuple[np.ndarray, np.ndarray]:
    positions: list[tuple[float, float, float]] = []
    sublattices: list[int] = []
    # Keep the two BCC sublattices in separate blocks for deterministic
    # displacement assignment, but do not tie magnetic labels to this order.
    for k in range(2):
        for j in range(2):
            for i in range(4):
                positions.append((i / 4, j / 2, k / 2))
                sublattices.append(0)
    for k in range(2):
        for j in range(2):
            for i in range(4):
                positions.append(((i + 0.5) / 4, (j + 0.5) / 2, (k + 0.5) / 2))
                sublattices.append(1)
    return np.asarray(positions, dtype=float), np.asarray(sublattices, dtype=int)


def spatially_mixed_labels(sublattices: np.ndarray) -> tuple[list[str], dict[str, object]]:
    """Put two copies of every direction on each BCC sublattice."""
    rng = np.random.default_rng(LABEL_SEED)
    labels = np.empty(NATOMS, dtype=object)
    per_sublattice: dict[str, dict[str, int]] = {}
    base = np.repeat(np.asarray(LABELS, dtype=object), 2)
    for sublattice in (0, 1):
        sites = np.flatnonzero(sublattices == sublattice)
        assigned = base.copy()
        rng.shuffle(assigned)
        labels[sites] = assigned
        per_sublattice[str(sublattice)] = {
            label: int(np.count_nonzero(assigned == label)) for label in LABELS
        }
    return labels.tolist(), {
        "label_seed": LABEL_SEED,
        "strategy": "independent deterministic shuffle; two copies of every direction on each BCC sublattice",
        "counts_per_sublattice": per_sublattice,
        "site_labels": labels.tolist(),
    }


def minimum_pair_distance(positions: np.ndarray, cell_lengths: np.ndarray) -> float:
    smallest = math.inf
    for index in range(len(positions) - 1):
        delta = positions[index + 1 :] - positions[index]
        delta -= np.rint(delta)
        smallest = min(smallest, float(np.linalg.norm(delta * cell_lengths, axis=1).min()))
    return smallest


def thermalized_positions() -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    if not THERMAL_SOURCE.is_file():
        raise FileNotFoundError(THERMAL_SOURCE)
    ideal, sublattices = ideal_positions_and_sublattices()
    with np.load(THERMAL_SOURCE, allow_pickle=False) as archive:
        source_positions = archive["positions"][THERMAL_FRAME_INDEX].astype(float)
        source_ideal = archive["initial_positions_alat"].astype(float)
        source_temperature = float(archive["temperature_K"][THERMAL_FRAME_INDEX])
        source_iteration = int(archive["iteration"][THERMAL_FRAME_INDEX])
        source_cell = np.diag(archive["input_cell_parameters"].astype(float))
    displacement_fractional = source_positions - source_ideal
    displacement_fractional -= np.rint(displacement_fractional)
    displacement_angstrom = displacement_fractional * source_cell
    selected = np.random.default_rng(DISPLACEMENT_SEED).choice(
        len(displacement_angstrom), NATOMS, replace=False
    )
    cell_lengths = np.asarray((4 * LATTICE_A, 2 * LATTICE_A, 2 * LATTICE_A))
    positions = (ideal + displacement_angstrom[selected] / cell_lengths) % 1.0
    applied = displacement_angstrom[selected]
    summary = {
        "source": str(THERMAL_SOURCE.relative_to(ROOT.parent.parent.parent)),
        "source_sha256": sha256(THERMAL_SOURCE),
        "source_frame_zero_based": THERMAL_FRAME_INDEX,
        "source_md_iteration": source_iteration,
        "source_temperature_K": source_temperature,
        "displacement_selection_seed": DISPLACEMENT_SEED,
        "selected_source_atom_indices_zero_based": selected.tolist(),
        "rms_applied_displacement_A": float(np.sqrt(np.mean(np.sum(applied**2, axis=1)))),
        "maximum_applied_displacement_A": float(np.linalg.norm(applied, axis=1).max()),
        "minimum_initial_pair_distance_A": minimum_pair_distance(positions, cell_lengths),
    }
    if summary["minimum_initial_pair_distance_A"] < 1.70:
        raise RuntimeError("Generated thermal structure contains an unexpectedly short pair distance")
    return positions, sublattices, summary


def maxwell_velocities() -> tuple[np.ndarray, dict[str, object]]:
    mass_amu = np.full(NATOMS, 55.845)
    sigma_si = np.sqrt(KB_SI * TARGET_TEMPERATURE / (mass_amu * AMU_SI))
    sigma_au = sigma_si / QE_VELOCITY_UNIT_M_PER_S
    velocities = np.random.default_rng(VELOCITY_SEED).normal(
        0.0, sigma_au[:, None], size=(NATOMS, 3)
    )
    velocities -= np.average(velocities, axis=0, weights=mass_amu)
    kinetic_joule = 0.5 * np.sum(
        (mass_amu[:, None] * AMU_SI)
        * (velocities * QE_VELOCITY_UNIT_M_PER_S) ** 2
    )
    measured = 2.0 * kinetic_joule / ((3 * NATOMS - 3) * KB_SI)
    velocities *= math.sqrt(TARGET_TEMPERATURE / measured)
    com = np.average(velocities, axis=0, weights=mass_amu)
    kinetic_joule = 0.5 * np.sum(
        (mass_amu[:, None] * AMU_SI)
        * (velocities * QE_VELOCITY_UNIT_M_PER_S) ** 2
    )
    measured = 2.0 * kinetic_joule / ((3 * NATOMS - 3) * KB_SI)
    return velocities, {
        "target_temperature_K": TARGET_TEMPERATURE,
        "measured_temperature_K": measured,
        "velocity_seed": VELOCITY_SEED,
        "center_of_mass_velocity_au": com.tolist(),
        "degrees_of_freedom": 3 * NATOMS - 3,
        "generator": "Maxwell-Boltzmann; COM removed; exactly rescaled",
    }


def system(lambda_value: float | None, degauss: float) -> list[str]:
    lines = [
        "&SYSTEM", "   ibrav=0", f"   nat={NATOMS}", "   ntyp=8",
        "   ecutwfc=71.0", "   ecutrho=496.0",
        "   occupations='smearing'", "   smearing='fd'",
        f"   degauss={degauss:.6f}", "   nosym=.true.",
        "   noncolin=.true.", "   lspinorb=.false.",
    ]
    if lambda_value is None:
        lines.extend(("   constrained_magnetization='none'", "   report=1"))
    else:
        lines.extend(("   constrained_magnetization='atomic'",
                      f"   lambda={lambda_value:.8f}", "   report=1"))
    for index, (theta, phi) in enumerate(ANGLES, 1):
        lines.extend((
            f"   starting_magnetization({index})=0.12500000",
            f"   angle1({index})={theta:.10f}",
            f"   angle2({index})={phi:.10f}",
        ))
    return lines + ["/"]


def electrons(
    read_state: bool, threshold: float, maxstep: int, mixing_beta: float
) -> list[str]:
    return [
        "&ELECTRONS",
        f"   startingpot='{'file' if read_state else 'atomic'}'",
        f"   startingwfc='{'file' if read_state else 'atomic+random'}'",
        f"   conv_thr={threshold:.1e}", f"   electron_maxstep={maxstep}",
        "   scf_must_converge=.true.", "   mixing_mode='plain'",
        f"   mixing_beta={mixing_beta:.3f}d0", "   mixing_ndim=40",
        "   diagonalization='cg'", "   diago_thr_init=1.0d-3",
        "   diago_cg_maxiter=200", "/",
    ]


def cards(
    positions: np.ndarray,
    labels: list[str],
    velocities: np.ndarray | None,
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
    lines.extend((
        "CELL_PARAMETERS angstrom",
        "9.560000000000 0.000000000000 0.000000000000",
        "0.000000000000 4.780000000000 0.000000000000",
        "0.000000000000 0.000000000000 4.780000000000",
        "K_POINTS automatic", "1 1 1 0 0 0", "",
    ))
    return lines


def control(calculation: str, prefix: str) -> list[str]:
    return [
        "&CONTROL", f"   calculation='{calculation}'", "   restart_mode='from_scratch'",
        f"   prefix='{prefix}'", "   pseudo_dir='../../pseudo'",
        "   outdir='./qe_tmp'", "   disk_io='low'", "   verbosity='high'",
        "   tstress=.true.", "   tprnfor=.true.",
    ]


def main() -> None:
    for name in ("cases", "logs", "pseudo"):
        (ROOT / name).mkdir(parents=True, exist_ok=True)
    source_pseudo = ROUND12 / "pseudo" / PSEUDO
    if not source_pseudo.is_file():
        raise FileNotFoundError(source_pseudo)
    shutil.copy2(source_pseudo, ROOT / "pseudo" / PSEUDO)

    positions, sublattices, thermal_summary = thermalized_positions()
    labels, label_summary = spatially_mixed_labels(sublattices)
    velocities, velocity_summary = maxwell_velocities()
    (ROOT / "thermalized_start_summary.json").write_text(
        json.dumps(thermal_summary, indent=2) + "\n"
    )
    (ROOT / "velocity_summary.json").write_text(
        json.dumps(velocity_summary, indent=2) + "\n"
    )
    directions = np.asarray([
        (math.sin(math.radians(theta)) * math.cos(math.radians(phi)),
         math.sin(math.radians(theta)) * math.sin(math.radians(phi)),
         math.cos(math.radians(theta)))
        for theta, phi in ANGLES
    ])
    magnetic_summary = {
        "unique_directions": len(ANGLES),
        "copies_per_direction": NATOMS // len(ANGLES),
        "target_local_moment_Bohr": 2.0,
        "initial_direction_vector_sum": (directions * (NATOMS // len(ANGLES))).sum(axis=0).tolist(),
        "initial_target_moment_vector_sum_Bohr": (
            2.0 * directions * (NATOMS // len(ANGLES))
        ).sum(axis=0).tolist(),
        "label_assignment": label_summary,
        "reason_for_repetition": "QE 7.5 ntypx=10; eight artificial Fe species preserve the existing patched executable",
    }
    (ROOT / "magnetic_seed_summary.json").write_text(
        json.dumps(magnetic_summary, indent=2) + "\n"
    )
    velocity_lines = ["ATOMIC_VELOCITIES { a.u }"] + [
        f"{label} {velocity[0]: .16e} {velocity[1]: .16e} {velocity[2]: .16e}"
        for label, velocity in zip(labels, velocities)
    ]
    (ROOT / "atomic_velocities_4000K.txt").write_text("\n".join(velocity_lines) + "\n")

    folder = ROOT / "cases" / CASE
    folder.mkdir(exist_ok=True)
    if list(folder.glob("*.out")) or (folder / "qe_tmp").exists():
        raise RuntimeError(f"Refusing to overwrite results in {folder}")
    prefix = "Fe_bcc32_r14_a2p39_4x2x2_m2_l0100_md400"
    for stage_index, (stem, lambda_value, threshold, maxstep, degauss, beta) in enumerate(STAGES):
        lines = control("scf", prefix) + ["/"]
        lines += system(lambda_value, degauss)
        lines += electrons(stage_index > 0, threshold, maxstep, beta)
        lines += cards(positions, labels, None)
        (folder / f"{stem}.in").write_text("\n".join(lines))
    md_lines = control("md", prefix) + ["   dt=20.6700d0", "   nstep=400", "/"]
    md_lines += system(0.100, 0.025334)
    md_lines += electrons(True, 4.0e-3, 500, 0.010)
    md_lines += [
        "&IONS", "   ion_dynamics='verlet'", "   pot_extrapolation='atomic'",
        "   wfc_extrapolation='none'", "   ion_temperature='svr'",
        "   ion_velocities='from_input'", "   tempw=4000.0", "   nraise=20", "/",
    ]
    md_lines += cards(positions, labels, velocities)
    (folder / "md_400steps.in").write_text("\n".join(md_lines))

    metadata = {
        "serial_index": 0,
        "case": CASE,
        "lattice_a_A": LATTICE_A,
        "target_temperature_K": TARGET_TEMPERATURE,
        "target_local_moment_Bohr": 2.0,
        "starting_magnetization": 0.125,
        "lambda_stage1": 0.005,
        "lambda_stage2": 0.020,
        "lambda_stage3_md": 0.100,
        "md_steps": 400,
        "velocity_seed": VELOCITY_SEED,
        "qe_requirement": "private patched QE 7.5 mcons=zv*starting_magnetization",
        "lambda_bootstrap": 0.001,
        "unconstrained_bootstrap": True,
        "unconstrained_bootstrap_degauss_Ry": 0.080000,
        "unconstrained_bootstrap_conv_thr_Ry": 0.020000,
        "unconstrained_bootstrap_mixing_beta": 0.050,
        "supercell": "4x2x2",
        "natoms": NATOMS,
        "nraise": 20,
        "thermal_source_frame": THERMAL_FRAME_INDEX,
        "displacement_selection_seed": DISPLACEMENT_SEED,
        "magnetic_label_seed": LABEL_SEED,
        "bootstrap_degauss_Ry": 0.050000,
        "production_degauss_Ry": 0.025334,
        "mixing_mode": "plain",
        "mixing_beta": 0.010,
        "mixing_ndim": 40,
    }
    (folder / "case_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    with (ROOT / "case_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metadata), lineterminator="\n")
        writer.writeheader()
        writer.writerow(metadata)
    (ROOT / "cases.txt").write_text(CASE + "\n")
    print("Generated one Round-14 unconstrained-bootstrap BCC 4x2x2 32-atom case")
    print(f"Thermal start: RMS displacement={thermal_summary['rms_applied_displacement_A']:.4f} A; "
          f"minimum pair distance={thermal_summary['minimum_initial_pair_distance_A']:.4f} A")
    print(f"Velocities: T={velocity_summary['measured_temperature_K']:.6f} K; seed={VELOCITY_SEED}")
    print("Initial magnetic target-vector sum:", magnetic_summary["initial_target_moment_vector_sum_Bohr"])
    print("Magnetic labels: two copies of every direction on each BCC sublattice")


if __name__ == "__main__":
    main()
