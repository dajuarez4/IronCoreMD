#!/usr/bin/env python3
"""Audit the returned Round-19 QE output and render a partial-results dashboard."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-round19-dashboard")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
CASE = "a2p39_6x2x2_m2_48species_lambda0p100_md1000"
DEFAULT_INPUT = ROOT / "cases" / CASE / "md_1000steps.in"
DEFAULT_OUTPUT = ROOT / "md_1000steps.out"
DEFAULT_PNG = ROOT / "round19_partial_dashboard.png"
DEFAULT_JSON = ROOT / "round19_analysis_summary.json"
DEFAULT_CSV = ROOT / "round19_md_records.csv"
NATOMS = 48
EXPECTED_STEPS = 1000
DT_AU = 20.67
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
AU_TIME_SECONDS = 2.4188843265864e-17
COLORS = {
    "navy": "#0b2d4d", "orange": "#f47c20", "purple": "#7b4ab5",
    "green": "#2a9d68", "blue": "#2878b5", "red": "#c44e52",
    "gray": "#687386", "light": "#d8dee7", "gold": "#d9a500",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--png", type=Path, default=DEFAULT_PNG)
    parser.add_argument("--summary", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def scalars(pattern: str, text: str, flags: int = re.I) -> np.ndarray:
    return np.asarray([number(value) for value in re.findall(pattern, text, flags)], dtype=float)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_positions_and_cell(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text().splitlines()
    position_start = lines.index("ATOMIC_POSITIONS crystal") + 1
    positions = np.asarray(
        [[number(value) for value in row.split()[1:4]]
         for row in lines[position_start:position_start + NATOMS]], dtype=float
    )
    cell_start = lines.index("CELL_PARAMETERS angstrom") + 1
    cell = np.asarray(
        [[number(value) for value in row.split()[:3]]
         for row in lines[cell_start:cell_start + 3]], dtype=float
    )
    return positions, cell


def position_blocks(text: str) -> np.ndarray:
    lines = text.splitlines()
    blocks: list[list[list[float]]] = []
    for index, line in enumerate(lines):
        if not line.startswith("ATOMIC_POSITIONS (crystal)"):
            continue
        rows: list[list[float]] = []
        for candidate in lines[index + 1:index + 1 + NATOMS]:
            fields = candidate.split()
            if len(fields) < 4:
                break
            rows.append([number(value) for value in fields[1:4]])
        if len(rows) == NATOMS:
            blocks.append(rows)
    return np.asarray(blocks, dtype=float).reshape(-1, NATOMS, 3)


def vector_frames(text: str, label: str) -> np.ndarray:
    rows = re.findall(
        rf"^\s*{re.escape(label)}\s*:\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        text, re.I | re.M,
    )
    values = np.asarray([[number(value) for value in row] for row in rows], dtype=float)
    complete = len(values) // NATOMS
    return values[:complete * NATOMS].reshape(complete, NATOMS, 3)


def converged_scf_history(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    magnetization: list[list[float]] = []
    accuracy: list[float] = []
    cycles: list[int] = []
    current_magnetization = [np.nan, np.nan, np.nan]
    current_accuracy = np.nan
    cycle_iterations = 0
    for line in text.splitlines():
        match = re.search(
            rf"total magnetization\s*=\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})", line
        )
        if match:
            current_magnetization = [number(value) for value in match.groups()]
        match = re.search(rf"estimated scf accuracy\s*<\s*({FLOAT})\s*Ry", line)
        if match:
            current_accuracy = number(match.group(1))
        if "iteration #" in line:
            cycle_iterations += 1
        if "convergence has been achieved" in line:
            magnetization.append(current_magnetization)
            accuracy.append(current_accuracy)
            cycles.append(cycle_iterations)
            cycle_iterations = 0
    return (
        np.asarray(magnetization, dtype=float), np.asarray(accuracy, dtype=float),
        np.asarray(cycles, dtype=int), cycle_iterations,
    )


def periodic_displacements(positions: np.ndarray, reference: np.ndarray, cell: np.ndarray) -> np.ndarray:
    delta = positions - reference[None, :, :]
    delta -= np.rint(delta)
    return delta @ cell


def minimum_separations(positions: np.ndarray, cell: np.ndarray) -> np.ndarray:
    result = []
    upper = np.triu_indices(NATOMS, 1)
    for frame in positions:
        delta = frame[:, None, :] - frame[None, :, :]
        delta -= np.rint(delta)
        distances = np.linalg.norm(delta @ cell, axis=2)
        result.append(float(np.min(distances[upper])))
    return np.asarray(result)


def align(values: np.ndarray, count: int, label: str) -> np.ndarray:
    if len(values) < count:
        raise ValueError(f"{label} has {len(values)} records; expected at least {count}")
    return values[:count]


def style_axis(axis) -> None:
    axis.grid(color=COLORS["light"], linewidth=0.8, alpha=0.85)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=9)


def finite(value: float) -> float | None:
    return float(value) if np.isfinite(value) else None


def main() -> int:
    args = arguments()
    text = args.output.read_text(errors="replace")
    initial_positions, cell = input_positions_and_cell(args.input)
    positions = position_blocks(text)
    total_magnetization, accuracy, cycles, unfinished_iterations = converged_scf_history(text)
    steps = len(total_magnetization)
    if steps == 0:
        raise ValueError("No converged Round-19 MD SCFs were found")

    positions = align(positions, steps, "positions")
    temperature = align(
        scalars(rf"^\s*temperature\s*=\s*({FLOAT})\s*K", text, re.I | re.M), steps, "temperature"
    )
    pressure_kbar = align(
        scalars(rf"(?m)(?:^|\s)P\s*=\s*({FLOAT})", text), steps, "pressure"
    )
    total_force = align(
        scalars(rf"Total force\s*=\s*({FLOAT})", text), steps, "total force"
    )
    force_correction = align(
        scalars(rf"Total SCF correction\s*=\s*({FLOAT})", text), steps, "SCF force correction"
    )
    local_frames = vector_frames(text, "magnetization")
    if len(local_frames) == steps + 1:
        local_frames = local_frames[1:]
    local_frames = align(local_frames, steps, "local magnetization")
    constrained_frames = vector_frames(text, "constrained moment")
    if len(constrained_frames) == steps + 1:
        constrained_frames = constrained_frames[1:]
    constrained_frames = align(constrained_frames, steps, "printed constrained moment")

    dt_ps = DT_AU * 2.0 * AU_TIME_SECONDS * 1.0e12
    time_ps = (np.arange(steps) + 1) * dt_ps
    local_magnitudes = np.linalg.norm(local_frames, axis=2)
    mean_local_moment = np.mean(local_magnitudes, axis=1)
    printed_constraint_magnitude = np.mean(np.linalg.norm(constrained_frames, axis=2), axis=1)
    total_magnitude = np.linalg.norm(total_magnetization, axis=1)
    displacements = periodic_displacements(positions, initial_positions, cell)
    rms_displacement = np.sqrt(np.mean(np.sum(displacements**2, axis=2), axis=1))
    max_displacement = np.max(np.linalg.norm(displacements, axis=2), axis=1)
    minimum_distance = minimum_separations(positions, cell)
    correction_ratio = force_correction / total_force

    job_done = "JOB DONE" in text
    fatal_errors = len(re.findall(r"Error in routine|convergence NOT achieved|Cholesky|cdiaghg", text, re.I))
    c_bands_warnings = len(re.findall(r"c_bands:.*not converged", text, re.I))
    large_force_warnings = len(re.findall(r"SCF correction compared to forces is large", text, re.I))
    complete = job_done and steps == EXPECTED_STEPS and fatal_errors == 0
    phonon_ready = complete and steps >= 1000 and float(np.max(correction_ratio)) <= 0.1
    constraint_scale_ok = bool(np.allclose(printed_constraint_magnitude, 2.0, rtol=0.01, atol=0.01))

    summary = {
        "analysis_scope": "partial returned MD output",
        "source": str(args.output.name),
        "source_sha256": sha256(args.output),
        "source_size_bytes": args.output.stat().st_size,
        "atoms": NATOMS,
        "expected_md_steps": EXPECTED_STEPS,
        "completed_md_steps": steps,
        "completion_fraction": steps / EXPECTED_STEPS,
        "simulated_time_ps": float(time_ps[-1]),
        "job_done": job_done,
        "complete": complete,
        "unfinished_scf_iterations_after_last_step": unfinished_iterations,
        "fatal_error_markers": fatal_errors,
        "c_bands_nonconvergence_warnings": c_bands_warnings,
        "large_scf_force_correction_warnings": large_force_warnings,
        "temperature_K": {
            "mean": float(np.mean(temperature)), "minimum": float(np.min(temperature)),
            "maximum": float(np.max(temperature)), "final": float(temperature[-1]),
        },
        "pressure_GPa": {
            "mean": float(np.mean(pressure_kbar) / 10.0),
            "minimum": float(np.min(pressure_kbar) / 10.0),
            "maximum": float(np.max(pressure_kbar) / 10.0),
            "final": float(pressure_kbar[-1] / 10.0),
        },
        "scf": {
            "mean_iterations": float(np.mean(cycles)), "maximum_iterations": int(np.max(cycles)),
            "mean_final_accuracy_Ry": float(np.mean(accuracy)),
            "maximum_final_accuracy_Ry": float(np.max(accuracy)),
        },
        "magnetization": {
            "mean_local_moment_Bohr_per_Fe": float(np.mean(mean_local_moment)),
            "final_local_moment_Bohr_per_Fe": float(mean_local_moment[-1]),
            "requested_constraint_magnitude_Bohr_per_Fe": 2.0,
            "mean_printed_constraint_magnitude": float(np.mean(printed_constraint_magnitude)),
            "constraint_scale_correct": constraint_scale_ok,
            "constraint_assessment": (
                "FAIL: QE printed a constrained-vector magnitude of 0.125 rather than 2.0. "
                "The expected patched mcons=Zval*starting_magnetization behavior was not active "
                "in the executable used for this run."
            ),
            "mean_total_vector_Bohr_per_cell": [float(value) for value in np.mean(total_magnetization, axis=0)],
            "mean_total_magnitude_Bohr_per_cell": float(np.mean(total_magnitude)),
            "maximum_total_magnitude_Bohr_per_cell": float(np.max(total_magnitude)),
        },
        "structure": {
            "maximum_rms_displacement_from_start_A": float(np.max(rms_displacement)),
            "final_rms_displacement_from_start_A": float(rms_displacement[-1]),
            "minimum_periodic_Fe_Fe_distance_A": float(np.min(minimum_distance)),
        },
        "force_quality": {
            "mean_total_force_Ry_per_Bohr": float(np.mean(total_force)),
            "mean_scf_correction_ratio": float(np.mean(correction_ratio)),
            "maximum_scf_correction_ratio": float(np.max(correction_ratio)),
        },
        "phonons": {
            "computed": False,
            "ready_for_tdep_fit": phonon_ready,
            "reason": (
                f"The returned trajectory has only {steps}/{EXPECTED_STEPS} completed MD steps "
                "and ends mid-SCF; "
                "it is too short for a converged TDEP force-constant fit or velocity-autocorrelation DOS. "
                "The run also used the wrong 0.125 constrained-moment scale, and its MD forces use "
                "conv_thr=4e-3 Ry with non-negligible SCF corrections."
            ),
        },
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")

    with args.csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "step", "time_ps", "temperature_K", "pressure_GPa", "scf_iterations",
            "final_accuracy_Ry", "Mx_Bohr_cell", "My_Bohr_cell", "Mz_Bohr_cell",
            "M_magnitude_Bohr_cell", "mean_local_moment_Bohr_Fe", "total_force_Ry_Bohr",
            "scf_force_correction_Ry_Bohr", "scf_force_correction_ratio",
            "rms_displacement_A", "max_displacement_A", "minimum_Fe_Fe_distance_A",
        ])
        for index in range(steps):
            writer.writerow([
                index + 1, time_ps[index], temperature[index], pressure_kbar[index] / 10.0,
                cycles[index], accuracy[index], *total_magnetization[index], total_magnitude[index],
                mean_local_moment[index], total_force[index], force_correction[index],
                correction_ratio[index], rms_displacement[index], max_displacement[index],
                minimum_distance[index],
            ])

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "axes.titleweight": "bold",
        "axes.titlecolor": COLORS["navy"], "axes.labelcolor": COLORS["navy"],
    })
    figure, axes = plt.subplots(3, 3, figsize=(18, 12), facecolor="white")
    figure.subplots_adjust(left=0.06, right=0.97, bottom=0.07, top=0.82, wspace=0.28, hspace=0.48)
    figure.suptitle("BCC Fe Round 19 — 48-atom noncollinear MD audit", fontsize=23,
                    fontweight="bold", color=COLORS["navy"], y=0.965)
    status_color = COLORS["green"] if complete else COLORS["red"]
    status_text = (
        "COMPLETE" if complete and constraint_scale_ok else
        "INVALID / INCOMPLETE — wrong constraint scale; output ends during next SCF"
    )
    figure.text(0.5, 0.925, status_text, ha="center", va="center", fontsize=15,
                fontweight="bold", color=status_color)
    figure.text(
        0.5, 0.892,
        f"{steps}/{EXPECTED_STEPS} MD steps ({100*steps/EXPECTED_STEPS:.1f}%)  |  "
        f"{time_ps[-1]:.4f} ps  |  38 MPI tasks  |  a=2.39 Å, 6×2×2  |  λ=0.10 Ry",
        ha="center", va="center", fontsize=11, color=COLORS["gray"],
    )

    axis = axes[0, 0]
    axis.plot(time_ps, temperature, color=COLORS["orange"], marker="o", ms=3, lw=1.8)
    axis.axhline(4000, color=COLORS["gray"], ls="--", lw=1.2, label="target")
    axis.set(title="Ionic temperature", xlabel="Time (ps)", ylabel="Temperature (K)")
    axis.legend(frameon=False, fontsize=8); style_axis(axis)

    axis = axes[0, 1]
    axis.plot(time_ps, pressure_kbar / 10.0, color=COLORS["purple"], marker="o", ms=3, lw=1.8)
    axis.set(title="Pressure", xlabel="Time (ps)", ylabel="Pressure (GPa)"); style_axis(axis)

    axis = axes[0, 2]
    for component, label, color in zip(total_magnetization.T, ("Mx", "My", "Mz"),
                                        (COLORS["red"], COLORS["green"], COLORS["blue"])):
        axis.plot(time_ps, component, label=label, color=color, lw=1.5)
    axis.plot(time_ps, total_magnitude, label="|M|", color=COLORS["navy"], lw=2.2)
    axis.axhline(0, color=COLORS["gray"], lw=0.8)
    axis.set(title="Total magnetization", xlabel="Time (ps)", ylabel="Bohr magneton/cell")
    axis.legend(frameon=False, ncol=4, fontsize=8); style_axis(axis)

    axis = axes[1, 0]
    axis.plot(time_ps, mean_local_moment, color=COLORS["purple"], marker="o", ms=3, lw=1.8)
    axis.plot(time_ps, printed_constraint_magnitude, color=COLORS["red"], ls="--", lw=1.4,
              label="printed constraint (wrong)")
    axis.axhline(2.0, color=COLORS["gray"], ls=":", lw=1.4, label="requested target")
    axis.set(title="Mean local moment", xlabel="Time (ps)", ylabel="Bohr magneton/Fe")
    axis.legend(frameon=False, fontsize=8); style_axis(axis)

    axis = axes[1, 1]
    axis.plot(np.arange(1, steps + 1), cycles, color=COLORS["blue"], marker="o", ms=3, lw=1.5)
    axis.set(title="Electronic convergence", xlabel="Completed MD SCF", ylabel="SCF iterations")
    twin = axis.twinx()
    twin.semilogy(np.arange(1, steps + 1), accuracy, color=COLORS["orange"], marker=".", lw=1.2)
    twin.axhline(4e-3, color=COLORS["gray"], ls="--", lw=1)
    twin.set_ylabel("Final accuracy (Ry)", color=COLORS["orange"], fontsize=9)
    style_axis(axis)

    axis = axes[1, 2]
    axis.plot(time_ps, total_force, label="total force", color=COLORS["blue"], lw=1.8)
    axis.plot(time_ps, force_correction, label="SCF correction", color=COLORS["red"], lw=1.8)
    axis.set(title="Force quality", xlabel="Time (ps)", ylabel="Ry/Bohr")
    axis.legend(frameon=False, fontsize=8); style_axis(axis)

    axis = axes[2, 0]
    axis.plot(time_ps, rms_displacement, label="RMS from start", color=COLORS["green"], lw=1.8)
    axis.plot(time_ps, max_displacement, label="maximum", color=COLORS["orange"], lw=1.5)
    axis.set(title="Structural displacement", xlabel="Time (ps)", ylabel="Displacement (Å)")
    axis.legend(frameon=False, fontsize=8); style_axis(axis)

    axis = axes[2, 1]
    axis.plot(time_ps, minimum_distance, color=COLORS["navy"], marker="o", ms=3, lw=1.8)
    axis.set(title="Closest periodic Fe–Fe pair", xlabel="Time (ps)", ylabel="Distance (Å)")
    style_axis(axis)

    axis = axes[2, 2]
    axis.set_facecolor("#fff5f5")
    axis.set_xticks([]); axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_color("#e8a0a0"); spine.set_linewidth(1.5)
    axis.text(0.5, 0.84, "PHONONS / TDEP", transform=axis.transAxes, ha="center", va="center",
              fontsize=16, fontweight="bold", color=COLORS["red"])
    axis.text(0.5, 0.67, "NOT COMPUTED", transform=axis.transAxes, ha="center", va="center",
              fontsize=19, fontweight="bold", color=COLORS["red"])
    axis.text(
        0.5, 0.37,
        f"Only {steps}/{EXPECTED_STEPS} configurations\n"
        f"Trajectory length: {time_ps[-1]:.4f} ps\n"
        f"Mean SCF-force/force ratio: {np.mean(correction_ratio):.3f}\n\n"
        "No defensible dispersion or DOS can be\n"
        "fit from this interrupted trajectory.",
        transform=axis.transAxes, ha="center", va="center", fontsize=11,
        color=COLORS["navy"], linespacing=1.45,
    )

    figure.text(
        0.06, 0.018,
        f"Mean T={np.mean(temperature):.0f} K | mean P={np.mean(pressure_kbar)/10:.1f} GPa | "
        f"mean local moment={np.mean(mean_local_moment):.3f} μB/Fe | printed constraint="
        f"{np.mean(printed_constraint_magnitude):.3f} instead of 2.000 μB | "
        f"min Fe–Fe={np.min(minimum_distance):.3f} Å | output SHA-256 {summary['source_sha256'][:12]}…",
        ha="left", va="bottom", fontsize=9.5, color=COLORS["gray"],
    )
    args.png.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.png, dpi=160, facecolor="white")
    plt.close(figure)

    print(args.png.resolve())
    print(args.summary.resolve())
    print(args.csv.resolve())
    print(json.dumps({
        "steps": steps, "complete": complete, "unfinished_scf_iterations": unfinished_iterations,
        "mean_temperature_K": finite(np.mean(temperature)),
        "mean_pressure_GPa": finite(np.mean(pressure_kbar) / 10.0),
        "mean_local_moment_Bohr_per_Fe": finite(np.mean(mean_local_moment)),
        "phonons_computed": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
