#!/usr/bin/env python3
"""Collect Round-20 QE outputs into one small, safe-loadable NPZ archive.

Run this in the Round-20 directory on Jakar after copying this script there.
Incomplete/failed trajectories are retained; use --require-complete to reject them.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
SCHEMA = "ironcoremd.qe.bcc16_dlm8_round20.v1"
NATOMS = 16
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
SOURCES = (("stage1_4e-3", "scf"), ("stage2_1e-3", "scf"),
           ("stage3_3e-4", "scf"), ("md_100steps", "md"))


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def values(pattern: str, text: str) -> list[float]:
    return [number(item) for item in re.findall(pattern, text, re.I)]


def uarray(items: list[str]) -> np.ndarray:
    width = max((len(item) for item in items), default=1)
    return np.asarray(items, dtype=f"U{width}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def flatten(series: list[list[float]]) -> tuple[np.ndarray, np.ndarray]:
    offsets = [0]
    merged: list[float] = []
    for item in series:
        merged.extend(item)
        offsets.append(len(merged))
    return np.asarray(merged, dtype=np.float64), np.asarray(offsets, dtype=np.int64)


def vector_frames(text: str, label: str) -> list[np.ndarray]:
    rows = re.findall(
        rf"^\s*{re.escape(label)}\s*:\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        text, re.I | re.M,
    )
    numeric = [[number(value) for value in row] for row in rows]
    return [np.asarray(numeric[i:i + NATOMS], dtype=np.float32)
            for i in range(0, len(numeric) - NATOMS + 1, NATOMS)]


def vector_means(text: str, label: str) -> list[float]:
    return [float(np.linalg.norm(frame, axis=1).mean())
            for frame in vector_frames(text, label)]


def position_blocks(text: str) -> list[tuple[str, np.ndarray]]:
    lines = text.splitlines()
    blocks: list[tuple[str, np.ndarray]] = []
    for index, line in enumerate(lines):
        if not line.strip().upper().startswith("ATOMIC_POSITIONS"):
            continue
        header = line.lower()
        unit = next((u for u in ("crystal", "angstrom", "bohr") if u in header), "alat")
        rows = []
        for candidate in lines[index + 1:index + 1 + NATOMS]:
            fields = candidate.split()
            try:
                rows.append([number(fields[1]), number(fields[2]), number(fields[3])])
            except (IndexError, ValueError):
                break
        if len(rows) == NATOMS:
            blocks.append((unit, np.asarray(rows, dtype=np.float32)))
    return blocks


def force_blocks(text: str) -> list[np.ndarray]:
    matches = re.findall(
        rf"atom\s+(\d+)\s+type\s+\d+\s+force\s*=\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        text, re.I,
    )
    frames: list[np.ndarray] = []
    current: dict[int, list[float]] = {}
    previous = 0
    for row in matches:
        atom = int(row[0])
        if current and atom <= previous:
            if len(current) == NATOMS:
                frames.append(np.asarray([current[i] for i in range(1, NATOMS + 1)], dtype=np.float32))
            current = {}
        previous = atom
        current[atom] = [number(value) for value in row[1:]]
    if len(current) == NATOMS:
        frames.append(np.asarray([current[i] for i in range(1, NATOMS + 1)], dtype=np.float32))
    return frames


def cycles(text: str) -> list[float]:
    result: list[float] = []
    count = 0
    for line in text.splitlines():
        if "iteration #" in line:
            count += 1
        if "convergence has been achieved" in line or "convergence NOT achieved" in line:
            if count:
                result.append(float(count))
                count = 0
    return result


def source_status(text: str, exists: bool, kind: str, expected_steps: int) -> str:
    if not exists:
        return "NOT_STARTED"
    fatal = bool(re.search(r"Error in routine|convergence NOT achieved|Cholesky|cdiaghg", text, re.I))
    if fatal:
        return "FAILED"
    if kind == "scf":
        return "CONVERGED" if "convergence has been achieved" in text else "INCOMPLETE"
    steps = len(re.findall(r"Entering Dynamics", text))
    if "JOB DONE" in text:
        return "DONE" if steps == expected_steps else "PARTIAL"
    return "RUNNING_OR_PARTIAL"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-terminal", action="store_true",
                        help="fail if any MD output is missing or still running")
    parser.add_argument("--require-complete", action="store_true",
                        help="fail unless every trajectory completed all 100 steps")
    parser.add_argument("-o", "--output", type=Path, default=ROOT / "round20_results_light.npz")
    args = parser.parse_args()

    with (ROOT / "case_manifest.csv").open(newline="") as handle:
        cases = list(csv.DictReader(handle))

    source_case: list[int] = []
    source_stem: list[str] = []
    source_kind: list[str] = []
    status: list[str] = []
    iteration_count: list[int] = []
    last_accuracy: list[float] = []
    md_steps: list[int] = []
    force_warnings: list[int] = []
    c_bands_warnings: list[int] = []
    errors: list[int] = []
    input_paths: list[str] = []
    input_texts: list[str] = []
    output_tails: list[str] = []
    series: dict[str, list[list[float]]] = {name: [] for name in (
        "accuracy_Ry", "scf_cycle_iterations", "temperature_K", "mean_local_moment_Bohr",
        "printed_constrained_moment_Bohr", "energy_Ry", "total_force_Ry_Bohr", "pressure_kbar")}
    position_source: list[int] = []
    position_units: list[str] = []
    positions: list[np.ndarray] = []
    atomic_force_source: list[int] = []
    atomic_forces: list[np.ndarray] = []
    local_vector_source: list[int] = []
    local_vectors: list[np.ndarray] = []
    constraint_vector_source: list[int] = []
    constraint_vectors: list[np.ndarray] = []

    for case_index, case in enumerate(cases):
        folder = ROOT / "cases" / case["case"]
        expected = int(case["md_steps"])
        for stem, kind in SOURCES:
            source_index = len(source_case)
            inp, out = folder / f"{stem}.in", folder / f"{stem}.out"
            text = out.read_text(errors="replace") if out.is_file() else ""
            accuracies = values(rf"estimated scf accuracy\s*<\s*({FLOAT})\s*Ry", text)
            source_case.append(case_index)
            source_stem.append(stem)
            source_kind.append(kind)
            status.append(source_status(text, out.is_file(), kind, expected))
            iteration_count.append(len(re.findall(r"iteration #", text)))
            last_accuracy.append(accuracies[-1] if accuracies else math.nan)
            md_steps.append(len(re.findall(r"Entering Dynamics", text)) if kind == "md" else 0)
            force_warnings.append(len(re.findall(r"SCF correction compared to forces is large", text, re.I)))
            c_bands_warnings.append(len(re.findall(r"c_bands:.*not converged", text, re.I)))
            errors.append(len(re.findall(r"Error in routine|Cholesky|cdiaghg", text, re.I)))
            input_paths.append(str(inp.relative_to(ROOT)))
            input_texts.append(inp.read_text(errors="replace") if inp.is_file() else "")
            output_tails.append("\n".join(text.splitlines()[-80:]))
            series["accuracy_Ry"].append(accuracies)
            series["scf_cycle_iterations"].append(cycles(text))
            series["temperature_K"].append(values(rf"temperature\s*=\s*({FLOAT})\s*K", text))
            series["mean_local_moment_Bohr"].append(vector_means(text, "magnetization"))
            series["printed_constrained_moment_Bohr"].append(vector_means(text, "constrained moment"))
            series["energy_Ry"].append(values(rf"!\s*total energy\s*=\s*({FLOAT})\s*Ry", text))
            series["total_force_Ry_Bohr"].append(values(rf"Total force\s*=\s*({FLOAT})", text))
            series["pressure_kbar"].append(values(rf"(?m)(?:^|\s)P\s*=\s*({FLOAT})", text))
            if kind == "md":
                for unit, block in position_blocks(text):
                    position_source.append(source_index); position_units.append(unit); positions.append(block)
                for block in force_blocks(text):
                    atomic_force_source.append(source_index); atomic_forces.append(block)
                for block in vector_frames(text, "magnetization"):
                    local_vector_source.append(source_index); local_vectors.append(block)
                for block in vector_frames(text, "constrained moment"):
                    constraint_vector_source.append(source_index); constraint_vectors.append(block)

    md_states = [status[i * len(SOURCES) + 3] for i in range(len(cases))]
    terminal_states = {"DONE", "PARTIAL", "FAILED"}
    if args.require_terminal and not all(state in terminal_states for state in md_states):
        raise SystemExit(f"Refusing nonterminal archive: MD states={md_states}")
    if args.require_complete and not all(state == "DONE" for state in md_states):
        raise SystemExit(f"Refusing incomplete archive: MD states={md_states}")

    payload: dict[str, np.ndarray] = {
        "schema_version": uarray([SCHEMA]),
        "created_utc": uarray([datetime.now(timezone.utc).isoformat()]),
        "case_name": uarray([case["case"] for case in cases]),
        "spin_seed": np.asarray([int(case["spin_seed"]) for case in cases], dtype=np.int32),
        "velocity_seed": np.asarray([int(case["velocity_seed"]) for case in cases], dtype=np.int32),
        "target_temperature_K": np.asarray([float(case["target_temperature_K"]) for case in cases]),
        "target_local_moment_Bohr": np.asarray([float(case["target_local_moment_Bohr"]) for case in cases]),
        "source_case_index": np.asarray(source_case, dtype=np.int16),
        "source_stem": uarray(source_stem), "source_kind": uarray(source_kind), "status": uarray(status),
        "iteration_count": np.asarray(iteration_count, dtype=np.int32),
        "last_accuracy_Ry": np.asarray(last_accuracy), "md_steps": np.asarray(md_steps, dtype=np.int16),
        "force_correction_warning_count": np.asarray(force_warnings, dtype=np.int16),
        "c_bands_warning_count": np.asarray(c_bands_warnings, dtype=np.int16),
        "qe_error_count": np.asarray(errors, dtype=np.int16),
        "input_relative_path": uarray(input_paths), "input_text": uarray(input_texts),
        "output_tail": uarray(output_tails),
        "position_source_index": np.asarray(position_source, dtype=np.int16),
        "position_unit": uarray(position_units),
        "positions": np.asarray(positions, dtype=np.float32).reshape(-1, NATOMS, 3),
        "atomic_force_source_index": np.asarray(atomic_force_source, dtype=np.int16),
        "atomic_forces_Ry_Bohr": np.asarray(atomic_forces, dtype=np.float32).reshape(-1, NATOMS, 3),
        "local_vector_source_index": np.asarray(local_vector_source, dtype=np.int16),
        "local_magnetization_vectors_Bohr": np.asarray(local_vectors, dtype=np.float32).reshape(-1, NATOMS, 3),
        "constraint_vector_source_index": np.asarray(constraint_vector_source, dtype=np.int16),
        "printed_constraint_vectors_Bohr": np.asarray(constraint_vectors, dtype=np.float32).reshape(-1, NATOMS, 3),
    }
    for name, item in series.items():
        payload[name], payload[f"{name}_offset"] = flatten(item)

    support_names = ("README.md", "case_manifest.csv", "cases.txt", "ensemble_summary.json",
                     "structure_summary.json", "run_round20_pilot_jakar.sbatch",
                     "run_round20_ensemble_jakar.sbatch", "check_round20_jakar.sh",
                     "generate_round20_jakar.py")
    support = [ROOT / name for name in support_names if (ROOT / name).is_file()]
    payload["support_relative_path"] = uarray([str(path.relative_to(ROOT)) for path in support])
    payload["support_text"] = uarray([path.read_text(errors="replace") for path in support])
    payload["support_sha256"] = uarray([sha256(path) for path in support])
    payload["metadata_json"] = uarray([json.dumps({
        "safe_load": "numpy.load(path, allow_pickle=False)",
        "full_outputs_included": False, "qe_tmp_included": False,
        "failed_and_partial_frames_retained": True, "natoms": NATOMS,
        "replicas": len(cases), "supercell": "2x2x2", "nraise": 5,
        "units": {"energy": "Ry", "force": "Ry/Bohr", "pressure": "kbar", "temperature": "K"},
    }, sort_keys=True)])

    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **payload)
    with np.load(output, allow_pickle=False) as checked:
        assert not any(checked[name].dtype.hasobject for name in checked.files)
        arrays = len(checked.files)
    print(output)
    print(f"MD states={md_states}")
    print(f"arrays={arrays} size_bytes={output.stat().st_size} sha256={sha256(output)}")


if __name__ == "__main__":
    main()
