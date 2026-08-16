#!/usr/bin/env python3
"""Collect Round-12 SCF/MD400 results into a lightweight safe NPZ archive."""

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
SCHEMA = "ironcoremd.qe.bcc32_patched_qe75_lambda0p100_md400_round12.v1"
NATOMS = 32
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
SOURCES = (("stage1_4e-3", "scf"), ("stage2_1e-3", "scf"),
           ("stage3_3e-4", "scf"), ("md_400steps", "md"))


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
        merged.extend(item); offsets.append(len(merged))
    return np.asarray(merged, dtype=np.float64), np.asarray(offsets, dtype=np.int64)


def vector_means(text: str, label: str) -> list[float]:
    vectors = re.findall(
        rf"^\s*{re.escape(label)}\s*:\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        text, re.I | re.M,
    )
    magnitudes = [math.sqrt(sum(number(x) ** 2 for x in row)) for row in vectors]
    return [sum(magnitudes[i:i + NATOMS]) / NATOMS
            for i in range(0, len(magnitudes) - NATOMS + 1, NATOMS)]


def vector_frames(text: str, label: str) -> list[np.ndarray]:
    vectors = re.findall(
        rf"^\s*{re.escape(label)}\s*:\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})",
        text, re.I | re.M,
    )
    numeric = [[number(value) for value in row] for row in vectors]
    return [np.asarray(numeric[i:i + NATOMS], dtype=np.float32)
            for i in range(0, len(numeric) - NATOMS + 1, NATOMS)]


def position_blocks(text: str) -> list[tuple[str, np.ndarray]]:
    lines = text.splitlines(); blocks: list[tuple[str, np.ndarray]] = []
    for index, line in enumerate(lines):
        if not line.strip().upper().startswith("ATOMIC_POSITIONS"):
            continue
        header = line.lower()
        unit = "crystal" if "crystal" in header else (
            "angstrom" if "angstrom" in header else "bohr" if "bohr" in header else "alat"
        )
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
    frames: list[np.ndarray] = []; current: dict[int, list[float]] = {}; previous = 0
    for row in matches:
        atom = int(row[0])
        if current and atom <= previous:
            if len(current) == NATOMS:
                frames.append(np.asarray([current[i] for i in range(1, NATOMS + 1)], dtype=np.float32))
            current = {}
        previous = atom; current[atom] = [number(value) for value in row[1:]]
    if len(current) == NATOMS:
        frames.append(np.asarray([current[i] for i in range(1, NATOMS + 1)], dtype=np.float32))
    return frames


def cycles(text: str) -> list[float]:
    result: list[float] = []; count = 0
    for line in text.splitlines():
        if "iteration #" in line:
            count += 1
        if "convergence has been achieved" in line or "convergence NOT achieved" in line:
            if count:
                result.append(float(count)); count = 0
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
    return "RUNNING"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-terminal", action="store_true")
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument("-o", "--output", type=Path, default=ROOT / "round12_results_light.npz")
    args = parser.parse_args()
    with (ROOT / "case_manifest.csv").open(newline="") as handle:
        cases = list(csv.DictReader(handle))

    source_case: list[int] = []; source_stem: list[str] = []; source_kind: list[str] = []
    status: list[str] = []; iteration_count: list[int] = []; last_accuracy: list[float] = []
    md_steps: list[int] = []; warnings: list[int] = []; cholesky: list[int] = []; errors: list[int] = []
    input_paths: list[str] = []; input_texts: list[str] = []; output_tails: list[str] = []
    accuracy_series: list[list[float]] = []; cycle_series: list[list[float]] = []
    temperature_series: list[list[float]] = []; moment_series: list[list[float]] = []; constraint_series: list[list[float]] = []
    energy_series: list[list[float]] = []; force_series: list[list[float]] = []; pressure_series: list[list[float]] = []
    position_source: list[int] = []; position_units: list[str] = []; positions: list[np.ndarray] = []
    atomic_force_source: list[int] = []; atomic_forces: list[np.ndarray] = []
    local_vector_source: list[int] = []; local_vectors: list[np.ndarray] = []
    constraint_vector_source: list[int] = []; constraint_vectors: list[np.ndarray] = []

    for case_index, case in enumerate(cases):
        folder = ROOT / "cases" / case["case"]
        expected = int(case["md_steps"])
        for stem, kind in SOURCES:
            source_index = len(source_case)
            inp = folder / f"{stem}.in"; out = folder / f"{stem}.out"
            text = out.read_text(errors="replace") if out.is_file() else ""
            accuracies = values(rf"estimated scf accuracy\s*<\s*({FLOAT})\s*Ry", text)
            temperatures = values(rf"temperature\s*=\s*({FLOAT})\s*K", text)
            energies = values(rf"!\s*total energy\s*=\s*({FLOAT})\s*Ry", text)
            forces = values(rf"Total force\s*=\s*({FLOAT})", text)
            pressures = values(rf"(?m)(?:^|\s)P\s*=\s*({FLOAT})", text)
            source_case.append(case_index); source_stem.append(stem); source_kind.append(kind)
            status.append(source_status(text, out.is_file(), kind, expected))
            iteration_count.append(len(re.findall(r"iteration #", text)))
            last_accuracy.append(accuracies[-1] if accuracies else math.nan)
            md_steps.append(len(re.findall(r"Entering Dynamics", text)) if kind == "md" else 0)
            warnings.append(len(re.findall(r"c_bands:.*not converged", text, re.I)))
            cholesky.append(len(re.findall(r"Cholesky|cdiaghg", text, re.I)))
            errors.append(len(re.findall(r"Error in routine", text, re.I)))
            input_paths.append(str(inp.relative_to(ROOT))); input_texts.append(inp.read_text())
            output_tails.append("\n".join(text.splitlines()[-80:]))
            accuracy_series.append(accuracies); cycle_series.append(cycles(text))
            temperature_series.append(temperatures); moment_series.append(vector_means(text, "magnetization"))
            constraint_series.append(vector_means(text, "constrained moment"))
            energy_series.append(energies); force_series.append(forces); pressure_series.append(pressures)
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
    terminal = all(state in ("DONE", "PARTIAL", "FAILED") for state in md_states)
    complete = all(state == "DONE" for state in md_states)
    if args.require_terminal and not terminal:
        raise SystemExit(f"Refusing nonterminal archive: MD states={md_states}")
    if args.require_complete and not complete:
        raise SystemExit(f"Refusing incomplete archive: MD states={md_states}")

    payload: dict[str, np.ndarray] = {
        "schema_version": uarray([SCHEMA]),
        "created_utc": uarray([datetime.now(timezone.utc).isoformat()]),
        "case_count": np.asarray([len(cases)], dtype=np.int32),
        "case_name": uarray([case["case"] for case in cases]),
        "final_lambda_Ry": np.asarray([float(case["lambda_stage3_md"]) for case in cases]),
        "target_local_moment_Bohr": np.asarray([float(case["target_local_moment_Bohr"]) for case in cases]),
        "source_case_index": np.asarray(source_case, dtype=np.int16),
        "source_stem": uarray(source_stem), "source_kind": uarray(source_kind), "status": uarray(status),
        "iteration_count": np.asarray(iteration_count, dtype=np.int32),
        "last_accuracy_Ry": np.asarray(last_accuracy), "md_steps": np.asarray(md_steps, dtype=np.int16),
        "c_bands_warning_count": np.asarray(warnings, dtype=np.int16),
        "cholesky_count": np.asarray(cholesky, dtype=np.int16), "qe_error_count": np.asarray(errors, dtype=np.int16),
        "input_relative_path": uarray(input_paths), "input_text": uarray(input_texts), "output_tail": uarray(output_tails),
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
    for name, series in (("accuracy_Ry", accuracy_series), ("scf_cycle_iterations", cycle_series),
                         ("temperature_K", temperature_series), ("mean_local_moment_Bohr", moment_series),
                         ("printed_constrained_moment_Bohr", constraint_series),
                         ("energy_Ry", energy_series), ("total_force_Ry_Bohr", force_series),
                         ("pressure_kbar", pressure_series)):
        payload[name], payload[f"{name}_offset"] = flatten(series)
    support = [path for path in (ROOT / "README.md", ROOT / "case_manifest.csv", ROOT / "velocity_summary.json",
                                  ROOT / "thermalized_start_summary.json",
                                  ROOT / "magnetic_seed_summary.json",
                                  ROOT / "run_round12_lonestar6.sbatch", ROOT / "check_round12.sh",
                                  ROOT / "plot_round12.sh", ROOT / "generate_round12.py") if path.is_file()]
    payload["support_relative_path"] = uarray([str(path.relative_to(ROOT)) for path in support])
    payload["support_text"] = uarray([path.read_text(errors="replace") for path in support])
    payload["support_sha256"] = uarray([sha256(path) for path in support])
    payload["metadata_json"] = uarray([json.dumps({
        "safe_load": "numpy.load(path, allow_pickle=False)", "full_outputs_included": False,
        "qe_tmp_included": False, "positions_thermalized": True, "velocity_seed": 320400,
        "natoms": NATOMS, "supercell": "4x2x2", "nraise": 20,
        "units": {"energy": "Ry", "force": "Ry/Bohr", "pressure": "kbar", "temperature": "K"},
    }, sort_keys=True)])
    output = args.output.expanduser().resolve()
    np.savez_compressed(output, **payload)
    with np.load(output, allow_pickle=False) as checked:
        assert not any(checked[name].dtype.hasobject for name in checked.files)
        arrays = len(checked.files)
    print(output)
    print(f"MD states={md_states} arrays={arrays} size_bytes={output.stat().st_size} sha256={sha256(output)}")


if __name__ == "__main__":
    main()
