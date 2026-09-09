#!/usr/bin/env python3
"""Summarize Round-20 completion and MD quality from the light NPZ."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def source_series(data: np.lib.npyio.NpzFile, name: str, source: int) -> np.ndarray:
    start, stop = data[f"{name}_offset"][source:source + 2]
    return data[name][start:stop]


def main() -> None:
    archive = ROOT / "round20_results_light.npz"
    rows = []
    with np.load(archive, allow_pickle=False) as data:
        for case_index, case_name in enumerate(data["case_name"]):
            source = case_index * 4 + 3
            temperature = source_series(data, "temperature_K", source)
            moment = source_series(data, "mean_local_moment_Bohr", source)
            energy = source_series(data, "energy_Ry", source)
            force = source_series(data, "total_force_Ry_Bohr", source)
            pressure = source_series(data, "pressure_kbar", source)
            steps = int(data["md_steps"][source])
            rows.append({
                "case": str(case_name), "status": str(data["status"][source]), "steps": steps,
                "last_temperature_K": float(temperature[-1]),
                "mean_temperature_K": float(temperature.mean()),
                "mean_last50_temperature_K": float(temperature[-50:].mean()),
                "mean_last50_local_moment_Bohr": float(moment[-50:].mean()),
                "mean_last50_energy_Ry": float(energy[-50:].mean()),
                "mean_last50_total_force_Ry_Bohr": float(force[-50:].mean()),
                "mean_last50_pressure_kbar": float(pressure[-50:].mean()),
                "scf_iterations": int(data["iteration_count"][source]),
                "mean_scf_iterations_per_step": float(data["iteration_count"][source] / max(steps, 1)),
                "force_correction_warnings": int(data["force_correction_warning_count"][source]),
                "c_bands_warnings": int(data["c_bands_warning_count"][source]),
                "last_accuracy_Ry": float(data["last_accuracy_Ry"][source]),
            })

    with (ROOT / "round20_quality.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    completed = [row for row in rows if row["status"] == "DONE"]
    summary = {
        "archive": archive.name,
        "replicas": len(rows), "completed_replicas": len(completed),
        "failed_replicas": [row["case"] for row in rows if row["status"] == "FAILED"],
        "usable_md_steps": sum(row["steps"] for row in rows),
        "completed_mean_last50_temperature_K": float(np.mean([row["mean_last50_temperature_K"] for row in completed])),
        "completed_std_last50_temperature_K_across_replicas": float(np.std([row["mean_last50_temperature_K"] for row in completed], ddof=1)),
        "completed_mean_last50_local_moment_Bohr": float(np.mean([row["mean_last50_local_moment_Bohr"] for row in completed])),
        "force_correction_warnings": sum(row["force_correction_warnings"] for row in rows),
        "assessment": "magnetically stable but electronic/force convergence unsuitable for longer production",
        "recommended_next_step": "controlled stabilization pilots on the failed dlm01 seed",
    }
    (ROOT / "round20_quality_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
