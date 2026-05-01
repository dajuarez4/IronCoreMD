#!/usr/bin/env python3
"""Plot combined BCC/HCP Fe volume vs pressure at 5000 K."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
from pathlib import Path
from types import ModuleType

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot BCC and HCP volume-pressure points with BM fits to AIMD mean pressures."
    )
    parser.add_argument(
        "--bcc-script",
        type=Path,
        default=Path("bcc/non-mag/codes/plot_volume_vs_pressure.py"),
        help="BCC plot_volume_vs_pressure.py script path, relative to dataset root.",
    )
    parser.add_argument(
        "--hcp-script",
        type=Path,
        default=Path("hcp/codes/plot_volume_vs_pressure.py"),
        help="HCP plot_volume_vs_pressure.py script path, relative to dataset root.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("volume_vs_pressure_5000K_bcc_hcp_compare.png"),
        help="Output PNG path, relative to dataset root.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("volume_vs_pressure_5000K_bcc_hcp_compare.csv"),
        help="Output CSV path, relative to dataset root.",
    )
    return parser.parse_args()


def resolve_path(dataset_dir: Path, path: Path) -> Path:
    return path if path.is_absolute() else dataset_dir / path


def load_module(path: Path, module_name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def collect_phase(
    script_path: Path,
    module_name: str,
) -> tuple[list[dict[str, float | str]], np.ndarray]:
    module = load_module(script_path, module_name)
    data_dir = script_path.parent.parent
    free_energy_csv = data_dir / "free_energy_vs_volume.csv"
    rows, _, dense_md_curve = module.collect_rows(data_dir, free_energy_csv)
    return rows, dense_md_curve


def write_csv(
    path: Path,
    bcc_rows: list[dict[str, float | str]],
    hcp_rows: list[dict[str, float | str]],
) -> None:
    fieldnames = [
        "phase",
        "folder",
        "total_volume_A3",
        "volume_per_atom_A3",
        "mean_md_pressure_GPa",
        "std_md_pressure_GPa",
        "pressure_from_md_bm_fit_GPa",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for phase, rows in (("BCC", bcc_rows), ("HCP", hcp_rows)):
            for row in rows:
                writer.writerow(
                    {
                        "phase": phase,
                        "folder": row["folder"],
                        "total_volume_A3": float(row["total_volume_A3"]),
                        "volume_per_atom_A3": float(row["volume_per_atom_A3"]),
                        "mean_md_pressure_GPa": float(row["mean_md_pressure_GPa"]),
                        "std_md_pressure_GPa": float(row["std_md_pressure_GPa"]),
                        "pressure_from_md_bm_fit_GPa": float(row["pressure_from_md_bm_fit_GPa"]),
                    }
                )


def plot_phase(
    ax: plt.Axes,
    rows: list[dict[str, float | str]],
    dense_md_curve: np.ndarray,
    *,
    point_face: str,
    point_edge: str,
    box_color: str,
    line_color: str,
    label: str,
) -> None:
    md_pressure = np.array([float(row["mean_md_pressure_GPa"]) for row in rows], dtype=float)
    md_std = np.array([float(row["std_md_pressure_GPa"]) for row in rows], dtype=float)
    volume = np.array([float(row["total_volume_A3"]) for row in rows], dtype=float)

    order_md = np.argsort(md_pressure)
    order_curve = np.argsort(dense_md_curve[:, 0])
    sorted_volume = volume[order_md]

    if sorted_volume.size > 1:
        box_height = 0.35 * float(np.median(np.abs(np.diff(sorted_volume))))
    else:
        box_height = 0.08

    ax.plot(
        dense_md_curve[order_curve, 0],
        dense_md_curve[order_curve, 1],
        ":",
        color=line_color,
        linewidth=1.6,
        label="_nolegend_",
    )
    ax.barh(
        sorted_volume,
        2.0 * md_std[order_md],
        left=md_pressure[order_md] - md_std[order_md],
        height=box_height,
        color=box_color,
        alpha=0.58,
        edgecolor="none",
        label="_nolegend_",
    )
    ax.plot(
        md_pressure[order_md],
        sorted_volume,
        "s",
        ms=7.4,
        mfc=point_face,
        mec=point_edge,
        mew=0.9,
        label=label,
    )


def plot(
    path: Path,
    bcc_rows: list[dict[str, float | str]],
    bcc_dense_md_curve: np.ndarray,
    hcp_rows: list[dict[str, float | str]],
    hcp_dense_md_curve: np.ndarray,
) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 6.2), constrained_layout=True)

    plot_phase(
        ax,
        bcc_rows,
        bcc_dense_md_curve,
        point_face="#1f5aa6",
        point_edge="#123b73",
        box_color="#d95c5c",
        line_color="#1f5aa6",
        label="Iron body centered cubic at 5000K",
    )
    plot_phase(
        ax,
        hcp_rows,
        hcp_dense_md_curve,
        point_face="#d8892b",
        point_edge="#b76c12",
        box_color="#f2c48f",
        line_color="#b76c12",
        label="Iron hexagonal close packed at 5000K",
    )

    ax.set_xlabel("Pressure (GPa)", fontsize=20)
    ax.set_ylabel(r"Volume ($\AA^3$)", fontsize=20)
    ax.set_xlim(150.0, 500.0)
    ax.set_ylim(12.0, 16.0)
    ax.tick_params(axis="both", labelsize=16)
    ax.grid(True, color="#d9c7f0", alpha=0.45)
    ax.legend(frameon=False, loc="upper right", fontsize=14)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    dataset_dir = Path(__file__).resolve().parent.parent

    bcc_script = resolve_path(dataset_dir, args.bcc_script)
    hcp_script = resolve_path(dataset_dir, args.hcp_script)
    output = resolve_path(dataset_dir, args.output)
    csv_path = resolve_path(dataset_dir, args.csv)

    bcc_rows, bcc_dense_md_curve = collect_phase(bcc_script, "bcc_plot_volume_vs_pressure")
    hcp_rows, hcp_dense_md_curve = collect_phase(hcp_script, "hcp_plot_volume_vs_pressure")

    write_csv(csv_path, bcc_rows, hcp_rows)
    plot(output, bcc_rows, bcc_dense_md_curve, hcp_rows, hcp_dense_md_curve)
    print(f"Wrote {csv_path}")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
