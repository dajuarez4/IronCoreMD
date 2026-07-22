#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


TICK_RE = re.compile(r'\(\s*"([^"]+)"\s+([-+0-9.eE]+)\s*\)')


def read_ticks(path: Path) -> tuple[list[str], np.ndarray]:
    matches = TICK_RE.findall(path.read_text())
    if len(matches) < 2:
        raise ValueError(f"Could not parse high-symmetry ticks from {path}")
    return [label for label, _ in matches], np.asarray([float(value) for _, value in matches])


def read_dispersion(folder: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    data = np.loadtxt(folder / "outfile.dispersion_relations")
    labels, ticks = read_ticks(folder / "outfile.dispersion_relations.gnuplot")
    return data, labels, ticks


def remap_path(x_values: np.ndarray, source_ticks: np.ndarray, target_ticks: np.ndarray) -> np.ndarray:
    if len(source_ticks) != len(target_ticks):
        raise ValueError("Cannot align dispersions with different high-symmetry paths")
    return np.interp(x_values, source_ticks, target_ticks)


def decorate_axis(axis, labels: list[str], ticks: np.ndarray, panel: str, title: str) -> None:
    axis.axhline(0.0, color="#222222", linewidth=0.9)
    for position in ticks:
        axis.axvline(position, color="#B8B8B8", linewidth=0.75, zorder=0)
    axis.set_xlim(ticks[0], ticks[-1])
    axis.set_xticks(ticks, labels)
    axis.set_ylabel("Frequency (THz)")
    axis.set_title(title, loc="left", fontweight="bold", pad=8)
    axis.text(-0.10, 1.02, panel, transform=axis.transAxes, fontsize=18, fontweight="bold")
    axis.grid(axis="y", color="#E2E2E2", linewidth=0.6, alpha=0.75)


def plot_single_phase(axis, folder: Path, color: str, panel: str, title: str) -> None:
    data, labels, ticks = read_dispersion(folder)
    for branch in range(1, data.shape[1]):
        axis.plot(data[:, 0], data[:, branch], color=color, linewidth=1.75, alpha=0.94)
    decorate_axis(axis, labels, ticks, panel, title)
    minimum = min(0.0, float(np.min(data[:, 1:])))
    maximum = float(np.max(data[:, 1:]))
    axis.set_ylim(minimum - 0.4, maximum * 1.08)


def plot_bcc_comparison(axis, magnetic_folder: Path, nonmagnetic_folder: Path) -> None:
    magnetic, labels, ticks = read_dispersion(magnetic_folder)
    nonmagnetic, nonmagnetic_labels, nonmagnetic_ticks = read_dispersion(nonmagnetic_folder)
    if labels != nonmagnetic_labels:
        raise ValueError("BCC magnetic and non-magnetic paths do not match")
    nonmagnetic_x = remap_path(nonmagnetic[:, 0], nonmagnetic_ticks, ticks)
    for branch in range(1, magnetic.shape[1]):
        axis.plot(magnetic[:, 0], magnetic[:, branch], color="#9B1B43", linewidth=1.9)
        axis.plot(
            nonmagnetic_x,
            nonmagnetic[:, branch],
            color="#315B7D",
            linewidth=1.55,
            linestyle="--",
        )
    axis.plot([], [], color="#9B1B43", linewidth=2.2, label="Collinear magnetic")
    axis.plot([], [], color="#315B7D", linewidth=1.8, linestyle="--", label="Non-magnetic")
    decorate_axis(axis, labels, ticks, "a", "BCC at 4000 K: effect of magnetic disorder")
    axis.legend(frameon=False, loc="upper right", ncol=2, handlelength=2.5)
    minimum = min(float(np.min(magnetic[:, 1:])), float(np.min(nonmagnetic[:, 1:])), 0.0)
    maximum = max(float(np.max(magnetic[:, 1:])), float(np.max(nonmagnetic[:, 1:])))
    axis.set_ylim(minimum - 0.3, maximum * 1.10)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Create a poster-ready phase-resolved iron phonon summary.")
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "dataset/figures/iron_phase_phonons_poster.png",
    )
    parser.add_argument("--dpi", type=int, default=400)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 14,
            "axes.titlesize": 16,
            "axes.labelsize": 15,
            "xtick.labelsize": 14,
            "ytick.labelsize": 13,
            "legend.fontsize": 12,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
        }
    )

    figure, axes = plt.subplots(3, 1, figsize=(8.0, 10.6), constrained_layout=True)
    figure.get_layout_engine().set(h_pad=0.06, rect=(0.0, 0.0, 1.0, 0.965))
    figure.suptitle("Phase-Resolved Finite-Temperature Phonons", fontsize=22, fontweight="bold")

    plot_bcc_comparison(
        axes[0],
        root / "dataset/bcc/magnetic-collinear/tdep_simulation",
        root / "dataset/bcc/non-mag/tdep_2.55_4000K",
    )
    plot_single_phase(
        axes[1],
        root / "dataset/fcc/non-mag/tdep_3.05_5000K",
        "#E2761B",
        "b",
        r"FCC at 5000 K: representative $a=3.05$ $\mathrm{\AA}$",
    )
    plot_single_phase(
        axes[2],
        root / "dataset/hcp/tdep_a_2.16_c_3.42_5000K",
        "#1E9D78",
        "c",
        r"HCP at 5000 K: representative $a=2.16$, $c=3.42$ $\mathrm{\AA}$",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf"):
        destination = args.output.with_suffix(suffix)
        figure.savefig(destination, dpi=args.dpi if suffix == ".png" else None, bbox_inches="tight")
        print(f"Wrote {destination}")
    plt.close(figure)


if __name__ == "__main__":
    main()
