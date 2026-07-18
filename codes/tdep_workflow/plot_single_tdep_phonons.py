#!/usr/bin/env python3
"""Plot one TDEP phonon dispersion and total DOS and write numerical diagnostics."""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--summary", type=Path, default=None)
    parser.add_argument("--title", default="BCC Fe Noncollinear TDEP Phonons")
    return parser.parse_args()


def read_ticks(path: Path) -> tuple[list[str], np.ndarray]:
    pattern = re.compile(r'\(\s*"([^"]+)"\s+([0-9Ee+\-.]+)\s*\)')
    matches = [pattern.search(line) for line in path.read_text().splitlines()]
    pairs = [(match.group(1), float(match.group(2))) for match in matches if match]
    if not pairs:
        raise ValueError(f"No path ticks found in {path}")
    return [label for label, _ in pairs], np.asarray([position for _, position in pairs], dtype=float)


def main() -> int:
    args = parse_args()
    folder = args.folder.resolve()
    output = (args.output or folder / "phonon_dispersion_and_dos.png").resolve()
    summary_path = (args.summary or folder / "phonon_summary.json").resolve()

    dispersion = np.loadtxt(folder / "outfile.dispersion_relations")
    dos = np.loadtxt(folder / "outfile.phonon_dos")
    labels, ticks = read_ticks(folder / "outfile.dispersion_relations.gnuplot")
    frequencies = dispersion[:, 1:]
    dos_frequency = dos[:, 0]
    total_dos = dos[:, 1]

    plt.rcParams.update({"font.family": "serif", "font.size": 14})
    fig, (axis, dos_axis) = plt.subplots(
        1,
        2,
        figsize=(11.5, 6.8),
        sharey=True,
        gridspec_kw={"width_ratios": [4.8, 1.4]},
        constrained_layout=True,
    )
    for branch in range(frequencies.shape[1]):
        axis.plot(dispersion[:, 0], frequencies[:, branch], color="#7b3294", linewidth=1.8)
    for tick in ticks:
        axis.axvline(tick, color="#888888", linewidth=0.8, alpha=0.55)
    axis.axhline(0.0, color="black", linewidth=0.9)
    axis.set_xlim(ticks[0], ticks[-1])
    axis.set_xticks(ticks)
    axis.set_xticklabels(labels)
    axis.set_ylabel("Frequency (THz)")
    axis.set_title(args.title)
    axis.grid(axis="y", alpha=0.22)

    dos_axis.plot(total_dos, dos_frequency, color="#008837", linewidth=1.8)
    dos_axis.axhline(0.0, color="black", linewidth=0.9)
    dos_axis.set_xlabel("DOS")
    dos_axis.set_title("Total DOS")
    dos_axis.grid(axis="y", alpha=0.22)
    dos_axis.tick_params(axis="y", left=False, labelleft=False)

    minimum = min(float(np.min(frequencies)), float(np.min(dos_frequency)))
    maximum = max(float(np.max(frequencies)), float(np.max(dos_frequency)))
    axis.set_ylim(min(0.0, minimum * 1.08), maximum * 1.04)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)

    summary = {
        "folder": str(folder),
        "minimum_frequency_THz": float(np.min(frequencies)),
        "maximum_frequency_THz": float(np.max(frequencies)),
        "negative_dispersion_points": int(np.count_nonzero(frequencies < -1.0e-8)),
        "minimum_dos_frequency_THz": float(np.min(dos_frequency)),
        "maximum_dos_frequency_THz": float(np.max(dos_frequency)),
        "path_labels": labels,
        "path_ticks": ticks.tolist(),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
