#!/usr/bin/env python3
"""Extract the successful BCC non-collinear trajectory as an HCP reference."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BCC = ROOT / "bcc" / "magnetic-non_coll" / "fe1.out"
OUT = Path(__file__).resolve().parent / "analysis"
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def main() -> None:
    records: list[dict[str, float]] = []
    current: dict[str, float] = {}
    iteration_re = re.compile(r"iteration\s*#\s*(\d+)")
    accuracy_re = re.compile(r"estimated scf accuracy\s*<\s*(" + FLOAT + r")")
    mag_re = re.compile(r"total magnetization\s*=\s*(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")")
    abs_re = re.compile(r"absolute magnetization\s*=\s*(" + FLOAT + r")")
    converged_re = re.compile(r"convergence has been achieved in\s+(\d+) iterations")
    dynamics_re = re.compile(r"Entering Dynamics:\s+iteration\s*=\s*(\d+)")

    for line in BCC.read_text(errors="replace").splitlines():
        if match := iteration_re.search(line):
            current["last_iteration"] = float(match.group(1))
        elif match := accuracy_re.search(line):
            current["accuracy_ry"] = float(match.group(1).replace("D", "E"))
        elif match := mag_re.search(line):
            current["mx"], current["my"], current["mz"] = map(float, match.groups())
        elif match := abs_re.search(line):
            current["absolute"] = float(match.group(1))
        elif match := converged_re.search(line):
            row = dict(current)
            row["scf_iterations"] = float(match.group(1))
            records.append(row)
            current = {}
        elif match := dynamics_re.search(line):
            if records and "md_step" not in records[-1]:
                records[-1]["md_step"] = float(match.group(1))

    # The cycle preceding Dynamics marker 67 converged, but the new positions
    # generated there immediately fail in cdiaghg.  Match the curated BCC
    # dataset definition and retain complete configurations 3--66 only.
    complete = [
        r for r in records
        if all(k in r for k in ("md_step", "scf_iterations", "mx", "my", "mz", "absolute", "accuracy_ry"))
        and r["md_step"] <= 66
    ]
    OUT.mkdir(exist_ok=True)
    with (OUT / "bcc_reference_results.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["md_step", "scf_iterations", "final_accuracy_ry", "mx_muB_cell", "my_muB_cell", "mz_muB_cell", "absolute_muB_cell", "absolute_muB_atom"])
        for r in complete:
            writer.writerow([int(r["md_step"]), int(r["scf_iterations"]), r["accuracy_ry"], r["mx"], r["my"], r["mz"], r["absolute"], r["absolute"] / 54.0])

    x = np.asarray([r["md_step"] for r in complete])
    niter = np.asarray([r["scf_iterations"] for r in complete])
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True, constrained_layout=True)
    for key, label, color in (("mx", r"$M_x$", "#0072B2"), ("my", r"$M_y$", "#D55E00"), ("mz", r"$M_z$", "#009E73")):
        axes[0].plot(x, [r[key] / 54.0 for r in complete], label=label, color=color, lw=1.2)
    axes[0].plot(x, [r["absolute"] / 54.0 for r in complete], label=r"absolute $M$", color="black", ls="--", lw=1.3)
    axes[0].set_ylabel(r"Magnetization ($\mu_B$/atom)")
    axes[0].legend(ncol=4, fontsize=8)
    axes[1].plot(x, niter, color="#7B3294", marker="o", ms=2.5, lw=1.0)
    axes[1].axhline(np.median(niter), color="black", ls="--", lw=1, label=f"median = {np.median(niter):.0f}")
    axes[1].set_ylabel("SCF iterations")
    axes[1].set_xlabel("BCC MD step")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("BCC non-collinear AIMD reference: 64 complete frames (4000 K, 54 atoms)")
    fig.savefig(OUT / "bcc_noncollinear_reference.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Extracted {len(complete)} converged BCC MD steps; median SCF iterations = {np.median(niter):.0f}")


if __name__ == "__main__":
    main()
