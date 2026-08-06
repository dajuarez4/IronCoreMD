#!/usr/bin/env python3
"""Compare the incomplete first SCF cycles of the two HCP non-collinear runs."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "analysis"
RUNS = {
    "try1": {"color": "#0072B2", "linestyle": "-", "label": "try1: m0=0.15, beta=0.10, TF"},
    "try2": {"color": "#D55E00", "linestyle": "--", "label": "try2: m0=0.05, beta=0.01, local-TF"},
}

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_output(path: Path) -> dict[str, np.ndarray]:
    records: list[dict[str, float]] = []
    current: dict[str, float] | None = None
    iteration_re = re.compile(r"iteration\s*#\s*(\d+).*?beta=\s*(" + FLOAT + r")")
    energy_re = re.compile(r"total energy\s*=\s*(" + FLOAT + r")\s*Ry")
    accuracy_re = re.compile(r"estimated scf accuracy\s*<\s*(" + FLOAT + r")\s*Ry")
    mag_re = re.compile(
        r"total magnetization\s*=\s*(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")"
    )
    abs_re = re.compile(r"absolute magnetization\s*=\s*(" + FLOAT + r")")

    for line in path.read_text(errors="replace").splitlines():
        if match := iteration_re.search(line):
            if current is not None and all(k in current for k in ("energy_ry", "accuracy_ry", "mx", "my", "mz", "absolute")):
                records.append(current)
            current = {"iteration": float(match.group(1)), "beta": number(match.group(2))}
        elif current is not None and (match := energy_re.search(line)):
            current["energy_ry"] = number(match.group(1))
        elif current is not None and (match := accuracy_re.search(line)):
            current["accuracy_ry"] = number(match.group(1))
        elif current is not None and (match := mag_re.search(line)):
            current["mx"], current["my"], current["mz"] = map(number, match.groups())
        elif current is not None and (match := abs_re.search(line)):
            current["absolute"] = number(match.group(1))
    if current is not None and all(k in current for k in ("energy_ry", "accuracy_ry", "mx", "my", "mz", "absolute")):
        records.append(current)
    if not records:
        raise RuntimeError(f"No complete SCF iteration records found in {path}")
    keys = ("iteration", "beta", "energy_ry", "accuracy_ry", "mx", "my", "mz", "absolute")
    data = {key: np.asarray([row[key] for row in records]) for key in keys}
    data["net"] = np.sqrt(data["mx"] ** 2 + data["my"] ** 2 + data["mz"] ** 2)
    return data


def dominant_period(values: np.ndarray, start: int = 20) -> float:
    y = values[start:] if values.size > start + 8 else values
    if y.size < 5:
        return math.nan
    x = np.arange(y.size)
    y = y - np.polyval(np.polyfit(x, y, 1), x)
    power = np.abs(np.fft.rfft(y)) ** 2
    power[0] = 0.0
    index = int(np.argmax(power))
    return y.size / index if index else math.nan


def lag1(values: np.ndarray, start: int = 20) -> float:
    y = values[start:] if values.size > start + 8 else values
    if y.size < 3 or np.std(y[:-1]) == 0 or np.std(y[1:]) == 0:
        return math.nan
    return float(np.corrcoef(y[:-1], y[1:])[0, 1])


def write_iteration_csv(all_data: dict[str, dict[str, np.ndarray]]) -> None:
    with (OUT / "scf_iteration_data.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["run", "scf_iteration", "mixing_beta", "energy_ry", "estimated_scf_accuracy_ry", "mx_muB_cell", "my_muB_cell", "mz_muB_cell", "net_magnetization_muB_cell", "absolute_magnetization_muB_cell", "absolute_magnetization_muB_atom"])
        for run, data in all_data.items():
            for i in range(data["iteration"].size):
                writer.writerow([run, int(data["iteration"][i]), data["beta"][i], data["energy_ry"][i], data["accuracy_ry"][i], data["mx"][i], data["my"][i], data["mz"][i], data["net"][i], data["absolute"][i], data["absolute"][i] / 128.0])


def write_parameters() -> None:
    rows = [
        ("Starting magnetization (all 128 types)", "0.15", "0.05"),
        ("mixing_beta", "0.100", "0.010"),
        ("mixing_mode", "TF", "local-TF"),
        ("mixing_ndim", "12", "20"),
        ("diago_thr_init", "QE default", "1.0e-6"),
        ("diago_cg_maxiter", "QE default", "200"),
        ("calculation / restart", "md / from_scratch", "md / from_scratch"),
        ("nat / ntyp", "128 / 128", "128 / 128"),
        ("ecutwfc / ecutrho (Ry)", "100 / 496", "100 / 496"),
        ("occupations / smearing / degauss (Ry)", "smearing / m-v / 0.04", "smearing / m-v / 0.04"),
        ("noncolin / lspinorb / nosym", "true / false / true", "true / false / true"),
        ("constraint / lambda (Ry)", "atomic direction / 0.2", "atomic direction / 0.2"),
        ("conv_thr (Ry) / electron_maxstep", "1.0e-4 / 500", "1.0e-4 / 500"),
        ("diagonalization", "cg", "cg"),
        ("MD nstep / dt (a.u.; fs)", "400 / 20.670; 0.500 fs", "400 / 20.670; 0.500 fs"),
        ("thermostat / T / nraise", "SVR / 5000 K / 20", "SVR / 5000 K / 20"),
        ("potential / wavefunction extrapolation", "second_order / second_order", "second_order / second_order"),
        ("cell vectors (Angstrom)", "(8.56,0,0); (4.28,7.41318,0); (0,0,13.68)", "same"),
        ("k-point grid / shift", "1x1x1 / 0 0 0", "1x1x1 / 0 0 0"),
        ("pseudopotential", "Fe PBE-spn KJPAW PSL 1.0.0", "same"),
    ]
    with (OUT / "input_parameters.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["parameter", "try1", "try2"])
        writer.writerows(rows)
    lines = ["# Input-parameter comparison", "", "| Parameter | try1 | try2 |", "|---|---:|---:|"]
    lines.extend(f"| {p} | {a} | {b} |" for p, a, b in rows)
    lines += ["", "The 128 polar and azimuthal spin-direction angles are identical in both inputs; only their common starting magnitude differs.", ""]
    (OUT / "input_parameters.md").write_text("\n".join(lines))


def style_axis(ax: plt.Axes) -> None:
    ax.grid(True, alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)


def make_plots(all_data: dict[str, dict[str, np.ndarray]]) -> None:
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 220, "font.size": 10})
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True, constrained_layout=True)
    components = (("mx", r"$M_x$"), ("my", r"$M_y$"), ("mz", r"$M_z$"))
    for ax, (key, ylabel) in zip(axes, components):
        for run, data in all_data.items():
            cfg = RUNS[run]
            ax.plot(data["iteration"], data[key], color=cfg["color"], ls=cfg["linestyle"], lw=1.25, label=cfg["label"])
        ax.axhline(0, color="0.35", lw=0.7)
        ax.set_ylabel(ylabel + r" ($\mu_B$/cell)")
        style_axis(ax)
    axes[0].legend(fontsize=8, ncol=2)
    axes[-1].set_xlabel("SCF iteration (first ionic step)")
    fig.suptitle("HCP Fe non-collinear total-magnetization components")
    fig.savefig(OUT / "magnetization_components.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True, constrained_layout=True)
    for run, data in all_data.items():
        cfg = RUNS[run]
        axes[0].plot(data["iteration"], data["absolute"], color=cfg["color"], ls=cfg["linestyle"], lw=1.35, label=cfg["label"])
        axes[1].plot(data["iteration"], data["absolute"] / 128.0, color=cfg["color"], ls=cfg["linestyle"], lw=1.35)
        axes[1].plot(data["iteration"], data["net"] / 128.0, color=cfg["color"], lw=0.9, ls=":", alpha=0.8)
    axes[0].set_ylabel(r"Absolute $M$ ($\mu_B$/cell)")
    axes[1].set_ylabel(r"Magnetization ($\mu_B$/atom)")
    axes[1].set_xlabel("SCF iteration (first ionic step)")
    axes[0].legend(fontsize=8, ncol=2)
    axes[1].legend(["absolute / atom", "net-vector magnitude / atom"], fontsize=8)
    for ax in axes:
        style_axis(ax)
    fig.suptitle("HCP Fe absolute and net magnetization")
    fig.savefig(OUT / "absolute_magnetization.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    for run, data in all_data.items():
        cfg = RUNS[run]
        mask = data["iteration"] >= 20
        ax.plot(data["iteration"][mask], data["absolute"][mask] / 128.0, color=cfg["color"], ls=cfg["linestyle"], lw=1.4, label=cfg["label"])
    ax.set_xlabel("SCF iteration (first ionic step; initialization omitted)")
    ax.set_ylabel(r"Absolute magnetization ($\mu_B$/atom)")
    ax.legend(fontsize=8, ncol=2)
    style_axis(ax)
    fig.suptitle("Late-iteration absolute-magnetization oscillations")
    fig.savefig(OUT / "absolute_magnetization_late_zoom.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True, constrained_layout=True)
    for run, data in all_data.items():
        cfg = RUNS[run]
        axes[0].semilogy(data["iteration"], data["accuracy_ry"], color=cfg["color"], ls=cfg["linestyle"], lw=1.3, label=cfg["label"])
        delta = np.abs(np.diff(data["energy_ry"], prepend=np.nan))
        axes[1].semilogy(data["iteration"][1:], np.maximum(delta[1:], 1e-8), color=cfg["color"], ls=cfg["linestyle"], lw=1.3)
    axes[0].axhline(1e-4, color="black", ls=":", lw=1.0, label=r"conv_thr = $10^{-4}$ Ry")
    axes[0].set_ylabel("Estimated SCF accuracy (Ry)")
    axes[1].set_ylabel(r"$|E_i-E_{i-1}|$ (Ry)")
    axes[1].set_xlabel("SCF iteration (first ionic step)")
    axes[0].legend(fontsize=8, ncol=2)
    for ax in axes:
        style_axis(ax)
    fig.suptitle("SCF convergence comparison")
    fig.savefig(OUT / "scf_convergence.png", bbox_inches="tight")
    plt.close(fig)


def write_analysis(all_data: dict[str, dict[str, np.ndarray]]) -> None:
    lines = ["# SCF and magnetization analysis", "", "Both output files stop during the **first electronic SCF cycle**. No ionic MD step was completed, so every plot uses SCF iteration—not MD time step.", ""]
    lines += ["| Quantity | try1 | try2 |", "|---|---:|---:|"]
    for label, function, fmt in [
        ("Complete printed SCF iterations", lambda d: d["iteration"].size, ".0f"),
        ("Last estimated accuracy (Ry)", lambda d: d["accuracy_ry"][-1], ".4g"),
        ("Minimum estimated accuracy (Ry)", lambda d: np.min(d["accuracy_ry"]), ".4g"),
        ("Last absolute M (muB/cell)", lambda d: d["absolute"][-1], ".3f"),
        ("Last absolute M (muB/atom)", lambda d: d["absolute"][-1] / 128, ".3f"),
        ("Late absolute-M mean (muB/atom)", lambda d: np.mean(d["absolute"][20:]) / 128, ".3f"),
        ("Late absolute-M std (muB/atom)", lambda d: np.std(d["absolute"][20:]) / 128, ".3f"),
        ("Late absolute-M lag-1 correlation", lambda d: lag1(d["absolute"]), ".3f"),
        ("Apparent absolute-M period (iterations)", lambda d: dominant_period(d["absolute"]), ".2f"),
    ]:
        values = [format(function(all_data[run]), fmt) for run in ("try1", "try2")]
        lines.append(f"| {label} | {values[0]} | {values[1]} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "- Neither run is converged: the requested threshold is 1e-4 Ry, while the final printed residuals remain many orders of magnitude larger. The ionic integrator therefore never starts.",
        "- The huge early absolute magnetization is an initialization transient. It collapses quickly, then fluctuates because charge/spin mixing is cycling among substantially different electronic states.",
        "- The component sign changes are not physical 5000 K spin dynamics: they occur inside one unconverged electronic minimization. The `atomic direction` penalty constrains directions but not fixed moment magnitudes, and the 128 independently typed atoms make the SCF problem unusually stiff.",
        "- try1's larger beta updates the density ten times more aggressively and shows stronger overshoot. try2 damps the update, but its residual stalls instead of approaching 1e-4 Ry. Thus the smaller beta reduces aggression but does not cure the underlying instability.",
        "- Any Fourier-derived period is only descriptive, not evidence of a stable physical oscillation: the traces are nonstationary, short, and unconverged.",
        "",
        "## Practical next checks",
        "",
        "Converge a static SCF at the initial geometry before MD; restart from its saved density/wavefunctions; test a staged constraint lambda and starting moment; and compare `plain`/`local-TF` mixing with smaller beta and larger history. Do not use forces, energies, or magnetization from these two outputs as MD data yet.",
        "",
    ]
    (OUT / "analysis.md").write_text("\n".join(lines))


def main() -> None:
    OUT.mkdir(exist_ok=True)
    all_data = {run: parse_output(ROOT / run / "fe1.out") for run in RUNS}
    write_iteration_csv(all_data)
    write_parameters()
    make_plots(all_data)
    write_analysis(all_data)
    print(f"Wrote analysis to {OUT}")


if __name__ == "__main__":
    main()
