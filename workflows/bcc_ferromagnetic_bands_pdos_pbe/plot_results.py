#!/usr/bin/env python3
"""Plot spin-resolved QE bands and projwfc DOS for both BCC Fe cases."""

from __future__ import annotations

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
CASES = {
    "a_2.50": {
        "title": r"Ferromagnetic BCC Fe, $a=2.50$ $\AA$",
        "pdos": "fe_a250.pdos_tot",
    },
    "a_2.78": {
        "title": r"Ferromagnetic BCC Fe, $a=2.78$ $\AA$",
        "pdos": "fe_a278.pdos_tot",
    },
}
PATH_LABELS = [r"$\Gamma$", "X", "M", r"$\Gamma$", "R", "X", "M", "R"]
FERMI_RE = re.compile(r"the Fermi energy is\s+([-+0-9.Ee]+)\s+ev", re.IGNORECASE)
HSYM_RE = re.compile(r"high-symmetry point:.*x coordinate\s+([-+0-9.Ee]+)", re.IGNORECASE)


def read_fermi_energy(path: Path) -> float:
    matches = FERMI_RE.findall(path.read_text(errors="replace"))
    if not matches:
        raise ValueError(f"Could not find the Fermi energy in {path}")
    return float(matches[-1])


def read_gnu_bands(path: Path) -> list[np.ndarray]:
    bands: list[np.ndarray] = []
    current: list[list[float]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2:
            current.append([float(fields[0]), float(fields[1])])
        elif current:
            bands.append(np.asarray(current))
            current = []
    if current:
        bands.append(np.asarray(current))
    if not bands:
        raise ValueError(f"No bands found in {path}")
    return bands


def read_path_ticks(path: Path) -> list[float]:
    values = [float(value) for value in HSYM_RE.findall(path.read_text(errors="replace"))]
    if len(values) < len(PATH_LABELS):
        raise ValueError(f"Expected at least {len(PATH_LABELS)} path points in {path}, found {len(values)}")
    return values[: len(PATH_LABELS)]


def plot_case(case_dir: Path, title: str, pdos_name: str) -> None:
    fermi_ev = read_fermi_energy(case_dir / "01_scf.out")
    up_bands = read_gnu_bands(case_dir / "bands_up.dat.gnu")
    down_bands = read_gnu_bands(case_dir / "bands_down.dat.gnu")
    ticks = read_path_ticks(case_dir / "03_bands_up.out")
    pdos = np.loadtxt(case_dir / pdos_name, comments="#")
    if pdos.ndim != 2 or pdos.shape[1] < 5:
        raise ValueError(f"Expected E DOSup DOSdw PDOSup PDOSdw in {case_dir / pdos_name}")

    figure, (band_axis, dos_axis) = plt.subplots(
        1,
        2,
        figsize=(10.5, 5.5),
        gridspec_kw={"width_ratios": [2.25, 1.0]},
        constrained_layout=True,
    )

    for index, band in enumerate(up_bands):
        band_axis.plot(
            band[:, 0],
            band[:, 1] - fermi_ev,
            color="#c62828",
            linewidth=1.0,
            label="spin up" if index == 0 else None,
        )
    for index, band in enumerate(down_bands):
        band_axis.plot(
            band[:, 0],
            band[:, 1] - fermi_ev,
            color="#1565c0",
            linewidth=1.0,
            label="spin down" if index == 0 else None,
        )

    for coordinate in ticks:
        band_axis.axvline(coordinate, color="0.75", linewidth=0.7)
    band_axis.axhline(0.0, color="black", linewidth=0.9, linestyle="--")
    band_axis.set_xticks(ticks, PATH_LABELS)
    band_axis.set_xlim(ticks[0], ticks[-1])
    band_axis.set_ylim(-10.0, 10.0)
    band_axis.set_ylabel(r"$E-E_F$ (eV)")
    band_axis.set_title("Electronic bands")
    band_axis.legend(frameon=False, loc="upper right")

    energy = pdos[:, 0] - fermi_ev
    dos_axis.plot(pdos[:, 1], energy, color="#c62828", linewidth=1.3, label="DOS up")
    dos_axis.plot(-pdos[:, 2], energy, color="#1565c0", linewidth=1.3, label="DOS down")
    dos_axis.axhline(0.0, color="black", linewidth=0.9, linestyle="--")
    dos_axis.axvline(0.0, color="0.5", linewidth=0.7)
    dos_axis.set_ylim(-10.0, 10.0)
    dos_axis.set_xlabel("DOS (states/eV)")
    dos_axis.set_yticklabels([])
    dos_axis.set_title("Spin-resolved DOS")
    dos_axis.legend(frameon=False, loc="upper right")

    figure.suptitle(title)
    output = case_dir / "ferromagnetic_bands_and_dos.png"
    figure.savefig(output, dpi=220)
    plt.close(figure)
    print(f"Wrote {output}")


def main() -> None:
    for directory, settings in CASES.items():
        plot_case(ROOT / directory, settings["title"], settings["pdos"])


if __name__ == "__main__":
    main()
