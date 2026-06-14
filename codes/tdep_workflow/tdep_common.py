#!/usr/bin/env python3
"""Shared helpers for the TDEP workflow scripts."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np

from tdep_phases import get_phase_spec

FOLDER_RE = re.compile(r"tdep_([0-9.]+)_([0-9]+)(?:K)?(?:[-_].+)?$")
NPZ_RE = re.compile(r"([0-9.]+)_([0-9]+K(?:[-_].+)?)$")
DEFAULT_COMPARISON_TEMPERATURES = ("4500", "5000", "5500")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_dataset_dir(phase: str = "bcc") -> Path:
    return get_phase_spec(phase).dataset_dir(repo_root())


def normalize_temperature_label(label: str | int | float) -> str:
    return str(label).removesuffix("K")


def resolve_path(dataset_dir: Path, path: Path) -> Path:
    return path if path.is_absolute() else dataset_dir / path


def folder_match(folder: Path) -> re.Match[str]:
    match = FOLDER_RE.fullmatch(folder.name)
    if not match:
        raise ValueError(f"Could not parse TDEP folder name: {folder}")
    return match


def npz_match(path: Path) -> re.Match[str]:
    match = NPZ_RE.fullmatch(path.stem)
    if not match:
        raise ValueError(f"Could not parse NPZ file name: {path}")
    return match


def lattice_parameter_from_folder(folder: Path) -> float:
    return float(folder_match(folder).group(1))


def lattice_parameter_from_npz(path: Path) -> float:
    return float(npz_match(path).group(1))


def folder_matches_temperature(folder: Path, temperature_label: str) -> bool:
    match = folder_match(folder)
    if "-disp" in folder.name:
        return False
    return match.group(2) == normalize_temperature_label(temperature_label)


def npz_matches_temperature(path: Path, temperature_label: str) -> bool:
    match = npz_match(path)
    return match.group(2).startswith(f"{normalize_temperature_label(temperature_label)}K")


def folder_preference(folder: Path) -> tuple[int, str]:
    match = folder_match(folder)
    exact_name = f"tdep_{match.group(1)}_{match.group(2)}K"
    has_suffix = folder.name != exact_name
    return (1 if has_suffix else 0, folder.name)


def prefer_unique_lattice_points(folders: list[Path]) -> list[Path]:
    selected: dict[float, Path] = {}
    for folder in folders:
        lattice = lattice_parameter_from_folder(folder)
        current = selected.get(lattice)
        if current is None or folder_preference(folder) > folder_preference(current):
            selected[lattice] = folder
    return sorted(selected.values(), key=lambda item: (lattice_parameter_from_folder(item), item.name))


def discover_npz_files(dataset_dir: Path, temperature_label: str | None = None) -> list[Path]:
    files = [path for path in dataset_dir.glob("*.npz") if NPZ_RE.fullmatch(path.stem)]
    if temperature_label is not None:
        files = [path for path in files if npz_matches_temperature(path, temperature_label)]
    return sorted(files, key=lambda item: (lattice_parameter_from_npz(item), item.stem))


def discover_bcc_npz_files(dataset_dir: Path, temperature_label: str | None = None) -> list[Path]:
    return discover_npz_files(dataset_dir, temperature_label)


def default_tdep_folders(dataset_dir: Path, temperature_label: str) -> list[Path]:
    folders = [
        path
        for path in dataset_dir.glob(f"tdep_*_{normalize_temperature_label(temperature_label)}*")
        if path.is_dir() and folder_matches_temperature(path, temperature_label)
    ]
    return prefer_unique_lattice_points(folders)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def read_free_energy(path: Path) -> tuple[float, float, float, float]:
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 4:
            return tuple(float(fields[index]) for index in range(4))
    raise ValueError(f"No free-energy row found in {path}")


def read_u0_second_order(path: Path) -> float:
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) >= 2:
            return float(fields[1])
    raise ValueError(f"No U0 row found in {path}")


def classify_free_energy(temperature: float, f_vib: float, entropy: float, cv: float, u0: float) -> str:
    values = np.array([temperature, f_vib, entropy, cv, u0, u0 + f_vib], dtype=float)
    if not bool(np.all(np.isfinite(values))):
        return "bad"
    if abs(f_vib) > 1.0e6:
        return "bad_free_energy_check_imaginary_modes"
    return "ok"


def read_uc_volume_per_atom(path: Path, phase: str = "bcc") -> tuple[float, float, float]:
    lines = path.read_text().splitlines()
    scale = float(lines[1].split()[0])
    cell = np.array([[float(value) for value in lines[index].split()[:3]] for index in range(2, 5)], dtype=float) * scale
    counts = [int(value) for value in lines[6].split()]
    natoms = sum(counts)
    volume = abs(float(np.linalg.det(cell)))
    volume_per_atom = volume / natoms
    spec = get_phase_spec(phase)
    lattice_a = float(spec.lattice_parameter_from_volume_per_atom(volume_per_atom))
    total_volume = spec.conventional_cell_volume(volume_per_atom)
    return lattice_a, volume_per_atom, total_volume


def source_npz_for_folder(folder: Path) -> Path:
    source_file = folder / "source_npz.txt"
    if source_file.exists():
        for line in source_file.read_text().splitlines():
            if line.startswith("source_npz:"):
                return Path(line.split(":", 1)[1].strip())
    return folder.parent / f"{folder.name.removeprefix('tdep_')}.npz"


def free_energy_csv_name(temperature_label: str) -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    if temperature_label == "5000":
        return Path("free_energy_vs_volume.csv")
    return Path(f"free_energy_vs_volume_{temperature_label}K.csv")


def free_energy_plot_name(temperature_label: str) -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    if temperature_label == "5000":
        return Path("free_energy_vs_volume.png")
    return Path(f"free_energy_vs_volume_{temperature_label}K.png")


def relative_free_energy_plot_name(temperature_label: str) -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    if temperature_label == "5000":
        return Path("relative_free_energy_vs_volume.png")
    return Path(f"relative_free_energy_vs_volume_{temperature_label}K.png")


def free_energy_lattice_plot_name(temperature_label: str) -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    if temperature_label == "5000":
        return Path("free_energy_vs_lattice.png")
    return Path(f"free_energy_vs_lattice_{temperature_label}K.png")


def relative_free_energy_lattice_plot_name(temperature_label: str) -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    if temperature_label == "5000":
        return Path("relative_free_energy_vs_lattice.png")
    return Path(f"relative_free_energy_vs_lattice_{temperature_label}K.png")


def pressure_plot_name(temperature_label: str, phase: str = "bcc") -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    return Path(f"volume_vs_pressure_{temperature_label}K_{get_phase_spec(phase).key}.png")


def pressure_csv_name(temperature_label: str, phase: str = "bcc") -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    return Path(f"volume_vs_pressure_{temperature_label}K_{get_phase_spec(phase).key}.csv")


def pressure_eos_plot_name(temperature_label: str, phase: str = "bcc") -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    return Path(f"volume_vs_pressure_{temperature_label}K_{get_phase_spec(phase).key}_eos_std.png")


def dispersion_plot_name(temperature_label: str) -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    if temperature_label == "5000":
        return Path("phonon_dispersion_overlay.png")
    return Path(f"phonon_dispersion_overlay_{temperature_label}K.png")


def comparison_output_name(prefix: str, temperatures: list[str], suffix: str) -> Path:
    labels = "_".join(f"{normalize_temperature_label(temp)}K" for temp in temperatures)
    return Path(f"{prefix}_{labels}{suffix}")


def discover_temperature_series(dataset_dir: Path, phase: str = "bcc") -> list[str]:
    available = [
        temperature
        for temperature in DEFAULT_COMPARISON_TEMPERATURES
        if resolve_path(dataset_dir, free_energy_csv_name(temperature)).exists()
        and resolve_path(dataset_dir, pressure_csv_name(temperature, phase)).exists()
    ]
    if not available:
        raise FileNotFoundError(f"No free-energy/pressure CSV pairs found in {dataset_dir}")
    return available


def find_tdep_root(explicit: Path | None = None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit.resolve())
    root = repo_root()
    candidates.extend(
        [
            root / "tdep" / "build" / "src",
            root.parent / "tdep" / "build" / "src",
        ]
    )
    for candidate in candidates:
        if (candidate / "extract_forceconstants" / "extract_forceconstants").exists() and (
            candidate / "phonon_dispersion_relations" / "phonon_dispersion_relations"
        ).exists():
            return candidate
    searched = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"Could not locate the TDEP build/src directory. Searched: {searched}")
