#!/usr/bin/env python3
"""Shared helpers for the TDEP workflow scripts."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np

from tdep_phases import get_phase_spec

CUBIC_FOLDER_RE = re.compile(r"tdep_([0-9.]+)_([0-9]+)(?:K)?(?:[-_].+)?$")
HCP_FOLDER_RE = re.compile(r"tdep_a_([0-9.]+)_c_([0-9.]+)_([0-9]+)(?:K)?(?:[-_].+)?$")
CUBIC_NPZ_RE = re.compile(r"([0-9.]+)_([0-9]+K(?:[-_].+)?)$")
HCP_NPZ_RE = re.compile(r"a_([0-9.]+)_c_([0-9.]+)_([0-9]+K(?:[-_].+)?)$")
DEFAULT_COMPARISON_TEMPERATURES = ("4500", "5000", "5500", "6000")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_dataset_dir(phase: str = "bcc") -> Path:
    return get_phase_spec(phase).dataset_dir(repo_root())


def normalize_temperature_label(label: str | int | float) -> str:
    return str(label).removesuffix("K")


def resolve_path(dataset_dir: Path, path: Path) -> Path:
    return path if path.is_absolute() else dataset_dir / path


def folder_match(folder: Path, phase: str = "bcc") -> re.Match[str]:
    pattern = HCP_FOLDER_RE if get_phase_spec(phase).key == "hcp" else CUBIC_FOLDER_RE
    match = pattern.fullmatch(folder.name)
    if not match:
        raise ValueError(f"Could not parse TDEP folder name: {folder}")
    return match


def npz_match(path: Path, phase: str = "bcc") -> re.Match[str]:
    pattern = HCP_NPZ_RE if get_phase_spec(phase).key == "hcp" else CUBIC_NPZ_RE
    match = pattern.fullmatch(path.stem)
    if not match:
        raise ValueError(f"Could not parse NPZ file name: {path}")
    return match


def lattice_parameter_from_folder(folder: Path, phase: str = "bcc") -> float:
    return float(folder_match(folder, phase).group(1))


def lattice_parameters_from_folder(folder: Path, phase: str = "bcc") -> tuple[float, float | None]:
    match = folder_match(folder, phase)
    if get_phase_spec(phase).key == "hcp":
        return float(match.group(1)), float(match.group(2))
    return float(match.group(1)), None


def lattice_parameter_from_npz(path: Path, phase: str = "bcc") -> float:
    return float(npz_match(path, phase).group(1))


def npz_sort_key(path: Path, phase: str = "bcc") -> tuple[float, ...] | tuple[float, str]:
    match = npz_match(path, phase)
    if get_phase_spec(phase).key == "hcp":
        return (float(match.group(1)), float(match.group(2)), path.stem)
    return (float(match.group(1)), path.stem)


def folder_sort_key(folder: Path, phase: str = "bcc") -> tuple[float, ...] | tuple[float, str]:
    match = folder_match(folder, phase)
    if get_phase_spec(phase).key == "hcp":
        return (float(match.group(1)), float(match.group(2)), folder.name)
    return (float(match.group(1)), folder.name)


def folder_case_label(folder: Path, phase: str = "bcc") -> str:
    a_lat, c_lat = lattice_parameters_from_folder(folder, phase)
    if c_lat is None:
        return f"a={a_lat:.2f} A"
    return f"a={a_lat:.2f}, c={c_lat:.2f} A"


def folder_matches_temperature(folder: Path, temperature_label: str, phase: str = "bcc") -> bool:
    match = folder_match(folder, phase)
    if "-disp" in folder.name:
        return False
    temp_group = 3 if get_phase_spec(phase).key == "hcp" else 2
    return match.group(temp_group) == normalize_temperature_label(temperature_label)


def npz_matches_temperature(path: Path, temperature_label: str, phase: str = "bcc") -> bool:
    match = npz_match(path, phase)
    temp_group = 3 if get_phase_spec(phase).key == "hcp" else 2
    return match.group(temp_group).startswith(f"{normalize_temperature_label(temperature_label)}K")


def folder_preference(folder: Path, phase: str = "bcc") -> tuple[int, str]:
    match = folder_match(folder, phase)
    if get_phase_spec(phase).key == "hcp":
        exact_name = f"tdep_a_{match.group(1)}_c_{match.group(2)}_{match.group(3)}K"
    else:
        exact_name = f"tdep_{match.group(1)}_{match.group(2)}K"
    has_suffix = folder.name != exact_name
    return (1 if has_suffix else 0, folder.name)


def prefer_unique_lattice_points(folders: list[Path], phase: str = "bcc") -> list[Path]:
    selected: dict[tuple[float, ...], Path] = {}
    for folder in folders:
        if get_phase_spec(phase).key == "hcp":
            a_lat, c_lat = lattice_parameters_from_folder(folder, phase)
            key = (a_lat, float(c_lat))
        else:
            key = (lattice_parameter_from_folder(folder, phase),)
        current = selected.get(key)
        if current is None or folder_preference(folder, phase) > folder_preference(current, phase):
            selected[key] = folder
    return sorted(selected.values(), key=lambda item: folder_sort_key(item, phase))


def discover_npz_files(dataset_dir: Path, temperature_label: str | None = None, phase: str = "bcc") -> list[Path]:
    pattern = HCP_NPZ_RE if get_phase_spec(phase).key == "hcp" else CUBIC_NPZ_RE
    files = [path for path in dataset_dir.glob("*.npz") if pattern.fullmatch(path.stem)]
    if temperature_label is not None:
        files = [path for path in files if npz_matches_temperature(path, temperature_label, phase)]
    return sorted(files, key=lambda item: npz_sort_key(item, phase))


def discover_bcc_npz_files(dataset_dir: Path, temperature_label: str | None = None) -> list[Path]:
    return discover_npz_files(dataset_dir, temperature_label, phase="bcc")


def default_tdep_folders(dataset_dir: Path, temperature_label: str, phase: str = "bcc") -> list[Path]:
    folders = [
        path
        for path in dataset_dir.glob(f"tdep_*_{normalize_temperature_label(temperature_label)}*")
        if path.is_dir() and folder_matches_temperature(path, temperature_label, phase)
    ]
    return prefer_unique_lattice_points(folders, phase)


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


def read_uc_geometry(path: Path, phase: str = "bcc") -> dict[str, float]:
    lines = path.read_text().splitlines()
    scale = float(lines[1].split()[0])
    cell = np.array([[float(value) for value in lines[index].split()[:3]] for index in range(2, 5)], dtype=float) * scale
    counts = [int(value) for value in lines[6].split()]
    natoms = sum(counts)
    volume = abs(float(np.linalg.det(cell)))
    volume_per_atom = volume / natoms
    spec = get_phase_spec(phase)
    if spec.key == "hcp":
        lattice_a = float(np.linalg.norm(cell[0]))
        lattice_c = float(np.linalg.norm(cell[2]))
        total_volume = volume
        return {
            "lattice_a_A": lattice_a,
            "lattice_c_A": lattice_c,
            "c_over_a": lattice_c / lattice_a,
            "volume_per_atom_A3": volume_per_atom,
            "total_volume_A3": total_volume,
        }
    lattice_a = float(spec.lattice_parameter_from_volume_per_atom(volume_per_atom))
    total_volume = spec.conventional_cell_volume(volume_per_atom)
    return {
        "lattice_a_A": lattice_a,
        "volume_per_atom_A3": volume_per_atom,
        "total_volume_A3": total_volume,
    }


def read_uc_volume_per_atom(path: Path, phase: str = "bcc") -> tuple[float, float, float]:
    geometry = read_uc_geometry(path, phase=phase)
    return float(geometry["lattice_a_A"]), float(geometry["volume_per_atom_A3"]), float(geometry["total_volume_A3"])


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


def dispersion_plot_name(temperature_label: str, phase: str = "bcc") -> Path:
    temperature_label = normalize_temperature_label(temperature_label)
    phase_key = get_phase_spec(phase).key
    if phase_key == "hcp":
        if temperature_label == "5000":
            return Path("phonon_dispersion_overlay_hcp.png")
        return Path(f"phonon_dispersion_overlay_{temperature_label}K_hcp.png")
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
