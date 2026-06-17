import argparse
import glob
import json
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write
from matplotlib.lines import Line2D

RY_TO_EV = 13.605693009
BOHR_TO_ANG = 0.529177210903
RY_PER_BOHR_TO_EV_PER_ANG = RY_TO_EV / BOHR_TO_ANG


def detect_repo_root():
    if "__file__" in globals():
        return Path(__file__).resolve().parents[1]
    return Path.cwd().resolve()


REPO_ROOT = detect_repo_root()


def dataset_pattern_map(dataset_root):
    return {
        "bcc": str(dataset_root / "bcc" / "non-mag" / "*.npz"),
        "fcc": str(dataset_root / "fcc" / "non-mag" / "*.npz"),
        "hcp": str(dataset_root / "hcp" / "*.npz"),
    }


def dataset_match_score(dataset_root):
    patterns = dataset_pattern_map(dataset_root)
    return sum(len(glob.glob(pattern)) for pattern in patterns.values())


def detect_dataset_root(repo_root):
    candidates = []
    for candidate in (repo_root / "dataset", repo_root.parent / "dataset"):
        candidate = candidate.resolve()
        if candidate.exists() and candidate not in candidates:
            candidates.append(candidate)

    if not candidates:
        return repo_root / "dataset"

    scored = sorted(
        ((dataset_match_score(candidate), candidate) for candidate in candidates),
        key=lambda item: item[0],
        reverse=True,
    )
    return scored[0][1]


DATASET_ROOT = detect_dataset_root(REPO_ROOT)
DEFAULT_DATASET_PATTERNS = dataset_pattern_map(DATASET_ROOT)


def scalar_string(value):
    arr = np.asarray(value)
    return str(arr.item() if arr.shape == () else value)


def slugify_name(text):
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._-")
    return slug or "dataset"


def infer_phase_name(file_path):
    parts = [part.lower() for part in Path(file_path).parts]
    for phase in ("bcc", "fcc", "hcp"):
        if phase in parts:
            return phase
    return "mixed"


def is_notebook_argv():
    argv0 = Path(sys.argv[0]).name.lower() if sys.argv else ""
    return "ipykernel_launcher" in argv0 or "ipykernel" in sys.modules


def parse_args(cli_args=None):
    parser = argparse.ArgumentParser(
        description="Convert QE NPZ archives into extxyz datasets with energies and forces for MLIP training."
    )
    parser.add_argument(
        "--phase",
        nargs="+",
        choices=("bcc", "fcc", "hcp"),
        default=["bcc", "fcc", "hcp"],
        help="Dataset phase groups to include when --inputs is not provided.",
    )
    parser.add_argument(
        "--inputs",
        nargs="*",
        default=None,
        help="Explicit NPZ file paths or glob patterns. Overrides the default phase-based selection.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Optional label for the current exported dataset.",
    )
    parser.add_argument(
        "--output-root",
        default=str(REPO_ROOT / "mlip-data"),
        help="Root directory where extxyz files, CSV metadata, and summaries are written.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed used for the train/val/test split.",
    )
    parser.add_argument(
        "--train-fraction",
        type=float,
        default=0.8,
        help="Fraction of accepted frames assigned to the training set.",
    )
    parser.add_argument(
        "--val-fraction",
        type=float,
        default=0.1,
        help="Fraction of accepted frames assigned to the validation set.",
    )
    parser.add_argument(
        "--max-frames-per-file",
        type=int,
        default=None,
        help="Optional cap on the number of accepted frames read from each input file.",
    )
    if cli_args is None and is_notebook_argv():
        return parser.parse_known_args()[0]
    return parser.parse_args(cli_args)


def resolve_input_files(args):
    patterns = args.inputs if args.inputs else [DEFAULT_DATASET_PATTERNS[phase] for phase in args.phase]

    files = []
    for pattern in patterns:
        expanded = sorted(glob.glob(os.path.expanduser(pattern)))
        if expanded:
            files.extend(expanded)
            continue

        path = os.path.expanduser(pattern)
        if os.path.isfile(path):
            files.append(path)

    files = sorted(dict.fromkeys(files))
    if not files:
        raise FileNotFoundError(
            "No input files matched the requested dataset selection. "
            "Use --inputs with explicit .npz paths or adjust --phase."
        )
    return files


def build_run_name(args, file_list):
    if args.run_name:
        return slugify_name(args.run_name)

    phases = sorted({infer_phase_name(path) for path in file_list})
    phase_label = "_".join(phases)
    if args.inputs and len(file_list) == 1:
        return slugify_name(Path(file_list[0]).stem)
    return slugify_name(f"{phase_label}_{len(file_list)}files_extxyz")


def build_output_dir(args, base_name):
    out_dir = Path(os.path.expanduser(args.output_root)) / base_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def fixed_cell_angstrom(data):
    cell = np.asarray(data["input_cell_parameters"], dtype=float)
    unit = scalar_string(data["input_cell_unit"]).lower()
    if unit in {"angstrom", "ang"}:
        return cell
    if unit in {"bohr", "a.u.", "au"}:
        return cell * BOHR_TO_ANG
    raise ValueError(f"Unsupported input_cell_unit={unit!r}")


def frame_cell_angstrom(data, frame_index, fallback_cell):
    if "cell_parameters" not in data.files:
        return fallback_cell

    frame_cells = np.asarray(data["cell_parameters"], dtype=float)
    if frame_cells.ndim != 3 or frame_index >= len(frame_cells):
        return fallback_cell

    cell = frame_cells[frame_index]
    if not np.isfinite(cell).all():
        return fallback_cell

    units = np.asarray(data["cell_parameters_unit"])
    unit = str(units[frame_index]).lower() if units.ndim > 0 else scalar_string(units).lower()
    if unit in {"none", ""}:
        return fallback_cell
    if unit in {"angstrom", "ang"}:
        return cell
    if unit in {"bohr", "a.u.", "au"}:
        return cell * BOHR_TO_ANG
    raise ValueError(f"Unsupported cell_parameters_unit={unit!r}")


def frame_positions_angstrom(data, frame_index, cell_ang):
    positions = np.asarray(data["positions"], dtype=float)[frame_index]
    units = np.asarray(data["positions_unit"])
    unit = str(units[frame_index]).lower() if units.ndim > 0 else scalar_string(units).lower()

    if unit == "crystal":
        return positions @ cell_ang
    if unit in {"angstrom", "ang"}:
        return positions
    if unit in {"bohr", "a.u.", "au"}:
        return positions * BOHR_TO_ANG

    raise ValueError(f"Unsupported positions_unit={unit!r}")


def optional_scalar(data, key, frame_index, default=np.nan):
    if key not in data.files:
        return default
    values = np.asarray(data[key], dtype=float)
    if values.ndim == 0:
        return float(values)
    if frame_index >= len(values):
        return default
    value = values[frame_index]
    return float(value) if np.isfinite(value) else default


def read_npz_frames(fp, max_frames=None):
    data = np.load(fp, allow_pickle=False)

    positions = np.asarray(data["positions"], dtype=float)
    energy_ry = np.asarray(data["energy_ry"], dtype=float)
    forces_ry_au = np.asarray(data["forces_ry_au"], dtype=float)
    symbols = [str(symbol) for symbol in np.asarray(data["symbols"])]
    fallback_cell = fixed_cell_angstrom(data)

    finite_mask = (
        np.isfinite(positions).all(axis=(1, 2))
        & np.isfinite(energy_ry)
        & np.isfinite(forces_ry_au).all(axis=(1, 2))
    )
    valid_indices = np.flatnonzero(finite_mask)
    if max_frames is not None:
        valid_indices = valid_indices[:max_frames]

    phase = infer_phase_name(fp)
    frames = []
    for idx in valid_indices:
        cell_ang = frame_cell_angstrom(data, idx, fallback_cell)
        if not np.isfinite(cell_ang).all():
            continue

        pos_ang = frame_positions_angstrom(data, idx, cell_ang)
        energy_ev = float(energy_ry[idx] * RY_TO_EV)
        forces_ev_ang = np.asarray(forces_ry_au[idx], dtype=float) * RY_PER_BOHR_TO_EV_PER_ANG
        if not np.isfinite(pos_ang).all() or not np.isfinite(forces_ev_ang).all():
            continue

        atoms = Atoms(symbols=symbols, positions=pos_ang, cell=cell_ang, pbc=True)
        atoms.calc = SinglePointCalculator(atoms, energy=energy_ev, forces=forces_ev_ang)
        atoms.info["source_file"] = str(fp)
        atoms.info["source_name"] = Path(fp).name
        atoms.info["frame_index"] = int(idx)
        atoms.info["phase"] = phase
        atoms.info["config_type"] = phase

        temperature_K = optional_scalar(data, "temperature_K", idx)
        pressure_GPa = optional_scalar(data, "pressure_GPa", idx)
        if not np.isfinite(pressure_GPa):
            pressure_kbar = optional_scalar(data, "pressure_kbar", idx)
            pressure_GPa = pressure_kbar * 0.1 if np.isfinite(pressure_kbar) else np.nan
        time_ps = optional_scalar(data, "time_ps", idx)
        iteration = optional_scalar(data, "iteration", idx)

        if np.isfinite(temperature_K):
            atoms.info["temperature_K"] = temperature_K
        if np.isfinite(pressure_GPa):
            atoms.info["pressure_GPa"] = pressure_GPa
        if np.isfinite(time_ps):
            atoms.info["time_ps"] = time_ps
        if np.isfinite(iteration):
            atoms.info["iteration"] = int(iteration)

        frames.append(atoms)

    return frames


def build_metadata_dataframe(all_atoms):
    rows = []
    for sample_index, atoms in enumerate(all_atoms):
        volume = abs(float(atoms.cell.volume)) if atoms.cell.rank == 3 else np.nan
        natoms = len(atoms)
        rows.append(
            {
                "sample_index": sample_index,
                "source_file": atoms.info.get("source_file", "unknown"),
                "source_name": atoms.info.get("source_name", "unknown"),
                "frame_index": atoms.info.get("frame_index", -1),
                "phase": atoms.info.get("phase", "mixed"),
                "natoms": natoms,
                "energy_eV": float(atoms.get_potential_energy()),
                "energy_eV_atom": float(atoms.get_potential_energy() / natoms),
                "volume_A3": volume,
                "volume_A3_atom": volume / natoms if np.isfinite(volume) else np.nan,
                "temperature_K": atoms.info.get("temperature_K", np.nan),
                "pressure_GPa": atoms.info.get("pressure_GPa", np.nan),
                "time_ps": atoms.info.get("time_ps", np.nan),
            }
        )
    return pd.DataFrame(rows)


def split_indices(num_samples, train_fraction, val_fraction, seed):
    if train_fraction <= 0 or val_fraction < 0 or train_fraction + val_fraction >= 1:
        raise ValueError("Require 0 < train_fraction and 0 <= val_fraction with train_fraction + val_fraction < 1.")

    rng = np.random.default_rng(seed)
    perm = rng.permutation(num_samples)

    n_train = int(num_samples * train_fraction)
    n_val = int(num_samples * val_fraction)
    n_test = num_samples - n_train - n_val
    if n_train == 0 or n_test == 0:
        raise ValueError("Split produced an empty train or test set. Adjust the split fractions.")

    train_idx = perm[:n_train]
    val_idx = perm[n_train:n_train + n_val]
    test_idx = perm[n_train + n_val:]
    return np.sort(train_idx), np.sort(val_idx), np.sort(test_idx)


def plot_dataset_preview(metadata_df, out_dir, base_name):
    split_colors = {"train": "#1f6f5b", "val": "#c98a00", "test": "#c24b2a"}
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), dpi=180)

    energy_all = metadata_df["energy_eV_atom"].to_numpy(dtype=float)
    bins = min(60, max(12, len(metadata_df) // 20))
    axes[0].hist(energy_all, bins=bins, alpha=0.30, color="#8a8a8a", label="all")

    for split in ("train", "val", "test"):
        subset = metadata_df.loc[metadata_df["split"] == split, "energy_eV_atom"].to_numpy(dtype=float)
        if len(subset) == 0:
            continue
        axes[0].hist(subset, bins=bins, alpha=0.65, color=split_colors[split], label=split)

    axes[0].set_xlabel("Energy (eV/atom)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Energy Distribution")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False)

    colors = [split_colors[split] for split in metadata_df["split"]]
    axes[1].scatter(metadata_df["sample_index"], metadata_df["energy_eV_atom"], c=colors, s=18, alpha=0.8)
    axes[1].set_xlabel("Sample index")
    axes[1].set_ylabel("Energy (eV/atom)")
    axes[1].set_title("Train/Val/Test Split")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor=split_colors["train"], label="train", markersize=8),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=split_colors["val"], label="val", markersize=8),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=split_colors["test"], label="test", markersize=8),
        ],
        frameon=False,
    )

    fig.suptitle(f"MLIP Dataset Preview: {base_name}", fontsize=16)
    fig.tight_layout()

    preview_png = out_dir / f"{base_name}_dataset_preview.png"
    fig.savefig(preview_png, bbox_inches="tight")
    plt.close(fig)
    print("saved dataset preview figure to:", preview_png)


def write_extxyz_subset(path, atoms_list):
    write(str(path), atoms_list, format="extxyz")
    print("wrote extxyz:", path)


def build_summary(metadata_df, file_list, args, out_dir, base_name):
    split_counts = metadata_df["split"].value_counts().to_dict()
    source_counts = metadata_df.groupby("source_name").size().sort_index().to_dict()
    phase_counts = metadata_df.groupby("phase").size().sort_index().to_dict()

    summary = {
        "run_name": base_name,
        "output_dir": str(out_dir),
        "input_files": [str(path) for path in file_list],
        "num_input_files": len(file_list),
        "num_frames_total": int(len(metadata_df)),
        "split_counts": {key: int(value) for key, value in split_counts.items()},
        "phase_counts": {key: int(value) for key, value in phase_counts.items()},
        "source_counts": {key: int(value) for key, value in source_counts.items()},
        "train_fraction": args.train_fraction,
        "val_fraction": args.val_fraction,
        "seed": args.seed,
        "max_frames_per_file": args.max_frames_per_file,
        "dataset_root": str(DATASET_ROOT),
        "repo_root": str(REPO_ROOT),
        "units": {
            "positions": "angstrom",
            "cell": "angstrom",
            "energy": "eV",
            "forces": "eV/angstrom",
            "pressure": "GPa",
            "time": "ps",
        },
    }
    return summary


def execute_export(args):
    print("----------------Begin NPZ -> extxyz---------------------")

    file_list = resolve_input_files(args)
    base_name = build_run_name(args, file_list)
    out_dir = build_output_dir(args, base_name)

    print("Files found:")
    for file_path in file_list:
        print("  ", file_path)
    print("run name:", base_name)
    print("output dir:", out_dir)

    all_atoms = []
    for file_path in file_list:
        print(f"\nReading: {file_path}")
        atoms_list = read_npz_frames(file_path, max_frames=args.max_frames_per_file)
        print(f"  usable force-bearing frames found: {len(atoms_list)}")
        all_atoms.extend(atoms_list)

    if not all_atoms:
        raise ValueError("No usable frames found in the selected NPZ archives.")

    metadata_df = build_metadata_dataframe(all_atoms)
    train_idx, val_idx, test_idx = split_indices(len(metadata_df), args.train_fraction, args.val_fraction, args.seed)
    metadata_df["split"] = "test"
    metadata_df.loc[train_idx, "split"] = "train"
    metadata_df.loc[val_idx, "split"] = "val"

    metadata_csv = out_dir / f"{base_name}_frames.csv"
    metadata_df.to_csv(metadata_csv, index=False)
    print("saved frame metadata csv to:", metadata_csv)

    plot_dataset_preview(metadata_df, out_dir, base_name)

    train_atoms = [all_atoms[index] for index in train_idx]
    val_atoms = [all_atoms[index] for index in val_idx]
    test_atoms = [all_atoms[index] for index in test_idx]

    all_path = out_dir / f"{base_name}_all.extxyz"
    train_path = out_dir / f"{base_name}_train.extxyz"
    val_path = out_dir / f"{base_name}_val.extxyz"
    test_path = out_dir / f"{base_name}_test.extxyz"

    write_extxyz_subset(all_path, all_atoms)
    write_extxyz_subset(train_path, train_atoms)
    if val_atoms:
        write_extxyz_subset(val_path, val_atoms)
    else:
        print("validation split is empty; skipping val.extxyz")
    write_extxyz_subset(test_path, test_atoms)

    summary = build_summary(metadata_df, file_list, args, out_dir, base_name)
    summary_json = out_dir / f"{base_name}_summary.json"
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("saved dataset summary to:", summary_json)

    print("\nSplit sizes:")
    print("  train =", len(train_atoms))
    print("  val   =", len(val_atoms))
    print("  test  =", len(test_atoms))
    print("  total =", len(all_atoms))

    return {
        "all_extxyz": all_path,
        "train_extxyz": train_path,
        "val_extxyz": val_path if val_atoms else None,
        "test_extxyz": test_path,
        "metadata_csv": metadata_csv,
        "summary_json": summary_json,
        "preview_png": out_dir / f"{base_name}_dataset_preview.png",
    }


def main(cli_args=None):
    args = parse_args(cli_args)
    return execute_export(args)


if __name__ == "__main__":
    main()
