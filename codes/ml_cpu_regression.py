import argparse
import glob
import os
import pickle
import re
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from ase import Atoms
from matplotlib.lines import Line2D
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

RY_TO_EV = 13.605693009
BOHR_TO_ANG = 0.529177210903


def detect_project_root():
    if "__file__" in globals():
        start_dir = Path(__file__).resolve().parent
    else:
        start_dir = Path.cwd().resolve()

    for candidate in (start_dir, *start_dir.parents):
        if (candidate / "dataset").exists() and (candidate / "IronCoreMD").exists():
            return candidate
        if (candidate / "codes").exists() and (candidate / "dataset").exists():
            return candidate

    return start_dir


PROJECT_ROOT = detect_project_root()


def dataset_pattern_map(dataset_root):
    return {
        "bcc": str(dataset_root / "bcc" / "non-mag" / "*.npz"),
        "fcc": str(dataset_root / "fcc" / "non-mag" / "*.npz"),
        "hcp": str(dataset_root / "hcp" / "*.npz"),
    }


def dataset_match_score(dataset_root):
    patterns = dataset_pattern_map(dataset_root)
    return sum(len(glob.glob(pattern)) for pattern in patterns.values())


def detect_dataset_root(project_root):
    candidates = []
    for candidate in (
        project_root / "dataset",
        project_root.parent / "dataset",
        project_root / "IronCoreMD" / "dataset",
    ):
        candidate = candidate.resolve()
        if candidate.exists() and candidate not in candidates:
            candidates.append(candidate)

    if not candidates:
        return project_root / "dataset"

    scored = sorted(
        ((dataset_match_score(candidate), candidate) for candidate in candidates),
        key=lambda item: item[0],
        reverse=True,
    )
    return scored[0][1]


DATASET_ROOT = detect_dataset_root(PROJECT_ROOT)
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
        description="CPU-only regression baseline for the NPZ/QE datasets without GraphDot or CUDA."
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
        help="Explicit file paths or glob patterns. Overrides the default phase-based dataset selection.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Optional label for the current dataset/model run.",
    )
    parser.add_argument(
        "--output-root",
        default=str(PROJECT_ROOT / "ml-results"),
        help="Root directory where preview files, trained models, and plots are written.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed used for the train/test split.",
    )
    parser.add_argument(
        "--preview-only",
        action="store_true",
        help="Build the feature dataset, write the preview/split files, and stop before fitting the CPU model.",
    )
    parser.add_argument(
        "--model",
        choices=("random_forest", "extra_trees"),
        default="random_forest",
        help="CPU-side scikit-learn regressor to fit.",
    )
    parser.add_argument(
        "--n-estimators",
        type=int,
        default=400,
        help="Number of trees in the ensemble model.",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=12,
        help="Number of nearest-neighbor shells to encode in the structural descriptor.",
    )
    parser.add_argument(
        "--min-samples-leaf",
        type=int,
        default=1,
        help="Minimum samples per leaf in the tree ensemble.",
    )
    parser.add_argument(
        "--max-frames-per-file",
        type=int,
        default=None,
        help="Optional cap on the number of frames read from each input file for fast tests.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Parallel CPU workers for the tree ensemble.",
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
    return slugify_name(f"{phase_label}_{len(file_list)}files_cpu")


def build_output_dir(args, base_name):
    out_dir = Path(os.path.expanduser(args.output_root)) / base_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def build_runtime_args(
    phase=None,
    inputs=None,
    run_name=None,
    output_root=None,
    seed=0,
    preview_only=False,
    model="random_forest",
    n_estimators=400,
    n_neighbors=12,
    min_samples_leaf=1,
    max_frames_per_file=None,
    n_jobs=-1,
):
    return argparse.Namespace(
        phase=list(phase) if phase is not None else ["bcc", "fcc", "hcp"],
        inputs=list(inputs) if inputs is not None else None,
        run_name=run_name,
        output_root=output_root if output_root is not None else str(PROJECT_ROOT / "ml-results"),
        seed=seed,
        preview_only=preview_only,
        model=model,
        n_estimators=n_estimators,
        n_neighbors=n_neighbors,
        min_samples_leaf=min_samples_leaf,
        max_frames_per_file=max_frames_per_file,
        n_jobs=n_jobs,
    )


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


def read_qe_converged_frames(fp):
    energy_pattern = re.compile(
        r"^\s*!\s+total energy\s*=\s*([\-0-9Ee+.]+)\s*Ry",
        re.IGNORECASE,
    )
    atpos_header_pattern = re.compile(
        r"^\s*ATOMIC_POSITIONS\s*(?:\(?([A-Za-z_]+)\)?)?",
        re.IGNORECASE,
    )
    atom_line_pattern = re.compile(
        r"^\s*([A-Z][a-z]?)\s+([\-0-9Ee+.]+)\s+([\-0-9Ee+.]+)\s+([\-0-9Ee+.]+)"
    )
    cell_header_pattern = re.compile(
        r"^\s*CELL_PARAMETERS\s*(?:\(?([A-Za-z_]+)\)?)?",
        re.IGNORECASE,
    )
    crystal_axis_pattern = re.compile(
        r"^\s*a\(\d\)\s*=\s*\(\s*([\-0-9Ee+.]+)\s+([\-0-9Ee+.]+)\s+([\-0-9Ee+.]+)\s*\)"
    )
    alat_pattern = re.compile(
        r"lattice parameter \(alat\)\s*=\s*([\-0-9Ee+.]+)\s*a\.u\.",
        re.IGNORECASE,
    )

    with open(fp, "r", encoding="utf-8", errors="ignore") as handle:
        lines = handle.readlines()

    alat = None
    for line in lines:
        match = alat_pattern.search(line)
        if match:
            alat = float(match.group(1)) * BOHR_TO_ANG
            break

    initial_cell = []
    for index, line in enumerate(lines):
        if "crystal axes:" in line.lower():
            for row in range(index + 1, min(index + 6, len(lines))):
                match = crystal_axis_pattern.match(lines[row].strip())
                if match:
                    initial_cell.append([
                        float(match.group(1)) * alat,
                        float(match.group(2)) * alat,
                        float(match.group(3)) * alat,
                    ])
            break

    initial_cell = np.array(initial_cell, dtype=float) if len(initial_cell) == 3 else None

    last_cell = initial_cell
    last_structure = None
    frames = []
    index = 0

    while index < len(lines):
        line = lines[index]

        match_cell = cell_header_pattern.match(line)
        if match_cell:
            units = match_cell.group(1).lower() if match_cell.group(1) else "angstrom"
            cell = []
            index += 1

            for _ in range(3):
                if index < len(lines):
                    parts = lines[index].split()
                    if len(parts) >= 3:
                        cell.append([float(parts[0]), float(parts[1]), float(parts[2])])
                    index += 1

            if len(cell) == 3:
                cell = np.array(cell, dtype=float)
                if units == "bohr":
                    cell *= BOHR_TO_ANG
                elif units == "alat":
                    cell *= alat
                last_cell = cell
            continue

        match_pos = atpos_header_pattern.match(line)
        if match_pos:
            units = match_pos.group(1).lower() if match_pos.group(1) else "unknown"
            atoms_block = []
            index += 1

            while index < len(lines):
                match_atom = atom_line_pattern.match(lines[index])
                if match_atom:
                    atoms_block.append((
                        match_atom.group(1),
                        float(match_atom.group(2)),
                        float(match_atom.group(3)),
                        float(match_atom.group(4)),
                    ))
                    index += 1
                else:
                    break

            last_structure = {
                "units": units,
                "atoms": atoms_block,
                "cell": last_cell.copy() if last_cell is not None else None,
            }
            continue

        match_energy = energy_pattern.match(line)
        if match_energy and last_structure is not None:
            energy_ev = float(match_energy.group(1)) * RY_TO_EV

            symbols = [atom[0] for atom in last_structure["atoms"]]
            pos = np.array([[atom[1], atom[2], atom[3]] for atom in last_structure["atoms"]], dtype=float)
            units = last_structure["units"]
            cell = last_structure["cell"]

            if units == "bohr":
                pos *= BOHR_TO_ANG
            elif units == "alat":
                pos *= alat
            elif units == "crystal":
                if cell is None:
                    index += 1
                    continue
                pos = pos @ cell

            atoms = Atoms(
                symbols=symbols,
                positions=pos,
                cell=cell,
                pbc=True if cell is not None else False,
            )
            atoms.info["energy"] = energy_ev
            atoms.info["energy_per_atom"] = energy_ev / len(atoms)
            atoms.info["source_file"] = str(fp)
            atoms.info["frame_index"] = len(frames)
            frames.append(atoms)

        index += 1

    return frames


def read_npz_frames(fp, max_frames=None):
    data = np.load(fp, allow_pickle=False)

    positions = np.asarray(data["positions"], dtype=float)
    energy_ry = np.asarray(data["energy_ry"], dtype=float)
    symbols = [str(symbol) for symbol in np.asarray(data["symbols"])]
    fallback_cell = fixed_cell_angstrom(data)

    finite_mask = np.isfinite(positions).all(axis=(1, 2)) & np.isfinite(energy_ry)
    valid_indices = np.flatnonzero(finite_mask)
    if max_frames is not None:
        valid_indices = valid_indices[:max_frames]

    frames = []
    for idx in valid_indices:
        cell_ang = frame_cell_angstrom(data, idx, fallback_cell)
        pos_ang = frame_positions_angstrom(data, idx, cell_ang)
        atoms = Atoms(symbols=symbols, positions=pos_ang, cell=cell_ang, pbc=True)
        energy_ev = float(energy_ry[idx] * RY_TO_EV)
        atoms.info["energy"] = energy_ev
        atoms.info["energy_per_atom"] = energy_ev / len(atoms)
        atoms.info["source_file"] = str(fp)
        atoms.info["frame_index"] = int(idx)
        frames.append(atoms)

    return frames


def read_dataset_frames(fp, max_frames=None):
    if str(fp).lower().endswith(".npz"):
        return read_npz_frames(fp, max_frames=max_frames)
    return read_qe_converged_frames(fp)


def phase_one_hot_columns(phases):
    categories = ["bcc", "fcc", "hcp", "mixed"]
    return [f"phase_{phase}" for phase in categories if phase in phases]


def nearest_neighbor_features(atoms, n_neighbors):
    n_atoms = len(atoms)
    feature_map = {}
    if n_atoms <= 1:
        for index in range(n_neighbors):
            feature_map[f"nn{index + 1}_mean"] = 0.0
            feature_map[f"nn{index + 1}_std"] = 0.0
        return feature_map

    distance_matrix = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distance_matrix, np.inf)
    sorted_distances = np.sort(distance_matrix, axis=1)
    usable_neighbors = min(n_neighbors, n_atoms - 1)
    neighbor_block = sorted_distances[:, :usable_neighbors]

    for index in range(n_neighbors):
        if index < usable_neighbors:
            shell = neighbor_block[:, index]
        else:
            shell = neighbor_block[:, -1]
        feature_map[f"nn{index + 1}_mean"] = float(np.mean(shell))
        feature_map[f"nn{index + 1}_std"] = float(np.std(shell))

    return feature_map


def row_from_atoms(sample_index, atoms, n_neighbors):
    natoms = len(atoms)
    volume = abs(float(atoms.cell.volume)) if atoms.cell.rank == 3 else np.nan
    cellpar = atoms.cell.cellpar() if atoms.cell.rank == 3 else np.array([np.nan] * 6, dtype=float)
    source_file = atoms.info.get("source_file", "unknown")
    phase = infer_phase_name(source_file)

    row = {
        "sample_index": sample_index,
        "Atoms": atoms,
        "E_per_atom": float(atoms.info["energy"] / natoms),
        "natoms": natoms,
        "source_file": source_file,
        "source_name": Path(source_file).name,
        "frame_index": int(atoms.info.get("frame_index", -1)),
        "phase": phase,
        "volume_A3_atom": volume / natoms if np.isfinite(volume) and natoms > 0 else np.nan,
        "density_atoms_A3": natoms / volume if np.isfinite(volume) and volume > 0 else np.nan,
        "a_A": float(cellpar[0]),
        "b_A": float(cellpar[1]),
        "c_A": float(cellpar[2]),
        "alpha_deg": float(cellpar[3]),
        "beta_deg": float(cellpar[4]),
        "gamma_deg": float(cellpar[5]),
    }
    row.update(nearest_neighbor_features(atoms, n_neighbors))
    return row


def build_feature_dataframe(all_atoms, n_neighbors):
    rows = [row_from_atoms(index, atoms, n_neighbors) for index, atoms in enumerate(all_atoms)]
    data = pd.DataFrame(rows)

    phase_dummies = pd.get_dummies(data["phase"], prefix="phase")
    data = pd.concat([data, phase_dummies], axis=1)

    feature_cols = [
        "natoms",
        "volume_A3_atom",
        "density_atoms_A3",
        "a_A",
        "b_A",
        "c_A",
        "alpha_deg",
        "beta_deg",
        "gamma_deg",
    ]
    feature_cols.extend([f"nn{index + 1}_mean" for index in range(n_neighbors)])
    feature_cols.extend([f"nn{index + 1}_std" for index in range(n_neighbors)])
    feature_cols.extend(sorted(phase_dummies.columns))

    return data, feature_cols


def split_dataset(data, seed):
    ntsteps = len(data)
    n_test = ntsteps // 3

    np.random.seed(seed)
    all_idx = np.arange(ntsteps)
    test_sel = np.random.choice(all_idx, n_test, replace=False)
    train_sel = np.setdiff1d(all_idx, test_sel)
    return train_sel, test_sel


def plot_dataset_preview(data, out_dir, base_name):
    split_colors = {"train": "#1f6f5b", "test": "#c24b2a"}
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), dpi=180)

    energy_all = data["E_per_atom"].to_numpy(dtype=float)
    energy_train = data.loc[data["split"] == "train", "E_per_atom"].to_numpy(dtype=float)
    energy_test = data.loc[data["split"] == "test", "E_per_atom"].to_numpy(dtype=float)
    bins = min(60, max(12, len(data) // 20))

    axes[0].hist(energy_all, bins=bins, alpha=0.35, color="#8a8a8a", label="all")
    axes[0].hist(energy_train, bins=bins, alpha=0.65, color=split_colors["train"], label="train")
    axes[0].hist(energy_test, bins=bins, alpha=0.65, color=split_colors["test"], label="test")
    axes[0].set_xlabel("Energy (eV/atom)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Energy Distribution")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False)

    colors = [split_colors[split] for split in data["split"]]
    axes[1].scatter(data["sample_index"], data["E_per_atom"], c=colors, s=18, alpha=0.8)
    axes[1].set_xlabel("Sample index")
    axes[1].set_ylabel("Energy (eV/atom)")
    axes[1].set_title("Train/Test Split")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor=split_colors["train"], label="train", markersize=8),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=split_colors["test"], label="test", markersize=8),
        ],
        frameon=False,
    )

    fig.suptitle(f"Dataset Preview: {base_name}", fontsize=16)
    fig.tight_layout()

    preview_png = out_dir / f"{base_name}_dataset_preview.png"
    fig.savefig(preview_png, bbox_inches="tight")
    plt.close(fig)
    print("saved dataset preview figure to:", preview_png)


def plot_parity(y_true, y_pred, out_dir, base_name, model_name):
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)

    lo = min(y_true.min(), y_pred.min())
    hi = max(y_true.max(), y_pred.max())
    pad = 0.03 * (hi - lo) if hi > lo else 1.0
    lo -= pad
    hi += pad

    fig, ax = plt.subplots(figsize=(9, 9), dpi=150)
    ax.scatter(y_true, y_pred, s=60, alpha=0.75, edgecolors="k", linewidths=0.4)
    ax.plot([lo, hi], [lo, hi], lw=3)

    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useOffset=False))
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useOffset=False))
    ax.yaxis.get_major_formatter().set_scientific(False)

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(axis="x", labelrotation=45)
    for label in ax.get_xticklabels():
        label.set_ha("right")

    ax.grid(True, which="major", alpha=0.25, linewidth=1.2)
    ax.minorticks_on()
    ax.grid(True, which="minor", alpha=0.12, linewidth=0.8)

    ax.set_xlabel("Ground truth energy (eV/atom)", labelpad=20)
    ax.set_ylabel("Predicted energy (eV/atom)", labelpad=20)
    ax.set_title(f"CPU {model_name.replace('_', ' ').title()} Parity Plot", pad=24)

    stats = f"N = {len(y_true)}\nMAE = {mae:.3g} eV/atom\nRMSE = {rmse:.3g} eV/atom\nR2 = {r2:.4f}"
    ax.text(
        0.02,
        0.98,
        stats,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=18,
        fontweight="bold",
        color="black",
        bbox=dict(
            boxstyle="round,pad=0.5",
            facecolor="white",
            edgecolor="black",
            linewidth=1.5,
            alpha=0.95,
        ),
    )

    plt.tight_layout()
    parity_png = out_dir / f"{base_name}_parity_plot.png"
    plt.savefig(parity_png, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.05)
    plt.close()
    print("saved parity plot to:", parity_png)


def plot_feature_importance(model, feature_cols, out_dir, base_name):
    if not hasattr(model, "feature_importances_"):
        return

    importance = pd.DataFrame({
        "feature": feature_cols,
        "importance": model.feature_importances_,
    }).sort_values("importance", ascending=False)

    importance_csv = out_dir / f"{base_name}_feature_importance.csv"
    importance.to_csv(importance_csv, index=False)
    print("saved feature importance csv to:", importance_csv)

    top = importance.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 6), dpi=180)
    ax.barh(top["feature"], top["importance"], color="#2f5597")
    ax.set_xlabel("Importance")
    ax.set_ylabel("Feature")
    ax.set_title("Top Feature Importances")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()

    importance_png = out_dir / f"{base_name}_feature_importance.png"
    fig.savefig(importance_png, bbox_inches="tight")
    plt.close(fig)
    print("saved feature importance plot to:", importance_png)


def save_dataset_preview(data, train_sel, test_sel, out_dir, base_name, feature_cols):
    preview_df = data.drop(columns=["Atoms"], errors="ignore").copy()
    preview_df["split"] = "train"
    preview_df.loc[test_sel, "split"] = "test"

    preview_csv = out_dir / f"{base_name}_dataset_preview.csv"
    train_csv = out_dir / f"{base_name}_train_split.csv"
    test_csv = out_dir / f"{base_name}_test_split.csv"
    features_csv = out_dir / f"{base_name}_feature_columns.txt"

    preview_df.to_csv(preview_csv, index=False)
    preview_df.loc[train_sel].to_csv(train_csv, index=False)
    preview_df.loc[test_sel].to_csv(test_csv, index=False)
    features_csv.write_text("\n".join(feature_cols) + "\n", encoding="utf-8")

    print("saved dataset preview csv to:", preview_csv)
    print("saved training split csv to:", train_csv)
    print("saved test split csv to:", test_csv)
    print("saved feature list to:", features_csv)

    plot_dataset_preview(preview_df, out_dir, base_name)
    return preview_df


def build_model(args):
    common = {
        "n_estimators": args.n_estimators,
        "random_state": args.seed,
        "n_jobs": args.n_jobs,
        "min_samples_leaf": args.min_samples_leaf,
    }
    if args.model == "random_forest":
        return RandomForestRegressor(**common)
    if args.model == "extra_trees":
        return ExtraTreesRegressor(**common)
    raise ValueError(f"Unsupported model={args.model!r}")


def training_cpu_all(file_list, args, base_name, out_dir):
    all_atoms = []

    for file_path in file_list:
        print(f"\nReading: {file_path}")
        atoms_list = read_dataset_frames(file_path, max_frames=args.max_frames_per_file)
        print(f"  usable frames found: {len(atoms_list)}")
        all_atoms.extend(atoms_list)

    if not all_atoms:
        print("No usable frames found in any file.")
        return None

    print("\nTotal usable frames read from all files:", len(all_atoms))

    data, feature_cols = build_feature_dataframe(all_atoms, args.n_neighbors)
    train_sel, test_sel = split_dataset(data, args.seed)
    preview_df = save_dataset_preview(data, train_sel, test_sel, out_dir, base_name, feature_cols)

    if args.preview_only:
        print("Preview-only mode enabled. Skipping CPU model fitting.")
        return None

    train_data = data.iloc[train_sel]
    test_data = data.iloc[test_sel]

    x_train = train_data[feature_cols].to_numpy(dtype=float)
    x_test = test_data[feature_cols].to_numpy(dtype=float)
    y_train = train_data["E_per_atom"].to_numpy(dtype=float)
    y_test = test_data["E_per_atom"].to_numpy(dtype=float)

    model = build_model(args)

    print("N_train =", len(train_data))
    print("N_test  =", len(test_data))
    print("model   =", args.model)
    print("n_estimators =", args.n_estimators)
    print("n_neighbors =", args.n_neighbors)

    start_time = time.time()
    model.fit(x_train, y_train)
    end_time = time.time()

    print(
        "the total time consumption for "
        + str(len(train_data))
        + " training steps is "
        + str((end_time - start_time) / 3600)
        + " hr."
    )

    model_path = out_dir / f"{base_name}_{args.model}_EperAtom_{len(train_data)}.pkl"
    with open(model_path, "wb") as handle:
        pickle.dump(
            {
                "model": model,
                "feature_cols": feature_cols,
                "args": vars(args),
            },
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    print("saved model to:", model_path)

    y_pred = model.predict(x_test)

    prediction_df = preview_df.loc[test_sel, [
        "sample_index",
        "source_file",
        "source_name",
        "frame_index",
        "phase",
        "E_per_atom",
        "split",
    ]].copy()
    prediction_df["predicted_E_per_atom"] = y_pred
    prediction_df["abs_error"] = np.abs(prediction_df["E_per_atom"] - prediction_df["predicted_E_per_atom"])
    prediction_csv = out_dir / f"{base_name}_test_predictions.csv"
    prediction_df.to_csv(prediction_csv, index=False)
    print("saved test predictions csv to:", prediction_csv)

    plot_parity(y_test, y_pred, out_dir, base_name, args.model)
    plot_feature_importance(model, feature_cols, out_dir, base_name)

    return model, y_pred, test_data, "E_per_atom", out_dir, base_name


def execute_cpu_ml(args):
    print("----------------Begin CPU ML---------------------")

    all_files = resolve_input_files(args)
    base_name = build_run_name(args, all_files)
    out_dir = build_output_dir(args, base_name)

    print("Files found:")
    for file_path in all_files:
        print("  ", file_path)

    print("run name:", base_name)
    print("output dir:", out_dir)

    return training_cpu_all(all_files, args, base_name, out_dir)


def run_cpu_ml(
    phase=None,
    inputs=None,
    run_name=None,
    output_root=None,
    seed=0,
    preview_only=False,
    model="random_forest",
    n_estimators=400,
    n_neighbors=12,
    min_samples_leaf=1,
    max_frames_per_file=None,
    n_jobs=-1,
):
    args = build_runtime_args(
        phase=phase,
        inputs=inputs,
        run_name=run_name,
        output_root=output_root,
        seed=seed,
        preview_only=preview_only,
        model=model,
        n_estimators=n_estimators,
        n_neighbors=n_neighbors,
        min_samples_leaf=min_samples_leaf,
        max_frames_per_file=max_frames_per_file,
        n_jobs=n_jobs,
    )
    return execute_cpu_ml(args)


def main(cli_args=None):
    args = parse_args(cli_args)
    return execute_cpu_ml(args)


if __name__ == "__main__":
    main()
