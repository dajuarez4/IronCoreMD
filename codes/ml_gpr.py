import argparse
import time
import os
import re
import glob
import shutil
import sys
import numpy as np
import pandas as pd
from ase import Atoms
from graphdot import Graph
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from graphdot.graph.adjacency import AtomicAdjacency
from graphdot.model.gaussian_process import GaussianProcessRegressor
from graphdot.kernel.fix import Normalization
from graphdot.kernel.molecular import Tang2019MolecularKernel as MolecularKernel
from pathlib import Path

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
DEFAULT_DATASET_PATTERNS = {
    "bcc": str(PROJECT_ROOT / "dataset" / "bcc" / "non-mag" / "*.npz"),
    "fcc": str(PROJECT_ROOT / "dataset" / "fcc" / "non-mag" / "*.npz"),
    "hcp": str(PROJECT_ROOT / "dataset" / "hcp" / "*.npz"),
}


def ensure_cuda_compiler_available():
    nvcc_path = shutil.which("nvcc")
    if nvcc_path:
        return nvcc_path

    raise RuntimeError(
        "GraphDot training requires CUDA and the 'nvcc' compiler, but 'nvcc' was not found in PATH. "
        "Load a CUDA module or activate an environment where nvcc is available, then rerun the training step. "
        "If you only want the dataset split/preview files, use preview_only=True or the plot_ml_dataset_split.py workflow."
    )


def plot(mu, test_data, target, out_dir, base_name):
    y_true =  np.asarray(test_data[target])
    y_pred =  np.asarray(mu)

    mae  = np.mean(np.abs(y_true - y_pred))
    rmse = np.sqrt(np.mean((y_true - y_pred)**2))

    lo = min(y_true.min(), y_pred.min())
    hi = max(y_true.max(), y_pred.max())
    pad = 0.03 * (hi - lo) if hi > lo else 1.0
    lo -= pad
    hi += pad

    fig, ax = plt.subplots(figsize=(9, 9), dpi=150)

    ax.scatter(y_true, y_pred, s=60, alpha=0.75, edgecolors="k", linewidths=0.4)
    ax.plot([lo, hi], [lo, hi], lw=3)

    ax.ticklabel_format(axis='x', style='plain', useOffset=False)
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useOffset=False))
    ax.xaxis.get_major_formatter().set_scientific(False)

    ax.ticklabel_format(axis='y', style='plain', useOffset=False)
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
    ax.set_title(f"Gausian Proccess Regression Parity Plot", pad=24)

    stats = f"N = {len(y_true)}\nMAE = {mae:.3g} eV/atom\nRMSE = {rmse:.3g} eV/atom"
    ax.text(
        0.02, 0.98, stats,
        transform=ax.transAxes,
        va="top", ha="left",
        fontsize=20,
        fontweight="bold",
        color="black",
        bbox=dict(
            boxstyle="round,pad=0.5",
            facecolor="white",
            edgecolor="black",
            linewidth=1.5,
            alpha=0.95
        )
    )

    plt.tight_layout()
    plt.savefig(
        os.path.join(out_dir, f"{base_name}_parity_plot.png"),
        dpi=300,
        transparent=True,
        bbox_inches="tight",
        pad_inches=0.05
    )
    plt.close()


def read_qe_converged_frames(fp):
    energy_pattern = re.compile(
        r'^\s*!\s+total energy\s*=\s*([\-0-9Ee+.]+)\s*Ry',
        re.IGNORECASE
    )
    atpos_header_pattern = re.compile(
        r'^\s*ATOMIC_POSITIONS\s*(?:\(?([A-Za-z_]+)\)?)?',
        re.IGNORECASE
    )
    atom_line_pattern = re.compile(
        r'^\s*([A-Z][a-z]?)\s+([\-0-9Ee+.]+)\s+([\-0-9Ee+.]+)\s+([\-0-9Ee+.]+)'
    )
    cell_header_pattern = re.compile(
        r'^\s*CELL_PARAMETERS\s*(?:\(?([A-Za-z_]+)\)?)?',
        re.IGNORECASE
    )
    crystal_axis_pattern = re.compile(
        r'^\s*a\(\d\)\s*=\s*\(\s*([\-0-9Ee+.]+)\s+([\-0-9Ee+.]+)\s+([\-0-9Ee+.]+)\s*\)'
    )
    alat_pattern = re.compile(
        r'lattice parameter \(alat\)\s*=\s*([\-0-9Ee+.]+)\s*a\.u\.',
        re.IGNORECASE
    )

    with open(fp, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    alat = None
    for line in lines:
        m = alat_pattern.search(line)
        if m:
            alat = float(m.group(1)) * BOHR_TO_ANG
            break

    initial_cell = []
    for i, line in enumerate(lines):
        if "crystal axes:" in line.lower():
            for j in range(i + 1, min(i + 6, len(lines))):
                m = crystal_axis_pattern.match(lines[j].strip())
                if m:
                    initial_cell.append([
                        float(m.group(1)) * alat,
                        float(m.group(2)) * alat,
                        float(m.group(3)) * alat
                    ])
            break

    initial_cell = np.array(initial_cell, dtype=float) if len(initial_cell) == 3 else None

    last_cell = initial_cell
    last_structure = None
    frames = []

    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]

        m_cell = cell_header_pattern.match(line)
        if m_cell:
            units = m_cell.group(1).lower() if m_cell.group(1) else "angstrom"
            cell = []
            i += 1

            for _ in range(3):
                if i < n:
                    parts = lines[i].split()
                    if len(parts) >= 3:
                        cell.append([float(parts[0]), float(parts[1]), float(parts[2])])
                    i += 1

            if len(cell) == 3:
                cell = np.array(cell, dtype=float)
                if units == "bohr":
                    cell *= BOHR_TO_ANG
                elif units == "alat":
                    cell *= alat
                last_cell = cell
            continue

        m_pos = atpos_header_pattern.match(line)
        if m_pos:
            units = m_pos.group(1).lower() if m_pos.group(1) else "unknown"
            atoms_block = []
            i += 1

            while i < n:
                m_atom = atom_line_pattern.match(lines[i])
                if m_atom:
                    atoms_block.append((
                        m_atom.group(1),
                        float(m_atom.group(2)),
                        float(m_atom.group(3)),
                        float(m_atom.group(4)),
                    ))
                    i += 1
                else:
                    break

            last_structure = {
                "units": units,
                "atoms": atoms_block,
                "cell": last_cell.copy() if last_cell is not None else None
            }
            continue

        m_energy = energy_pattern.match(line)
        if m_energy and last_structure is not None:
            energy_ev = float(m_energy.group(1)) * RY_TO_EV

            symbols = [a[0] for a in last_structure["atoms"]]
            pos = np.array([[a[1], a[2], a[3]] for a in last_structure["atoms"]], dtype=float)
            units = last_structure["units"]
            cell = last_structure["cell"]

            if units == "bohr":
                pos *= BOHR_TO_ANG
            elif units == "alat":
                pos *= alat
            elif units == "crystal":
                if cell is None:
                    i += 1
                    continue
                pos = pos @ cell

            atoms = Atoms(
                symbols=symbols,
                positions=pos,
                cell=cell,
                pbc=True if cell is not None else False
            )
            atoms.info["energy"] = energy_ev
            atoms.info["energy_per_atom"] = energy_ev / len(atoms)
            frames.append(atoms)

        i += 1

    return frames


def scalar_string(value):
    arr = np.asarray(value)
    return str(arr.item() if arr.shape == () else value)


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


def read_npz_frames(fp):
    data = np.load(fp, allow_pickle=False)

    positions = np.asarray(data["positions"], dtype=float)
    energy_ry = np.asarray(data["energy_ry"], dtype=float)
    symbols = [str(symbol) for symbol in np.asarray(data["symbols"])]
    fallback_cell = fixed_cell_angstrom(data)

    finite_mask = np.isfinite(positions).all(axis=(1, 2)) & np.isfinite(energy_ry)
    valid_indices = np.flatnonzero(finite_mask)

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


def read_dataset_frames(fp):
    if str(fp).lower().endswith(".npz"):
        return read_npz_frames(fp)
    return read_qe_converged_frames(fp)


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
        description="Train the GraphDot GPR model from selected QE output files or NPZ archives."
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
        help="Explicit file paths or glob patterns. Overrides the default phase-based NPZ selection.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Optional label for the current dataset/model run. Default: derived from the selected files.",
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
        help="Build the dataset, write the preview/split files, and stop before fitting the GPR model.",
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
    return slugify_name(f"{phase_label}_{len(file_list)}files")


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
):
    return argparse.Namespace(
        phase=list(phase) if phase is not None else ["bcc", "fcc", "hcp"],
        inputs=list(inputs) if inputs is not None else None,
        run_name=run_name,
        output_root=output_root if output_root is not None else str(PROJECT_ROOT / "ml-results"),
        seed=seed,
        preview_only=preview_only,
    )


def build_dataset_dataframe(all_atoms):
    rows = []
    for sample_index, atoms in enumerate(all_atoms):
        rows.append(
            {
                "sample_index": sample_index,
                "Atoms": atoms,
                "E_per_atom": atoms.info["energy"] / len(atoms),
                "natoms": len(atoms),
                "source_file": atoms.info.get("source_file", "unknown"),
                "source_name": Path(atoms.info.get("source_file", "unknown")).name,
                "frame_index": atoms.info.get("frame_index", -1),
                "phase": infer_phase_name(atoms.info.get("source_file", "")),
            }
        )
    return pd.DataFrame(rows)


def split_dataset(data, seed):
    ntsteps = len(data)
    n_test = ntsteps // 3

    np.random.seed(seed)
    all_idx = np.arange(ntsteps)
    test_sel = np.random.choice(all_idx, n_test, replace=False)
    train_sel = np.setdiff1d(all_idx, test_sel)
    return train_sel, test_sel


def save_dataset_preview(data, train_sel, test_sel, out_dir, base_name):
    preview_df = data.drop(columns=["Atoms", "Graphs"], errors="ignore").copy()
    preview_df["split"] = "train"
    preview_df.loc[test_sel, "split"] = "test"

    preview_csv = out_dir / f"{base_name}_dataset_preview.csv"
    train_csv = out_dir / f"{base_name}_train_split.csv"
    test_csv = out_dir / f"{base_name}_test_split.csv"
    preview_df.to_csv(preview_csv, index=False)
    preview_df.loc[train_sel].to_csv(train_csv, index=False)
    preview_df.loc[test_sel].to_csv(test_csv, index=False)

    fig, axes = plt.subplots(1, 3, figsize=(20, 5.8), dpi=180)

    energy_all = preview_df["E_per_atom"].to_numpy(dtype=float)
    energy_train = preview_df.loc[train_sel, "E_per_atom"].to_numpy(dtype=float)
    energy_test = preview_df.loc[test_sel, "E_per_atom"].to_numpy(dtype=float)
    bins = min(60, max(12, len(preview_df) // 20))

    axes[0].hist(energy_all, bins=bins, alpha=0.35, color="#8a8a8a", label="all")
    axes[0].hist(energy_train, bins=bins, alpha=0.65, color="#1f6f5b", label="train")
    axes[0].hist(energy_test, bins=bins, alpha=0.65, color="#c24b2a", label="test")
    axes[0].set_xlabel("Energy (eV/atom)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Energy Distribution")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False)

    split_colors = {"train": "#1f6f5b", "test": "#c24b2a"}
    colors = [split_colors[split] for split in preview_df["split"]]
    axes[1].scatter(preview_df["sample_index"], preview_df["E_per_atom"], c=colors, s=18, alpha=0.8)
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

    counts = preview_df.groupby(["source_name", "split"]).size().unstack(fill_value=0)
    counts = counts.sort_index()
    x = np.arange(len(counts.index))
    train_counts = counts["train"].to_numpy(dtype=float) if "train" in counts else np.zeros(len(x))
    test_counts = counts["test"].to_numpy(dtype=float) if "test" in counts else np.zeros(len(x))
    axes[2].bar(x, train_counts, color=split_colors["train"], label="train")
    axes[2].bar(x, test_counts, bottom=train_counts, color=split_colors["test"], label="test")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(counts.index, rotation=60, ha="right")
    axes[2].set_ylabel("Frames")
    axes[2].set_title("Frames Per Source File")
    axes[2].grid(True, axis="y", alpha=0.25)
    axes[2].legend(frameon=False)

    fig.suptitle(f"Dataset Preview: {base_name}", fontsize=16)
    fig.tight_layout()
    preview_png = out_dir / f"{base_name}_dataset_preview.png"
    fig.savefig(preview_png, bbox_inches="tight")
    plt.close(fig)

    print("saved dataset preview csv to:", preview_csv)
    print("saved training split csv to:", train_csv)
    print("saved test split csv to:", test_csv)
    print("saved dataset preview figure to:", preview_png)


# def training_graphdot_all(file_list):
#     gpr = GaussianProcessRegressor(
#         kernel=Normalization(MolecularKernel()),
#         alpha=1e-4,
#         optimizer=True,
#         normalize_y=True,
#     )
# #         gpr = gpr.fit(train_data["Graph"], train_data[target], repeat=6,tol=1e-40 , verbose=True)

#     all_atoms = []

#     for file_path in file_list:
#         print(f"\nReading: {file_path}")
#         atoms_list = read_qe_converged_frames(file_path)
#         print(f"  converged frames found: {len(atoms_list)}")
#         all_atoms.extend(atoms_list)

#     if len(all_atoms) == 0:
#         print("No converged frames found in any file.")
#         return None

#     print("\nTotal converged frames read from all files:", len(all_atoms))

#     graphs = [
#         Graph.from_ase(a, adjacency=AtomicAdjacency(shape="compactbell3,2"))
#         for a in all_atoms
#     ]
#     energy_gt = [a.info["energy"] for a in all_atoms]

#     data = pd.DataFrame({
#         "Graphs": graphs,
#         "Pot_Energy": energy_gt
#     })

#     ntsteps = len(all_atoms)
#     target = "Pot_Energy"
    
#     N_test = ntsteps // 3
#     N_train = ntsteps

#     print("N_train =", N_train)
#     print("N_test  =", N_test)

#     np.random.seed(0)

#     all_idx = np.arange(ntsteps)
#     test_sel = np.random.choice(all_idx, N_test, replace=False)
#     train_sel = np.setdiff1d(all_idx, test_sel)

#     train_data = data.iloc[train_sel]
#     test_data = data.iloc[test_sel]

#     start_time = time.time()
#     gpr = gpr.fit(train_data["Graphs"], train_data[target], repeat=6,tol=1e-25 , verbose=True)
#     end_time = time.time()

#     print(
#         "the total time consumption for "
#         + str(len(train_data))
#         + " training steps is "
#         + str((end_time - start_time) / 3600)
#         + " hr."
#     )

#     base_name = "all_sqs_combined"
#     out_dir = "/work/dajuarez4/sqs_files/ml-results/all_sqs_combined"
#     os.makedirs(out_dir, exist_ok=True)

#     fname = f"{base_name}_gpr_DFT_PotEng{len(train_data)}.pkl"
#     gpr.save(out_dir, filename=fname, overwrite=True)

#     print("saved model to:", os.path.join(out_dir, fname))

#     mu = gpr.predict(test_data["Graphs"])
#     return gpr, mu, test_data, target, out_dir, base_name
def training_graphdot_all(file_list, args, base_name, out_dir):
    gpr = GaussianProcessRegressor(
        kernel=Normalization(MolecularKernel()),
        alpha=1e-4,
        optimizer=True,
        normalize_y=True,
    )

    all_atoms = []

    for file_path in file_list:
        print(f"\nReading: {file_path}")
        atoms_list = read_dataset_frames(file_path)
        print(f"  usable frames found: {len(atoms_list)}")
        all_atoms.extend(atoms_list)

    if len(all_atoms) == 0:
        print("No usable frames found in any file.")
        return None

    print("\nTotal converged frames read from all files:", len(all_atoms))

    data = build_dataset_dataframe(all_atoms)

    ntsteps = len(all_atoms)
    target = "E_per_atom"

    graphs = [
        Graph.from_ase(a, adjacency=AtomicAdjacency(shape="compactbell3,2"))
        for a in data["Atoms"]
    ]
    data["Graphs"] = graphs

    train_sel, test_sel = split_dataset(data, args.seed)
    save_dataset_preview(data, train_sel, test_sel, out_dir, base_name)

    if args.preview_only:
        print("Preview-only mode enabled. Skipping model fitting.")
        return None

    nvcc_path = ensure_cuda_compiler_available()
    print("nvcc found at:", nvcc_path)

    N_test = len(test_sel)
    N_train = len(train_sel)

    print("N_train =", N_train)
    print("N_test  =", N_test)

    train_data = data.iloc[train_sel]
    test_data = data.iloc[test_sel]

    start_time = time.time()
    gpr = gpr.fit(train_data["Graphs"], train_data[target], repeat=6, tol=1e-25, verbose=True)
    end_time = time.time()

    print(
        "the total time consumption for "
        + str(len(train_data))
        + " training steps is "
        + str((end_time - start_time) / 3600)
        + " hr."
    )

    fname = f"{base_name}_gpr_DFT_EperAtom_{len(train_data)}.pkl"
    gpr.save(str(out_dir), filename=fname, overwrite=True)

    print("saved model to:", os.path.join(out_dir, fname))

    mu = gpr.predict(test_data["Graphs"])
    return gpr, mu, test_data, target, out_dir, base_name


def execute_ml_gpr(args):
    print('----------------Begin ML---------------------')

    all_sqs_files = resolve_input_files(args)
    base_name = build_run_name(args, all_sqs_files)
    out_dir = build_output_dir(args, base_name)

    print("Files found:")
    for file_path in all_sqs_files:
        print("  ", file_path)

    print("run name:", base_name)
    print("output dir:", out_dir)

    result = training_graphdot_all(all_sqs_files, args, base_name, out_dir)

    if result is not None:
        gpr, mu, test_data, target, out_dir, base_name = result
        plot(mu, test_data, target, out_dir, base_name)
        print(f"saved parity plot to: {os.path.join(out_dir, f'{base_name}_parity_plot.png')}")

    return result


def run_ml_gpr(
    phase=None,
    inputs=None,
    run_name=None,
    output_root=None,
    seed=0,
    preview_only=False,
):
    args = build_runtime_args(
        phase=phase,
        inputs=inputs,
        run_name=run_name,
        output_root=output_root,
        seed=seed,
        preview_only=preview_only,
    )
    return execute_ml_gpr(args)


def main(cli_args=None):
    args = parse_args(cli_args)
    return execute_ml_gpr(args)


if __name__ == "__main__":
    main()
