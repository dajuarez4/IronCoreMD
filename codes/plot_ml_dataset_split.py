import argparse
import glob
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase import Atoms

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
PHASE_STYLES = {
    "bcc": {"color": "#1f6f5b", "marker": "o"},
    "fcc": {"color": "#c24b2a", "marker": "s"},
    "hcp": {"color": "#2f5597", "marker": "^"},
    "mixed": {"color": "#6b6b6b", "marker": "D"},
}


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
        description="Visualize the pre-model train/test dataset split as energy vs volume."
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
        help="Optional label for the current dataset preview run.",
    )
    parser.add_argument(
        "--output-root",
        default=str(PROJECT_ROOT / "ml-results"),
        help="Root directory where the preview files are written.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed used for the train/test split.",
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
    return slugify_name(f"{phase_label}_{len(file_list)}files_pre_model")


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
):
    return argparse.Namespace(
        phase=list(phase) if phase is not None else ["bcc", "fcc", "hcp"],
        inputs=list(inputs) if inputs is not None else None,
        run_name=run_name,
        output_root=output_root if output_root is not None else str(PROJECT_ROOT / "ml-results"),
        seed=seed,
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


def record_from_atoms(atoms, source_file, frame_index):
    volume = abs(float(atoms.cell.volume)) if atoms.cell.rank == 3 else np.nan
    natoms = len(atoms)
    return {
        "source_file": str(source_file),
        "source_name": Path(source_file).name,
        "phase": infer_phase_name(source_file),
        "frame_index": int(frame_index),
        "natoms": natoms,
        "energy_eV": float(atoms.info["energy"]),
        "energy_eV_atom": float(atoms.info["energy"] / natoms),
        "volume_A3": volume,
        "volume_A3_atom": volume / natoms if np.isfinite(volume) else np.nan,
    }


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

            atoms = Atoms(symbols=symbols, positions=pos, cell=cell, pbc=True if cell is not None else False)
            atoms.info["energy"] = energy_ev
            frames.append(atoms)

        index += 1

    return frames


def read_npz_samples(fp):
    data = np.load(fp, allow_pickle=False)

    positions = np.asarray(data["positions"], dtype=float)
    energy_ry = np.asarray(data["energy_ry"], dtype=float)
    symbols = [str(symbol) for symbol in np.asarray(data["symbols"])]
    fallback_cell = fixed_cell_angstrom(data)

    finite_mask = np.isfinite(positions).all(axis=(1, 2)) & np.isfinite(energy_ry)
    valid_indices = np.flatnonzero(finite_mask)

    records = []
    for idx in valid_indices:
        cell_ang = frame_cell_angstrom(data, idx, fallback_cell)
        pos_ang = frame_positions_angstrom(data, idx, cell_ang)
        atoms = Atoms(symbols=symbols, positions=pos_ang, cell=cell_ang, pbc=True)
        atoms.info["energy"] = float(energy_ry[idx] * RY_TO_EV)
        records.append(record_from_atoms(atoms, fp, idx))
    return records


def read_dataset_samples(fp):
    if str(fp).lower().endswith(".npz"):
        return read_npz_samples(fp)

    frames = read_qe_converged_frames(fp)
    return [record_from_atoms(atoms, fp, idx) for idx, atoms in enumerate(frames)]


def split_dataset(data, seed):
    ntsteps = len(data)
    n_test = ntsteps // 3
    np.random.seed(seed)
    all_idx = np.arange(ntsteps)
    test_sel = np.random.choice(all_idx, n_test, replace=False)
    train_sel = np.setdiff1d(all_idx, test_sel)
    return train_sel, test_sel


def plot_split(data, out_dir, base_name):
    split_colors = {"train": "#1f6f5b", "test": "#c24b2a"}
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=180, sharex=True, sharey=True)

    for axis, split in zip(axes, ("train", "test")):
        subset = data[data["split"] == split]
        for phase in ("bcc", "fcc", "hcp", "mixed"):
            phase_subset = subset[subset["phase"] == phase]
            if phase_subset.empty:
                continue
            style = PHASE_STYLES[phase]
            axis.scatter(
                phase_subset["volume_A3_atom"],
                phase_subset["energy_eV_atom"],
                s=26,
                alpha=0.82,
                color=style["color"],
                marker=style["marker"],
                edgecolors="black",
                linewidths=0.3,
                label=phase,
            )
        axis.set_title(f"{split.capitalize()} Set ({len(subset)} frames)")
        axis.set_xlabel(r"Volume ($\AA^3$/atom)")
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False)

    axes[0].set_ylabel("Energy (eV/atom)")
    fig.suptitle(f"Pre-Model Energy vs Volume Split: {base_name}", fontsize=16)
    fig.tight_layout()

    output_png = out_dir / f"{base_name}_energy_vs_volume_train_test.png"
    fig.savefig(output_png, bbox_inches="tight")
    plt.close(fig)
    print("saved split plot to:", output_png)


def execute_pre_model_plot(args):
    print('----------------Begin Pre-Model Dataset Plot---------------------')

    input_files = resolve_input_files(args)
    base_name = build_run_name(args, input_files)
    out_dir = build_output_dir(args, base_name)

    print("Files found:")
    for file_path in input_files:
        print("  ", file_path)
    print("run name:", base_name)
    print("output dir:", out_dir)

    records = []
    for file_path in input_files:
        print(f"\nReading: {file_path}")
        file_records = read_dataset_samples(file_path)
        print(f"  usable frames found: {len(file_records)}")
        records.extend(file_records)

    if not records:
        raise ValueError("No usable frames found in the selected dataset.")

    data = pd.DataFrame(records)
    data.insert(0, "sample_index", np.arange(len(data), dtype=int))

    train_sel, test_sel = split_dataset(data, args.seed)
    data["split"] = "train"
    data.loc[test_sel, "split"] = "test"

    preview_csv = out_dir / f"{base_name}_energy_volume_dataset.csv"
    train_csv = out_dir / f"{base_name}_energy_volume_train.csv"
    test_csv = out_dir / f"{base_name}_energy_volume_test.csv"
    data.to_csv(preview_csv, index=False)
    data.loc[train_sel].to_csv(train_csv, index=False)
    data.loc[test_sel].to_csv(test_csv, index=False)

    print("\nTotal usable frames read from all files:", len(data))
    print("N_train =", len(train_sel))
    print("N_test  =", len(test_sel))
    print("saved dataset csv to:", preview_csv)
    print("saved training csv to:", train_csv)
    print("saved test csv to:", test_csv)

    plot_split(data, out_dir, base_name)
    return data


def run_pre_model_split(
    phase=None,
    inputs=None,
    run_name=None,
    output_root=None,
    seed=0,
):
    args = build_runtime_args(
        phase=phase,
        inputs=inputs,
        run_name=run_name,
        output_root=output_root,
        seed=seed,
    )
    return execute_pre_model_plot(args)


def main(cli_args=None):
    args = parse_args(cli_args)
    return execute_pre_model_plot(args)


if __name__ == "__main__":
    main()
