"""Export HELD BCC force constants to Phonopy and calculate dispersions."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from phonopy import Phonopy
from phonopy.file_IO import write_FORCE_CONSTANTS
from phonopy.interface.vasp import write_vasp
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from phonopy.structure.atoms import PhonopyAtoms
from scipy.optimize import linear_sum_assignment


HERE = Path(__file__).resolve().parent
IRONCORE_ROOT = HERE.parent
WORKSPACE_ROOT = IRONCORE_ROOT.parent
DATASET_ROOT = IRONCORE_ROOT / "dataset" / "bcc" / "non-mag"
HELD_ROOT = WORKSPACE_ROOT / "HELD"
RESULTS_ROOT = HERE / "results"
PHONOPY_ROOT = HERE / "phonopy"
PATH_LABELS = ("G", "H", "N", "G", "P", "H")
PATH_POINTS = np.asarray(
    [
        [0.0, 0.0, 0.0],
        [0.5, -0.5, 0.5],
        [0.0, 0.0, 0.5],
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
        [0.5, -0.5, 0.5],
    ],
    dtype=float,
)

if str(HELD_ROOT) not in sys.path:
    sys.path.insert(0, str(HELD_ROOT))

from held.model import build_model_from_npz  # noqa: E402


def resolve_source(result) -> Path:
    source = Path(str(result["source_path"].item()))
    if source.is_file():
        return source.resolve()
    fallback = DATASET_ROOT / str(result["source_relative_path"].item())
    if fallback.is_file():
        return fallback.resolve()
    raise FileNotFoundError(f"Cannot locate source trajectory: {source} or {fallback}")


def select_coefficients(result, md_step: int | None) -> tuple[np.ndarray, str, int | None]:
    if md_step is None:
        return np.asarray(result["held_mean_coefficients"], dtype=float), "mean", None
    step_ids = np.asarray(result["step_ids"], dtype=int)
    matches = np.flatnonzero(step_ids == int(md_step))
    if len(matches) != 1:
        raise KeyError(f"Expected one stored MD step {md_step}; found {len(matches)}")
    coefficients = np.asarray(result["held_coefficients_per_frame"][matches[0]], dtype=float)
    return coefficients, f"md_step_{int(md_step):06d}", int(md_step)


def full_dataset_force_constants(model, coefficients: np.ndarray) -> np.ndarray:
    """Expand primitive HELD tensors to the source-supercell atom ordering."""
    offsite, _onsite = model.primitive_force_constants_from_coefficients(coefficients)
    force_constants = np.zeros((model.n_supercell, model.n_supercell, 3, 3), dtype=float)
    for (atom_i, atom_j), key in model.full_pair_key_map.items():
        force_constants[atom_i, atom_j] = offsite[model.offsite_index[key]]
    for atom_i in range(model.n_supercell):
        force_constants[atom_i, atom_i] = -np.sum(force_constants[atom_i], axis=0)
    return force_constants


def make_phonopy(model) -> Phonopy:
    unitcell = PhonopyAtoms(
        symbols=["Fe"],
        cell=model.uc_cell,
        scaled_positions=model.uc_frac,
    )
    return Phonopy(
        unitcell,
        supercell_matrix=model.supercell_transform,
        primitive_matrix=np.eye(3),
        symprec=model.symprec,
    )


def phonopy_to_dataset_map(phonon: Phonopy, model) -> np.ndarray:
    """Map Phonopy supercell order to the AIMD/HELD ideal-supercell order."""
    phonopy_cart = np.asarray(phonon.supercell.scaled_positions) @ np.asarray(phonon.supercell.cell)
    phonopy_frac_in_dataset = phonopy_cart @ np.linalg.inv(model.ss_cell)
    delta = phonopy_frac_in_dataset[:, None, :] - model.ss_frac[None, :, :]
    delta -= np.rint(delta)
    distance = np.linalg.norm(delta @ model.ss_cell, axis=2)
    rows, cols = linear_sum_assignment(distance)
    if not np.array_equal(rows, np.arange(model.n_supercell)):
        raise RuntimeError("Unexpected incomplete Phonopy-to-dataset atom assignment")
    maximum = float(distance[rows, cols].max())
    if maximum > 1.0e-5:
        raise ValueError(f"Phonopy/dataset atom mapping error is {maximum:.3e} angstrom")
    return cols.astype(int)


def band_conf_text(transform: np.ndarray, points_per_segment: int) -> str:
    dim = " ".join(str(int(value)) for value in np.asarray(transform).reshape(-1))
    band = "  ".join(" ".join(f"{value:g}" for value in point) for point in PATH_POINTS)
    return (
        "ATOM_NAME = Fe\n"
        f"DIM = {dim}\n"
        "PRIMITIVE_AXES = 1 0 0  0 1 0  0 0 1\n"
        f"BAND = {band}\n"
        f"BAND_LABELS = {' '.join(PATH_LABELS)}\n"
        f"BAND_POINTS = {int(points_per_segment)}\n"
        "FORCE_CONSTANTS = READ\n"
    )


def _plot_dispersion(distances, frequencies, output: Path, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(8.5, 5.8), constrained_layout=True)
    boundaries = [float(distances[0][0])]
    for segment_distance, segment_frequency in zip(distances, frequencies):
        for band in range(segment_frequency.shape[1]):
            axis.plot(segment_distance, segment_frequency[:, band], color="#5e2ca5", linewidth=1.4)
        boundaries.append(float(segment_distance[-1]))
    for boundary in boundaries:
        axis.axvline(boundary, color="0.65", linewidth=0.7)
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xticks(boundaries)
    axis.set_xticklabels(["Γ" if label == "G" else label for label in PATH_LABELS])
    axis.set_xlim(boundaries[0], boundaries[-1])
    axis.set_xlabel("Wave-vector path")
    axis.set_ylabel("Frequency (THz)")
    axis.set_title(title)
    axis.grid(axis="y", alpha=0.2)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def export_and_run(
    result_npz: str | Path,
    md_step: int | None = None,
    output_root: str | Path = PHONOPY_ROOT,
    points_per_segment: int = 101,
    validate_against_held: bool = True,
) -> dict[str, object]:
    """Create a Phonopy-ready folder and calculate one HELD dispersion."""
    result_npz = Path(result_npz).resolve()
    with np.load(result_npz, allow_pickle=False) as result:
        if str(result["schema_version"].item()) != "held-bcc-5nn-bvk14-v2":
            raise ValueError(f"Unsupported HELD result schema in {result_npz}")
        case_id = str(result["case_id"].item())
        source = resolve_source(result)
        coefficients, selection, selected_step = select_coefficients(result, md_step)

    model, _dataset, metadata = build_model_from_npz("bcc", source, num_shells=5)
    dataset_fc = full_dataset_force_constants(model, coefficients)
    phonon = make_phonopy(model)
    mapping = phonopy_to_dataset_map(phonon, model)
    phonopy_fc = dataset_fc[np.ix_(mapping, mapping)]
    phonon.force_constants = phonopy_fc

    output = Path(output_root) / case_id / selection
    output.mkdir(parents=True, exist_ok=True)
    write_vasp(str(output / "POSCAR"), phonon.unitcell, direct=True)
    write_FORCE_CONSTANTS(phonopy_fc, filename=str(output / "FORCE_CONSTANTS"))
    (output / "band.conf").write_text(
        band_conf_text(model.supercell_transform, points_per_segment), encoding="utf-8"
    )

    path_endpoints = [[PATH_POINTS[index], PATH_POINTS[index + 1]] for index in range(len(PATH_POINTS) - 1)]
    paths, connections = get_band_qpoints_and_path_connections(
        path_endpoints, npoints=points_per_segment
    )
    phonon.run_band_structure(
        paths,
        path_connections=connections,
        labels=list(PATH_LABELS),
        with_eigenvectors=False,
    )
    phonon.write_yaml_band_structure(filename=str(output / "band.yaml"))
    bands = phonon.get_band_structure_dict()
    qpoints = [np.asarray(values, dtype=float) for values in bands["qpoints"]]
    distances = [np.asarray(values, dtype=float) for values in bands["distances"]]
    frequencies = [np.asarray(values, dtype=float) for values in bands["frequencies"]]

    held_max_abs_delta = float("nan")
    if validate_against_held:
        held_segments = [model.dispersion_thz_from_reduced_path(coefficients, values) for values in qpoints]
        differences = [
            np.max(np.abs(np.sort(phonopy_values, axis=1) - np.sort(held_values, axis=1)))
            for phonopy_values, held_values in zip(frequencies, held_segments)
        ]
        held_max_abs_delta = float(max(differences))
        if held_max_abs_delta > 1.0e-5:
            raise ValueError(
                f"Phonopy and HELD dispersions disagree by {held_max_abs_delta:.6e} THz"
            )

    np.savez_compressed(
        output / "phonon_dispersion.npz",
        case_id=np.asarray(case_id),
        selection=np.asarray(selection),
        md_step=np.asarray(-1 if selected_step is None else selected_step, dtype=int),
        path_labels=np.asarray(PATH_LABELS, dtype="U2"),
        qpoints=np.stack(qpoints),
        distances=np.stack(distances),
        frequencies_THz=np.stack(frequencies),
        held_max_abs_delta_THz=np.asarray(held_max_abs_delta),
    )
    rows = []
    for segment, (q_segment, d_segment, f_segment) in enumerate(zip(qpoints, distances, frequencies)):
        for point, (qpoint, distance, frequency) in enumerate(zip(q_segment, d_segment, f_segment)):
            rows.append([segment, point, distance, *qpoint, *frequency])
    header = "segment point distance q1 q2 q3 frequency_1_THz frequency_2_THz frequency_3_THz"
    np.savetxt(output / "phonon_dispersion.dat", np.asarray(rows), header=header, fmt="%.10f")
    _plot_dispersion(
        distances,
        frequencies,
        output / "phonon_dispersion.png",
        f"BCC Fe HELD→Phonopy: {case_id} ({selection})",
    )
    summary = {
        "case_id": case_id,
        "source_result_npz": str(result_npz),
        "source_trajectory_npz": str(source),
        "selection": selection,
        "md_step": selected_step,
        "phase": metadata.phase,
        "symbol": metadata.symbol,
        "natoms_supercell": model.n_supercell,
        "supercell_matrix": model.supercell_transform.tolist(),
        "phonopy_to_dataset_atom_map": mapping.tolist(),
        "points_per_segment": points_per_segment,
        "path_labels": list(PATH_LABELS),
        "frequency_units": "THz",
        "force_constant_units": "eV/angstrom^2",
        "held_max_abs_delta_THz": held_max_abs_delta,
    }
    (output / "metadata.json").write_text(json.dumps(summary, indent=2) + "\n")
    return {**summary, "output_directory": str(output)}
