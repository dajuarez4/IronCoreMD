"""Batch-fit BCC trajectories with HELD through the fifth neighbor shell.

The complete HELD basis and force-constant tensors are retained.  The 14
traditional monoatomic-BCC Born-von Karman elements (onsite plus 1NN--5NN)
are extracted after rotating every symmetry-equivalent bond tensor to its
canonical cubic orientation.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


HERE = Path(__file__).resolve().parent
IRONCORE_ROOT = HERE.parent
WORKSPACE_ROOT = IRONCORE_ROOT.parent
DATASET_ROOT = IRONCORE_ROOT / "dataset" / "bcc" / "non-mag"
HELD_ROOT = WORKSPACE_ROOT / "HELD"
RESULTS_ROOT = HERE / "results"
BVK_LABELS = np.asarray(
    [
        "alpha_0",
        "alpha_1", "beta_1",
        "alpha_2", "beta_2",
        "alpha_3", "beta_3", "gamma_3",
        "alpha_4", "beta_4", "gamma_4", "delta_4",
        "alpha_5", "beta_5",
    ],
    dtype="U8",
)
REQUIRED_KEYS = {
    "input_cell_parameters",
    "input_cell_unit",
    "positions",
    "positions_unit",
    "initial_positions_alat",
    "initial_cell_alat",
    "forces_ry_au",
    "symbols",
}
EXCLUDED_NAMES = {
    "simulation-ideal.npz",
    "held_heatmap_steps_all_complete.npz",
    "round19_held_phonon_frames.npz",
}

if str(HELD_ROOT) not in sys.path:
    sys.path.insert(0, str(HELD_ROOT))

from held.model import build_model_from_npz  # noqa: E402


@dataclass(frozen=True)
class FitConfig:
    skip: int = 0
    every: int = 1
    max_frames: int = 0
    aggregate: str = "mean"
    num_shells: int = 5
    overwrite: bool = False


def is_held_trajectory(path: Path) -> tuple[bool, str]:
    """Return whether an NPZ has the HELD trajectory schema and usable arrays."""
    if path.name in EXCLUDED_NAMES or "ForceConstants" in path.parts:
        return False, "excluded generated or ideal-reference file"
    try:
        with np.load(path, allow_pickle=False) as data:
            missing = sorted(REQUIRED_KEYS.difference(data.files))
            if missing:
                return False, "missing keys: " + ", ".join(missing)
            positions = np.asarray(data["positions"])
            forces = np.asarray(data["forces_ry_au"])
            symbols = np.asarray(data["symbols"])
            if positions.ndim != 3 or positions.shape[-1] != 3:
                return False, f"invalid positions shape {positions.shape}"
            if forces.shape != positions.shape:
                return False, f"force shape {forces.shape} != positions shape {positions.shape}"
            if positions.shape[0] == 0:
                return False, "zero trajectory frames"
            if symbols.shape != (positions.shape[1],):
                return False, f"symbols shape {symbols.shape} != ({positions.shape[1]},)"
    except Exception as exc:
        return False, f"cannot read NPZ: {type(exc).__name__}: {exc}"
    return True, "ok"


def discover_bcc_npz(dataset_root: Path = DATASET_ROOT) -> tuple[list[Path], list[dict[str, str]]]:
    """Discover valid input trajectories and record every skipped NPZ."""
    valid: list[Path] = []
    skipped: list[dict[str, str]] = []
    for path in sorted(Path(dataset_root).rglob("*.npz")):
        accepted, reason = is_held_trajectory(path)
        if accepted:
            valid.append(path.resolve())
        else:
            skipped.append({"path": str(path.resolve()), "reason": reason})
    return valid, skipped


def case_id(path: Path, dataset_root: Path = DATASET_ROOT) -> str:
    """Make a readable, collision-safe identifier from a dataset-relative path."""
    try:
        relative = path.resolve().relative_to(Path(dataset_root).resolve())
    except ValueError:
        relative = path.resolve()
    stem = re.sub(r"[^A-Za-z0-9]+", "_", str(relative.with_suffix(""))).strip("_")
    digest = hashlib.sha1(str(relative).encode("utf-8")).hexdigest()[:8]
    return f"{stem}__{digest}"


def _shell_indices(model) -> list[np.ndarray]:
    distances = np.linalg.norm(model.offsite_separations_cart, axis=1)
    tolerance = 1.0e-5 * min(np.linalg.norm(vector) for vector in model.uc_cell)
    indices = [
        np.flatnonzero(np.isclose(distances, shell_distance, atol=tolerance, rtol=0.0))
        for shell_distance in model.selected_shell_distances
    ]
    if len(indices) != 5 or any(len(index) == 0 for index in indices):
        raise ValueError("HELD did not identify five populated BCC neighbor shells")
    return indices


def _canonical_rotation(vector: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Return a signed permutation Q such that Q @ vector points to target."""
    vector = np.asarray(vector, dtype=float)
    target = np.asarray(target, dtype=float)
    rotation = np.zeros((3, 3), dtype=float)
    unused = set(range(3))
    for canonical_axis, target_value in enumerate(target):
        candidates = [
            axis for axis in unused
            if np.isclose(abs(vector[axis]), abs(target_value), rtol=1.0e-5, atol=1.0e-8)
        ]
        if not candidates:
            raise ValueError(f"Cannot map BCC bond {vector} to canonical target {target}")
        source_axis = min(candidates)
        unused.remove(source_axis)
        if np.isclose(target_value, 0.0):
            sign = 1.0
        else:
            sign = np.sign(target_value / vector[source_axis])
        rotation[canonical_axis, source_axis] = sign
    if not np.allclose(rotation @ vector, target, rtol=1.0e-5, atol=1.0e-8):
        raise ValueError(f"Signed permutation failed for BCC bond {vector}")
    return rotation


def _canonical_shell_tensors(model, offsite: np.ndarray) -> list[np.ndarray]:
    """Average every shell after mapping bonds to standard BCC directions."""
    lattice_a = float(model.selected_shell_distances[1])
    canonical = 0.5 * lattice_a * np.asarray(
        [[1, 1, 1], [2, 0, 0], [0, 2, 2], [3, 1, 1], [2, 2, 2]],
        dtype=float,
    )
    shell_tensors: list[np.ndarray] = []
    for shell_number, indices in enumerate(_shell_indices(model)):
        transformed = []
        for pair_index in indices:
            rotation = _canonical_rotation(
                model.offsite_separations_cart[pair_index], canonical[shell_number]
            )
            transformed.append(
                np.einsum(
                    "ij,fjk,lk->fil",
                    rotation,
                    offsite[:, pair_index],
                    rotation,
                    optimize=True,
                )
            )
        shell_tensors.append(np.mean(np.stack(transformed, axis=1), axis=1))
    return shell_tensors


def bvk_14_constants(model, coefficient_rows: np.ndarray) -> np.ndarray:
    """Extract [alpha_0, alpha_1, beta_1, ..., alpha_5, beta_5]."""
    coefficients = np.asarray(coefficient_rows, dtype=float)
    was_vector = coefficients.ndim == 1
    if was_vector:
        coefficients = coefficients[None, :]
    offsite = np.tensordot(coefficients, model.basis_offsite_mats, axes=(1, 0))
    onsite = np.tensordot(coefficients, model.basis_onsite_mats, axes=(1, 0))
    shell1, shell2, shell3, shell4, shell5 = _canonical_shell_tensors(model, offsite)

    def diagonal_mean(matrix: np.ndarray) -> np.ndarray:
        return np.diagonal(matrix, axis1=1, axis2=2).mean(axis=1)

    def symmetric(matrix: np.ndarray, axis_a: int, axis_b: int) -> np.ndarray:
        return 0.5 * (matrix[:, axis_a, axis_b] + matrix[:, axis_b, axis_a])

    values = np.column_stack(
        [
            np.diagonal(onsite, axis1=2, axis2=3).mean(axis=(1, 2)),
            diagonal_mean(shell1),
            np.mean(
                np.column_stack([symmetric(shell1, 0, 1), symmetric(shell1, 0, 2), symmetric(shell1, 1, 2)]),
                axis=1,
            ),
            shell2[:, 0, 0],
            0.5 * (shell2[:, 1, 1] + shell2[:, 2, 2]),
            shell3[:, 0, 0],
            0.5 * (shell3[:, 1, 1] + shell3[:, 2, 2]),
            symmetric(shell3, 1, 2),
            shell4[:, 0, 0],
            0.5 * (shell4[:, 1, 1] + shell4[:, 2, 2]),
            0.5 * (symmetric(shell4, 0, 1) + symmetric(shell4, 0, 2)),
            symmetric(shell4, 1, 2),
            diagonal_mean(shell5),
            np.mean(
                np.column_stack([symmetric(shell5, 0, 1), symmetric(shell5, 0, 2), symmetric(shell5, 1, 2)]),
                axis=1,
            ),
        ]
    )
    return values[0] if was_vector else values


def _output_path(source: Path, results_root: Path, dataset_root: Path) -> Path:
    return Path(results_root) / f"{case_id(source, dataset_root)}_held_5nn.npz"


def prepare_iron_input(source: Path, results_root: Path) -> tuple[Path, list[str]]:
    """Validate that the selected non-magnetic trajectory is elemental Fe."""
    source = Path(source).resolve()
    notes: list[str] = []
    with np.load(source, allow_pickle=False) as data:
        symbols = np.asarray(data["symbols"]).astype(str)
        unique = sorted(set(symbols.tolist()))
        if unique != ["Fe"]:
            raise ValueError(
                f"Expected elemental non-magnetic Fe; found symbols {unique}"
            )
    return source, notes


def fit_trajectory(
    source: Path,
    config: FitConfig = FitConfig(),
    results_root: Path = RESULTS_ROOT,
    dataset_root: Path = DATASET_ROOT,
    verbose: bool = True,
) -> dict[str, object]:
    """Fit one trajectory and write a self-describing compressed NPZ result."""
    source = Path(source).resolve()
    output = _output_path(source, results_root, dataset_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    cache_is_current = False
    if output.exists():
        with np.load(output, allow_pickle=False) as saved:
            cache_is_current = (
                "schema_version" in saved.files
                and str(saved["schema_version"].item()) == "held-bcc-5nn-bvk14-v2"
            )
    if output.exists() and not config.overwrite and cache_is_current:
        with np.load(output, allow_pickle=False) as saved:
            return {
                "case_id": str(saved["case_id"].item()),
                "source_path": str(saved["source_path"].item()),
                "output_path": str(output),
                "status": "cached",
                "n_frames": int(saved["n_frames"].item()),
                "natoms": int(saved["natoms"].item()),
                "fc_mean": np.asarray(saved["fc_mean"], dtype=float),
                "shell_distances_ang": np.asarray(saved["shell_distances_ang"], dtype=float),
                "elapsed_s": 0.0,
            }

    started = time.perf_counter()
    fit_source, preparation_notes = prepare_iron_input(source, results_root)
    model, dataset, metadata = build_model_from_npz(
        phase="bcc", npz_path=fit_source, num_shells=config.num_shells
    )
    frame_indices = dataset.select_frame_indices(
        skip=config.skip, every=config.every, max_frames=config.max_frames
    )
    positions = dataset.positions_frac[frame_indices]
    forces = dataset.forces_ev_ang[frame_indices]
    steps = dataset.step_ids[frame_indices]
    result = model.fit_series(
        positions, forces, steps, aggregate=config.aggregate, verbose=verbose
    )

    fc_per_md_step = bvk_14_constants(model, result.step_values)
    fc_mean = bvk_14_constants(model, result.mean_values)
    offsite_mean, onsite_mean = model.primitive_force_constants_from_coefficients(
        result.mean_values
    )
    offsite_per_frame = np.tensordot(
        result.step_values, model.basis_offsite_mats, axes=(1, 0)
    )
    onsite_per_frame = np.tensordot(
        result.step_values, model.basis_onsite_mats, axes=(1, 0)
    )
    shell_indices = _shell_indices(model)
    shell_id_per_pair = np.zeros(len(model.offsite_keys), dtype=int)
    for shell_number, indices in enumerate(shell_indices, 1):
        shell_id_per_pair[indices] = shell_number

    relative_source = str(source.relative_to(dataset_root.resolve()))
    metadata_json = json.dumps(
        {
            "definition": (
                "14 monoatomic-BCC Born-von Karman elements through 5NN; "
                "alpha_0 is onsite and is derived by HELD's acoustic sum rule"
            ),
            "independent_fc_elements": 13,
            "stored_fc_elements": 14,
            "fc_units": "eV/angstrom^2",
            "tensor_units": "eV/angstrom^2",
            "phase": metadata.phase,
            "symbol": metadata.symbol,
            "mass_amu": metadata.mass_amu,
            "fit_config": asdict(config),
            "input_preparation": preparation_notes,
        },
        sort_keys=True,
    )
    np.savez_compressed(
        output,
        schema_version=np.asarray("held-bcc-5nn-bvk14-v2"),
        case_id=np.asarray(case_id(source, dataset_root)),
        source_path=np.asarray(str(source)),
        source_relative_path=np.asarray(relative_source),
        phase=np.asarray("bcc"),
        symbol=np.asarray(metadata.symbol),
        natoms=np.asarray(dataset.natoms, dtype=int),
        n_frames=np.asarray(len(frame_indices), dtype=int),
        frame_indices=frame_indices,
        step_ids=np.asarray(result.step_ids, dtype=int),
        fc_labels=BVK_LABELS,
        fc_mean=fc_mean,
        fc_per_md_step=fc_per_md_step,
        fc_units=np.asarray("eV/angstrom^2"),
        shell_distances_ang=np.asarray(metadata.selected_shell_distances, dtype=float),
        held_basis_labels=np.asarray(result.labels, dtype="U16"),
        held_mean_coefficients=result.mean_values,
        held_fitted_sigmas=result.fitted_sigmas,
        held_coefficients_per_frame=result.step_values,
        offsite_pair_keys=np.asarray(
            [[key[0], key[1], *key[2]] for key in model.offsite_keys], dtype=int
        ),
        offsite_shell_id=shell_id_per_pair,
        offsite_separations_ang=model.offsite_separations_cart,
        offsite_fc_mean=offsite_mean,
        onsite_fc_mean=onsite_mean,
        offsite_fc_per_frame=offsite_per_frame,
        onsite_fc_per_frame=onsite_per_frame,
        tensor_units=np.asarray("eV/angstrom^2"),
        metadata_json=np.asarray(metadata_json),
    )
    elapsed = time.perf_counter() - started
    return {
        "case_id": case_id(source, dataset_root),
        "source_path": str(source),
        "output_path": str(output),
        "status": "completed",
        "n_frames": len(frame_indices),
        "natoms": dataset.natoms,
        "fc_mean": fc_mean,
        "shell_distances_ang": np.asarray(metadata.selected_shell_distances),
        "elapsed_s": elapsed,
    }


def write_master_index(
    records: Iterable[dict[str, object]],
    skipped: Iterable[dict[str, str]],
    output: Path = RESULTS_ROOT / "bcc_held_5nn_index.npz",
) -> Path:
    """Write a pickle-free index spanning successful, cached, and failed cases."""
    rows = list(records)
    skipped_rows = list(skipped)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fc_values = np.full((len(rows), 14), np.nan, dtype=float)
    distances = np.full((len(rows), 5), np.nan, dtype=float)
    for index, row in enumerate(rows):
        if "fc_mean" in row:
            fc_values[index] = np.asarray(row["fc_mean"], dtype=float)
        if "shell_distances_ang" in row:
            distances[index] = np.asarray(row["shell_distances_ang"], dtype=float)
    np.savez_compressed(
        output,
        schema_version=np.asarray("held-bcc-5nn-bvk14-index-v2"),
        fc_labels=BVK_LABELS,
        case_id=np.asarray([str(row.get("case_id", "")) for row in rows], dtype="U512"),
        source_path=np.asarray([str(row.get("source_path", "")) for row in rows], dtype="U2048"),
        output_path=np.asarray([str(row.get("output_path", "")) for row in rows], dtype="U2048"),
        status=np.asarray([str(row.get("status", "")) for row in rows], dtype="U32"),
        error=np.asarray([str(row.get("error", "")) for row in rows], dtype="U4096"),
        n_frames=np.asarray([int(row.get("n_frames", 0)) for row in rows], dtype=int),
        natoms=np.asarray([int(row.get("natoms", 0)) for row in rows], dtype=int),
        elapsed_s=np.asarray([float(row.get("elapsed_s", 0.0)) for row in rows]),
        fc_mean=fc_values,
        fc_units=np.asarray("eV/angstrom^2"),
        shell_distances_ang=distances,
        skipped_path=np.asarray([row["path"] for row in skipped_rows], dtype="U2048"),
        skipped_reason=np.asarray([row["reason"] for row in skipped_rows], dtype="U4096"),
    )
    return output


def run_dataset(
    paths: Iterable[Path] | None = None,
    config: FitConfig = FitConfig(),
    dataset_root: Path = DATASET_ROOT,
    results_root: Path = RESULTS_ROOT,
    continue_on_error: bool = True,
    verbose: bool = True,
) -> tuple[list[dict[str, object]], list[dict[str, str]], Path]:
    """Fit all discovered BCC trajectories and always write a master index."""
    discovered, skipped = discover_bcc_npz(dataset_root)
    selected = discovered if paths is None else [Path(path).resolve() for path in paths]
    records: list[dict[str, object]] = []
    for number, path in enumerate(selected, 1):
        if verbose:
            print(f"[dataset {number}/{len(selected)}] {path.relative_to(dataset_root)}")
        try:
            records.append(
                fit_trajectory(
                    path,
                    config=config,
                    results_root=results_root,
                    dataset_root=dataset_root,
                    verbose=verbose,
                )
            )
        except Exception as exc:
            record = {
                "case_id": case_id(path, dataset_root),
                "source_path": str(path),
                "output_path": str(_output_path(path, results_root, dataset_root)),
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
                "n_frames": 0,
                "natoms": 0,
                "elapsed_s": 0.0,
            }
            records.append(record)
            print(f"  FAILED: {record['error']}")
            if not continue_on_error:
                raise
    index = write_master_index(
        records,
        skipped,
        output=Path(results_root) / "bcc_held_5nn_index.npz",
    )
    return records, skipped, index


def summarize_index(path: Path = RESULTS_ROOT / "bcc_held_5nn_index.npz") -> None:
    with np.load(path, allow_pickle=False) as data:
        statuses, counts = np.unique(data["status"], return_counts=True)
        print("Results:", dict(zip(statuses.tolist(), counts.tolist())))
        print("Skipped NPZ files:", len(data["skipped_path"]))
        print("FC labels:", data["fc_labels"].tolist())
        print("FC units:", data["fc_units"].item())
        print("index:", Path(path).resolve())
