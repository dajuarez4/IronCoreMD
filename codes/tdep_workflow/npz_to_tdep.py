#!/usr/bin/env python3
"""Convert QE AIMD NPZ archives into phase-aware TDEP input folders."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

from tdep_common import default_dataset_dir, discover_npz_files
from tdep_phases import PHASE_SPECS, get_phase_spec, write_qpoints_dispersion

RY_TO_EV = 13.605693122994
BOHR_TO_ANG = 0.529177210903
RY_PER_BOHR_TO_EV_PER_ANG = RY_TO_EV / BOHR_TO_ANG


def parse_args(argv: list[str] | None = None, *, default_phase: str = "bcc") -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert QE NPZ files into TDEP folders.")
    parser.add_argument("npz", nargs="*", type=Path, help="NPZ files to convert. Default: all matching NPZ files.")
    parser.add_argument(
        "--phase",
        choices=sorted(PHASE_SPECS),
        default=default_phase,
        help=f"Crystal phase. Default: {default_phase}.",
    )
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=None,
        help="Directory containing the QE NPZ archives. Default: <repo>/dataset/<phase>.",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=None,
        help="Directory where tdep_* folders are written. Default: dataset-dir.",
    )
    parser.add_argument("--supercell", type=int, nargs=3, default=None, metavar=("NX", "NY", "NZ"))
    parser.add_argument("--temperature-K", type=float, default=None, help="Override the temperature written to TDEP files.")
    parser.add_argument("--skip", type=int, default=0, help="Skip this many initial frames.")
    parser.add_argument("--every", type=int, default=1, help="Use every Nth frame after --skip.")
    parser.add_argument("--max-frames", type=int, default=0, help="Maximum frames after filtering. 0 means all.")
    parser.add_argument(
        "--keep-invalid",
        action="store_true",
        help="Keep frames containing NaN or Inf values. Not recommended for TDEP fitting.",
    )
    return parser.parse_args(argv)


def sanitize_stem(stem: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", stem).strip("_")


def scalar_string(value) -> str:
    arr = np.asarray(value)
    return str(arr.item() if arr.shape == () else value)


def positions_alat_to_fractional(positions_alat: np.ndarray, cell_alat: np.ndarray) -> np.ndarray:
    return np.asarray(positions_alat, dtype=float) @ np.linalg.inv(np.asarray(cell_alat, dtype=float))


def write_poscar(path: Path, comment: str, cell: np.ndarray, positions: np.ndarray, symbol: str = "Fe") -> None:
    lines = [comment, "1.0"]
    for row in cell:
        lines.append(f"{row[0]: .16f} {row[1]: .16f} {row[2]: .16f}")
    lines.append(symbol)
    lines.append(str(len(positions)))
    lines.append("Direct")
    for pos in positions:
        wrapped = pos - np.floor(pos)
        lines.append(f"{wrapped[0]: .16f} {wrapped[1]: .16f} {wrapped[2]: .16f}")
    path.write_text("\n".join(lines) + "\n")


def matching_ideal_npz(path: Path) -> Path:
    if not path.name.endswith("-disp.npz"):
        return path
    candidate = path.with_name(path.name.replace("-disp.npz", ".npz"))
    if not candidate.exists():
        raise FileNotFoundError(f"Missing non-displaced reference for {path.name}: {candidate}")
    return candidate


def fixed_cell_angstrom(data) -> np.ndarray:
    cell = np.asarray(data["input_cell_parameters"], dtype=float)
    unit = scalar_string(data["input_cell_unit"]).lower()
    if unit in {"angstrom", "ang"}:
        return cell
    if unit in {"bohr", "a.u.", "au"}:
        return cell * BOHR_TO_ANG
    raise ValueError(f"Unsupported input_cell_unit={unit!r}")


def frame_positions_fractional(data, cell_ang: np.ndarray, selected: np.ndarray) -> np.ndarray:
    positions = np.asarray(data["positions"], dtype=float)[selected]
    units = np.asarray(data["positions_unit"])[selected]
    unique_units = {str(unit).lower() for unit in units}
    if unique_units == {"crystal"}:
        return positions
    if unique_units == {"angstrom"}:
        return positions @ np.linalg.inv(cell_ang)
    raise ValueError(f"Unsupported or mixed positions_unit values: {sorted(unique_units)}")


def infer_dt_fs(time_ps: np.ndarray, selected: np.ndarray) -> float:
    times = np.asarray(time_ps, dtype=float)[selected]
    if len(times) >= 2 and np.isfinite(times[:2]).all():
        return float((times[1] - times[0]) * 1000.0)
    return 1.0


def finite_frame_mask(data) -> np.ndarray:
    positions = np.asarray(data["positions"], dtype=float)
    forces = np.asarray(data["forces_ry_au"], dtype=float)
    energy = np.asarray(data["energy_ry"], dtype=float)
    internal = np.asarray(data["internal_energy_ry"], dtype=float)
    ekin = np.asarray(data["ekin_ry"], dtype=float)
    temp = np.asarray(data["temperature_K"], dtype=float)
    pressure = np.asarray(data["pressure_GPa"], dtype=float)
    return (
        np.isfinite(positions).all(axis=(1, 2))
        & np.isfinite(forces).all(axis=(1, 2))
        & np.isfinite(energy)
        & np.isfinite(internal)
        & np.isfinite(ekin)
        & np.isfinite(temp)
        & np.isfinite(pressure)
    )


def resolve_supercell(
    phase: str,
    cell_ang: np.ndarray,
    natoms: int,
    requested: tuple[int, int, int] | None,
) -> tuple[int, int, int]:
    spec = get_phase_spec(phase)
    inferred = spec.infer_supercell(cell_ang, natoms)
    if requested is not None and tuple(requested) != inferred:
        raise ValueError(f"Requested supercell {tuple(requested)} does not match the inferred {phase} cell shape {inferred}")
    return inferred


def write_tdep_folder(
    npz_path: Path,
    outdir: Path,
    phase: str,
    supercell: tuple[int, int, int] | None,
    temperature_override: float | None,
    skip: int,
    every: int,
    max_frames: int,
    keep_invalid: bool,
) -> None:
    spec = get_phase_spec(phase)
    data = np.load(npz_path, allow_pickle=False)
    ideal_path = matching_ideal_npz(npz_path)
    ideal_data = np.load(ideal_path, allow_pickle=False)

    natoms = int(np.asarray(data["positions"]).shape[1])
    all_indices = np.arange(np.asarray(data["positions"]).shape[0])
    selected = all_indices[skip::every]
    dropped_invalid_frames = 0
    if not keep_invalid:
        selected = selected[finite_frame_mask(data)[selected]]
        dropped_invalid_frames = int(len(all_indices[skip::every]) - len(selected))
    if max_frames > 0:
        selected = selected[:max_frames]
    if len(selected) == 0:
        raise ValueError(f"No valid frames selected for {npz_path}")

    cell_ang = fixed_cell_angstrom(data)
    ideal_cell_ang = fixed_cell_angstrom(ideal_data)
    if not np.allclose(cell_ang, ideal_cell_ang, rtol=1e-6, atol=1e-6):
        raise ValueError(f"Cell mismatch between {npz_path.name} and {ideal_path.name}")

    resolved_supercell = resolve_supercell(phase, cell_ang, natoms, supercell)

    ideal_frac = positions_alat_to_fractional(
        np.asarray(ideal_data["initial_positions_alat"], dtype=float),
        np.asarray(ideal_data["initial_cell_alat"], dtype=float),
    )
    if ideal_frac.shape != (natoms, 3):
        raise ValueError(f"Ideal reference has shape {ideal_frac.shape}, expected {(natoms, 3)}")

    positions_frac = frame_positions_fractional(data, cell_ang, selected)
    forces_ev_ang = np.asarray(data["forces_ry_au"], dtype=float)[selected] * RY_PER_BOHR_TO_EV_PER_ANG
    energies_ev = np.asarray(data["energy_ry"], dtype=float)[selected] * RY_TO_EV
    internal_ev = np.asarray(data["internal_energy_ry"], dtype=float)[selected] * RY_TO_EV
    ekin_ev = np.asarray(data["ekin_ry"], dtype=float)[selected] * RY_TO_EV
    temperatures = np.asarray(data["temperature_K"], dtype=float)[selected]
    pressure_gpa = np.asarray(data["pressure_GPa"], dtype=float)[selected]
    iterations = np.asarray(data["iteration"], dtype=int)[selected]
    time_ps = np.asarray(data["time_ps"], dtype=float)
    times_fs = time_ps[selected] * 1000.0
    if not np.isfinite(times_fs).all():
        times_fs = np.arange(1, len(selected) + 1, dtype=float) * infer_dt_fs(time_ps, selected)

    temperature_meta = float(temperature_override) if temperature_override is not None else float(np.nanmean(temperatures))
    temperatures_for_stat = np.full(len(selected), temperature_meta) if temperature_override is not None else temperatures
    dt_fs = infer_dt_fs(time_ps, selected)

    uc_cell, uc_positions = spec.primitive_cell(cell_ang, natoms)

    outdir.mkdir(parents=True, exist_ok=True)
    write_poscar(outdir / "infile.ucposcar", f"Fe {spec.key} primitive unit cell", uc_cell, uc_positions)
    write_poscar(outdir / "infile.ssposcar", f"Fe {spec.key} ideal supercell", cell_ang, ideal_frac)
    write_qpoints_dispersion(outdir / "infile.qpoints_dispersion", phase)

    with (outdir / "infile.positions").open("w") as handle:
        for frame in positions_frac:
            for pos in frame:
                wrapped = pos - np.floor(pos)
                handle.write(f"{wrapped[0]: .16f} {wrapped[1]: .16f} {wrapped[2]: .16f}\n")

    with (outdir / "infile.forces").open("w") as handle:
        for frame in forces_ev_ang:
            for force in frame:
                handle.write(f"{force[0]: .16f} {force[1]: .16f} {force[2]: .16f}\n")

    with (outdir / "infile.stat").open("w") as handle:
        for iframe in range(len(selected)):
            handle.write(
                f"{times_fs[iframe]: .6f} "
                f"{internal_ev[iframe]: .12f} "
                f"{ekin_ev[iframe]: .12f} "
                f"{energies_ev[iframe]: .12f} "
                f"{temperatures_for_stat[iframe]: .6f} "
                f"{pressure_gpa[iframe]: .8f} "
                "0.0 0.0 0.0 0.0 0.0 0.0\n"
            )

    meta_lines = [
        f"{natoms}",
        f"{len(selected)}",
        f"{dt_fs:.8f}",
        f"{temperature_meta:.8f}",
        f"phase: {spec.key}",
        f"supercell: {' '.join(str(value) for value in resolved_supercell)}",
        f"mean_temperature_K_selected: {float(np.nanmean(temperatures))}",
        f"mean_pressure_GPa_selected: {float(np.nanmean(pressure_gpa))}",
        f"skip: {skip}",
        f"every: {every}",
        f"dropped_invalid_frames: {dropped_invalid_frames}",
    ]
    (outdir / "infile.meta").write_text("\n".join(meta_lines) + "\n")
    (outdir / "source_npz.txt").write_text(
        "\n".join(
            [
                f"phase: {spec.key}",
                f"source_npz: {npz_path.resolve()}",
                f"ideal_reference_npz: {ideal_path.resolve()}",
                f"supercell: {' '.join(str(value) for value in resolved_supercell)}",
                f"frames: {len(selected)}",
                f"skip: {skip}",
                f"every: {every}",
                f"dropped_invalid_frames: {dropped_invalid_frames}",
                f"mean_temperature_K_selected: {float(np.nanmean(temperatures))}",
                f"mean_pressure_GPa_selected: {float(np.nanmean(pressure_gpa))}",
                f"first_iteration: {int(iterations[0])}",
                f"last_iteration: {int(iterations[-1])}",
            ]
        )
        + "\n"
    )


def main(argv: list[str] | None = None, *, default_phase: str = "bcc") -> None:
    args = parse_args(argv, default_phase=default_phase)
    dataset_dir = (args.dataset_dir.resolve() if args.dataset_dir is not None else default_dataset_dir(args.phase).resolve())
    out_root = args.out_root.resolve() if args.out_root is not None else dataset_dir
    requested_supercell = tuple(args.supercell) if args.supercell is not None else None
    npz_files = [path.resolve() if path.is_absolute() else (dataset_dir / path).resolve() for path in args.npz]
    if not npz_files:
        npz_files = discover_npz_files(dataset_dir, None)

    for npz_path in npz_files:
        outdir = out_root / f"tdep_{sanitize_stem(npz_path.stem)}"
        write_tdep_folder(
            npz_path=npz_path,
            outdir=outdir,
            phase=args.phase,
            supercell=requested_supercell,
            temperature_override=args.temperature_K,
            skip=args.skip,
            every=args.every,
            max_frames=args.max_frames,
            keep_invalid=args.keep_invalid,
        )
        print(f"[OK] {npz_path.name} -> {outdir}")


if __name__ == "__main__":
    main()
