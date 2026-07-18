#!/usr/bin/env python3
"""Extract the initial and first completed MD configurations from a partial noncollinear QE output."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

from data_compress import (
    BOHR_TO_ANG,
    FLOAT_RE,
    KBAR_TO_GPA,
    parse_atomic_positions_block,
    parse_forces_block,
    parse_initial_info,
    parse_qe_float,
    re_abs_mag,
    re_ekin,
    re_internal_energy,
    re_pressure,
    re_temperature,
    re_total_energy,
)


KB_EV_K = 8.617333262145e-5
RY_TO_EV = 13.605693009
RE_STARTING_TEMPERATURE = re.compile(r"Starting temperature\s*=\s*(" + FLOAT_RE + r")\s*K")
RE_TOTAL_MAG_VECTOR = re.compile(
    r"total magnetization\s*=\s*("
    + FLOAT_RE
    + r")\s+("
    + FLOAT_RE
    + r")\s+("
    + FLOAT_RE
    + r")\s+Bohr mag/cell"
)
RE_TIME_PS = re.compile(r"time\s*=\s*(" + FLOAT_RE + r")\s*pico-seconds")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qe_output", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ideal-reference", type=Path, required=True)
    return parser.parse_args()


def last_match(lines: list[str], start: int, stop: int, pattern: re.Pattern[str]) -> re.Match[str]:
    for index in range(stop - 1, start - 1, -1):
        match = pattern.search(lines[index])
        if match:
            return match
    raise ValueError(f"Could not find {pattern.pattern!r} between lines {start + 1} and {stop}")


def first_match(lines: list[str], start: int, stop: int, pattern: re.Pattern[str]) -> re.Match[str]:
    for index in range(start, stop):
        match = pattern.search(lines[index])
        if match:
            return match
    raise ValueError(f"Could not find {pattern.pattern!r} between lines {start + 1} and {stop}")


def match_float(match: re.Match[str], group: int = 1) -> float:
    return parse_qe_float(match.group(group))


def main() -> int:
    args = parse_args()
    lines = args.qe_output.read_text(errors="replace").splitlines(keepends=True)
    natoms, alat_bohr, initial_cell_alat, initial_positions_alat, _symbols = parse_initial_info(lines)

    force_starts = [index for index, line in enumerate(lines) if "Forces acting on atoms" in line]
    dynamics_starts = [index for index, line in enumerate(lines) if "Entering Dynamics:" in line]
    if len(force_starts) < 2 or len(dynamics_starts) < 2:
        raise ValueError("The output does not contain two complete pre-iteration force evaluations.")

    first_atomic_positions = next(
        index
        for index in range(dynamics_starts[0], force_starts[1])
        if "ATOMIC_POSITIONS" in lines[index]
    )
    first_md_positions, _labels, position_unit, _ = parse_atomic_positions_block(
        lines,
        first_atomic_positions,
        natoms,
    )
    if position_unit != "crystal":
        raise ValueError(f"Expected crystal coordinates, found {position_unit!r}")

    initial_forces, _ = parse_forces_block(lines, force_starts[0], natoms)
    first_md_forces, _ = parse_forces_block(lines, force_starts[1], natoms)
    initial_fractional = initial_positions_alat @ np.linalg.inv(initial_cell_alat)
    positions = np.stack([initial_fractional, first_md_positions])
    forces = np.stack([initial_forces, first_md_forces])

    energy = np.array(
        [
            match_float(last_match(lines, 0, force_starts[0], re_total_energy)),
            match_float(last_match(lines, dynamics_starts[0], force_starts[1], re_total_energy)),
        ],
        dtype=np.float64,
    )
    internal_energy = np.array(
        [
            match_float(last_match(lines, 0, force_starts[0], re_internal_energy)),
            match_float(last_match(lines, dynamics_starts[0], force_starts[1], re_internal_energy)),
        ],
        dtype=np.float64,
    )
    temperature = np.array(
        [
            match_float(first_match(lines, force_starts[0], dynamics_starts[0], RE_STARTING_TEMPERATURE)),
            match_float(last_match(lines, dynamics_starts[0], force_starts[1], re_temperature)),
        ],
        dtype=np.float64,
    )
    pressure_kbar = np.array(
        [
            match_float(first_match(lines, force_starts[0], dynamics_starts[0], re_pressure)),
            match_float(first_match(lines, force_starts[1], dynamics_starts[1], re_pressure)),
        ],
        dtype=np.float64,
    )
    ekin = np.array(
        [
            1.5 * (natoms - 1) * KB_EV_K * temperature[0] / RY_TO_EV,
            match_float(last_match(lines, dynamics_starts[0], force_starts[1], re_ekin)),
        ],
        dtype=np.float64,
    )
    total_mag_vector = np.array(
        [
            [match_float(last_match(lines, 0, force_starts[0], RE_TOTAL_MAG_VECTOR), group) for group in (1, 2, 3)],
            [
                match_float(last_match(lines, dynamics_starts[0], force_starts[1], RE_TOTAL_MAG_VECTOR), group)
                for group in (1, 2, 3)
            ],
        ],
        dtype=np.float64,
    )
    absolute_mag = np.array(
        [
            match_float(last_match(lines, 0, force_starts[0], re_abs_mag)),
            match_float(last_match(lines, dynamics_starts[0], force_starts[1], re_abs_mag)),
        ],
        dtype=np.float64,
    )
    first_time = match_float(first_match(lines, dynamics_starts[0], first_atomic_positions, RE_TIME_PS))
    time_ps = np.array([0.0, first_time], dtype=np.float64)

    input_cell_ang = np.asarray(initial_cell_alat, dtype=np.float64) * alat_bohr * BOHR_TO_ANG
    symbols = np.full(natoms, "Fe", dtype="U2")
    frame_valid = np.ones(2, dtype=bool)
    payload = {
        "input_cell_parameters": input_cell_ang.astype(np.float32),
        "input_cell_unit": np.asarray("angstrom", dtype="U16"),
        "symbols": symbols,
        "species": symbols,
        "source_labels": np.asarray(_symbols, dtype="U8"),
        "initial_positions_alat": np.asarray(initial_positions_alat, dtype=np.float32),
        "initial_cell_alat": np.asarray(initial_cell_alat, dtype=np.float32),
        "positions": positions.astype(np.float32),
        "positions_unit": np.asarray(["crystal", "crystal"], dtype="U16"),
        "cell_parameters": np.repeat(input_cell_ang[None, :, :], 2, axis=0).astype(np.float32),
        "cell_parameters_unit": np.asarray(["angstrom", "angstrom"], dtype="U16"),
        "iteration": np.asarray([0, 1], dtype=np.int32),
        "time_ps": time_ps.astype(np.float32),
        "forces_ry_au": forces.astype(np.float32),
        "force_frame_valid": frame_valid,
        "position_frame_valid": frame_valid,
        "frame_valid": frame_valid,
        "energy_ry": energy,
        "internal_energy_ry": internal_energy,
        "temperature_K": temperature.astype(np.float32),
        "pressure_kbar": pressure_kbar.astype(np.float32),
        "pressure_GPa": (pressure_kbar * KBAR_TO_GPA).astype(np.float32),
        "mag_total_vector_Bohr": total_mag_vector.astype(np.float32),
        "mag_total_Bohr": np.linalg.norm(total_mag_vector, axis=1).astype(np.float32),
        "abs_mag_total_Bohr": absolute_mag.astype(np.float32),
        "ekin_ry": ekin.astype(np.float32),
        "parse_warnings": np.asarray(
            ["The unfinished MD iteration 2 positions were excluded because no completed force block follows them."],
            dtype="U512",
        ),
        "metadata_json": np.asarray(
            json.dumps(
                {
                    "source_file": str(args.qe_output.resolve()),
                    "natoms": natoms,
                    "nsteps": 2,
                    "alat_bohr": alat_bohr,
                    "frame_definition": "initial configuration plus completed MD iteration 1",
                }
            ),
            dtype="U",
        ),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **payload)

    with np.load(args.ideal_reference, allow_pickle=False) as ideal:
        ideal_payload = dict(payload)
        ideal_payload["initial_positions_alat"] = np.asarray(ideal["initial_positions_alat"], dtype=np.float32)
        ideal_payload["initial_cell_alat"] = np.asarray(ideal["initial_cell_alat"], dtype=np.float32)
        ideal_payload["metadata_json"] = np.asarray(
            json.dumps(
                {
                    "source_file": str(args.qe_output.resolve()),
                    "ideal_reference": str(args.ideal_reference.resolve()),
                    "natoms": natoms,
                    "nsteps": 2,
                }
            ),
            dtype="U",
        )
    ideal_output = args.output.with_name(f"{args.output.stem}-ideal.npz")
    np.savez_compressed(ideal_output, **ideal_payload)
    print(f"Wrote {args.output}")
    print(f"Wrote {ideal_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
