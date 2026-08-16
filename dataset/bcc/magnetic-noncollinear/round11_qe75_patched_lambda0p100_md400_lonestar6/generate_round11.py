#!/usr/bin/env python3
"""Generate the patched-QE 7.5 noncollinear 400-step Round-11 run."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ROUND10 = ROOT.parent / "round10_qe75_patched_lambda0p100_md50_lonestar6"
PSEUDO = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
LABELS = [f"Fe{i:02d}" for i in range(1, 9)]
ANGLES = (
    (72.8157590495, 72.7795960677),
    (73.5303944933, 92.3546343919),
    (33.6579186653, 331.8844103957),
    (95.5478199306, 0.9655109597),
    (107.1842409505, 252.7795960677),
    (106.4696055067, 272.3546343919),
    (146.3420813347, 151.8844103957),
    (84.4521800694, 180.9655109597),
)
POSITIONS = """ATOMIC_POSITIONS crystal
Fe01 0.012027047575 0.073702611029 0.968028113246
Fe02 0.284393727779 0.288436353207 0.473334848881
Fe03 0.984382264316 0.580798864365 0.031577926129
Fe04 0.243689253926 0.798996508121 0.501007080078
Fe05 0.466723054647 0.925974182785 0.995024635457
Fe06 0.772856652737 0.200137659907 0.519648075104
Fe07 0.513922214508 0.431229859591 0.018644046038
Fe08 0.722005784512 0.700724005699 0.492735266685""".splitlines()
VELOCITIES = """ATOMIC_VELOCITIES { a.u }
Fe01 -5.8056229450392154e-04  6.5703246584198729e-05  1.1226770978898905e-03
Fe02 -7.8499856503635166e-04 -8.8175025581256074e-04  4.8300684775749253e-04
Fe03 -6.1115577772391939e-04  9.1417844027661028e-04 -8.2963795741283018e-05
Fe04 -1.1900374259849299e-04 -2.6646984989837588e-04 -7.4874180322197312e-04
Fe05 -2.8546093095084936e-04 -4.6743783269815096e-04  5.1918724377899097e-04
Fe06  3.2015656969003526e-04  4.4415703705978881e-04 -4.0938205212582174e-04
Fe07  1.9123073094099092e-03  3.8614250271665397e-04 -6.2671987999268156e-04
Fe08  1.4871743171359028e-04 -1.9452328822816396e-04 -2.5706365834461437e-04""".splitlines()
CELL = """CELL_PARAMETERS angstrom
4.780000000000 0.000000000000 0.000000000000
0.000000000000 4.780000000000 0.000000000000
0.000000000000 0.000000000000 2.390000000000""".splitlines()
CASES = (
    ("a2p39_m2_lambda0p100_md400", (0.005, 0.020, 0.100)),
)
STAGES = (
    ("stage1_4e-3", 4.0e-3, 500),
    ("stage2_1e-3", 1.0e-3, 500),
    ("stage3_3e-4", 3.0e-4, 600),
)


def system(lambda_value: float) -> list[str]:
    lines = [
        "&SYSTEM", "   ibrav=0", "   nat=8", "   ntyp=8",
        "   ecutwfc=71.0", "   ecutrho=496.0",
        "   occupations='smearing'", "   smearing='fd'",
        "   degauss=0.025334", "   nosym=.true.",
        "   noncolin=.true.", "   lspinorb=.false.",
        "   constrained_magnetization='atomic'",
        f"   lambda={lambda_value:.8f}", "   report=1",
    ]
    for index, (theta, phi) in enumerate(ANGLES, 1):
        lines.extend((
            f"   starting_magnetization({index})=0.12500000",
            f"   angle1({index})={theta:.10f}",
            f"   angle2({index})={phi:.10f}",
        ))
    return lines + ["/"]


def electrons(read_state: bool, threshold: float, maxstep: int) -> list[str]:
    return [
        "&ELECTRONS",
        f"   startingpot='{'file' if read_state else 'atomic'}'",
        f"   startingwfc='{'file' if read_state else 'atomic+random'}'",
        f"   conv_thr={threshold:.1e}", f"   electron_maxstep={maxstep}",
        "   scf_must_converge=.true.", "   mixing_mode='local-TF'",
        "   mixing_beta=0.010d0", "   mixing_ndim=20",
        "   diagonalization='cg'", "   diago_thr_init=1.0d-3",
        "   diago_cg_maxiter=200", "/",
    ]


def cards(include_velocities: bool) -> list[str]:
    lines = ["ATOMIC_SPECIES"]
    lines.extend(f"{label} 55.845 {PSEUDO}" for label in LABELS)
    lines.extend(POSITIONS)
    if include_velocities:
        lines.extend(VELOCITIES)
    return lines + CELL + ["K_POINTS automatic", "1 1 1 0 0 0", ""]


def scf_input(prefix: str, path: tuple[float, ...], stage_index: int) -> str:
    _, threshold, maxstep = STAGES[stage_index]
    lines = [
        "&CONTROL", "   calculation='scf'", "   restart_mode='from_scratch'",
        f"   prefix='{prefix}'", "   pseudo_dir='../../pseudo'",
        "   outdir='./qe_tmp'", "   disk_io='low'", "   verbosity='high'",
        "   tstress=.true.", "   tprnfor=.true.", "/",
        *system(path[stage_index]),
        *electrons(stage_index > 0, threshold, maxstep),
        *cards(False),
    ]
    return "\n".join(lines)


def md_input(prefix: str, final_lambda: float) -> str:
    lines = [
        "&CONTROL", "   calculation='md'", "   restart_mode='from_scratch'",
        f"   prefix='{prefix}'", "   pseudo_dir='../../pseudo'",
        "   outdir='./qe_tmp'", "   disk_io='low'", "   verbosity='high'",
        "   tstress=.true.", "   tprnfor=.true.",
        "   dt=20.6700d0", "   nstep=400", "/",
        *system(final_lambda),
        *electrons(True, 4.0e-3, 500),
        "&IONS", "   ion_dynamics='verlet'", "   pot_extrapolation='atomic'",
        "   wfc_extrapolation='none'", "   ion_temperature='svr'",
        "   ion_velocities='from_input'", "   tempw=4000.0",
        "   nraise=20", "/", *cards(True),
    ]
    return "\n".join(lines)


def main() -> None:
    for name in ("cases", "logs", "pseudo"):
        (ROOT / name).mkdir(parents=True, exist_ok=True)
    source_pseudo = ROUND10 / "pseudo" / PSEUDO
    if not source_pseudo.is_file():
        raise FileNotFoundError(source_pseudo)
    shutil.copy2(source_pseudo, ROOT / "pseudo" / PSEUDO)
    (ROOT / "atomic_velocities_4000K.txt").write_text("\n".join(VELOCITIES) + "\n")
    velocity_summary = {
        "target_temperature_K": 4000.0,
        "measured_temperature_K": 4000.0,
        "velocity_seed": 400054,
        "center_of_mass_velocity_au": [-7.765795845540367e-21, 4.65947750732422e-20, 5.436057091878257e-20],
        "degrees_of_freedom": 21,
        "source": "established thermalized BCC8 velocity set from Rounds 5-6",
    }
    (ROOT / "velocity_summary.json").write_text(json.dumps(velocity_summary, indent=2) + "\n")

    rows = []
    for serial_index, (case, path) in enumerate(CASES):
        folder = ROOT / "cases" / case
        folder.mkdir(exist_ok=True)
        if list(folder.glob("*.out")) or (folder / "qe_tmp").exists():
            raise RuntimeError(f"Refusing to overwrite results in {folder}")
        prefix = f"Fe_bcc8_r11_{case}"
        for stage_index, (stem, _, _) in enumerate(STAGES):
            (folder / f"{stem}.in").write_text(scf_input(prefix, path, stage_index))
        (folder / "md_400steps.in").write_text(md_input(prefix, path[-1]))
        row = {
            "serial_index": serial_index, "case": case, "lattice_a_A": 2.39,
            "target_temperature_K": 4000.0, "target_local_moment_Bohr": 2.0,
            "starting_magnetization": 0.125, "lambda_stage1": path[0],
            "lambda_stage2": path[1], "lambda_stage3_md": path[2],
            "md_steps": 400, "velocity_seed": 400054,
            "qe_requirement": "private patched QE 7.5 mcons=zv*starting_magnetization",
        }
        rows.append(row)
        (folder / "case_metadata.json").write_text(json.dumps(row, indent=2) + "\n")
    with (ROOT / "case_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (ROOT / "cases.txt").write_text("\n".join(row["case"] for row in rows) + "\n")
    print(f"Generated {len(rows)} Round-11 400-step case with thermalized positions and velocities")


if __name__ == "__main__":
    main()
