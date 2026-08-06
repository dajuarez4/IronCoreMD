#!/usr/bin/env python3
"""Generate ideal ferromagnetic BCC 4x4x4 QE inputs (128 Fe atoms)."""

from pathlib import Path


ROOT = Path(__file__).resolve().parent
CASES = {
    "a_2.50": (2.50, "a250"),
    "a_2.78": (2.78, "a278"),
}


def positions() -> str:
    lines = []
    for ix in range(4):
        for iy in range(4):
            for iz in range(4):
                lines.append(
                    f"Fe {ix / 4:.10f} {iy / 4:.10f} {iz / 4:.10f}"
                )
                lines.append(
                    f"Fe {(ix + 0.5) / 4:.10f} "
                    f"{(iy + 0.5) / 4:.10f} {(iz + 0.5) / 4:.10f}"
                )
    assert len(lines) == 128
    return "\n".join(lines)


def structure(a: float) -> str:
    length = 4.0 * a
    return f"""ATOMIC_SPECIES
Fe 55.845 Fe.pbe-spn-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS crystal
{positions()}

CELL_PARAMETERS angstrom
{length:.10f} 0.0000000000 0.0000000000
0.0000000000 {length:.10f} 0.0000000000
0.0000000000 0.0000000000 {length:.10f}
"""


def scf(a: float, tag: str) -> str:
    return f"""&CONTROL
    calculation = 'scf'
    prefix = 'fe_bcc128_fm_{tag}'
    pseudo_dir = '../pseudo'
    outdir = './tmp'
    tstress = .true.
    tprnfor = .true.
    disk_io = 'low'
/
&SYSTEM
    ibrav = 0
    nat = 128
    ntyp = 1
    ecutwfc = 71.0
    ecutrho = 496.0
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.02
    nspin = 2
    starting_magnetization(1) = 0.20
    nbnd = 1280
/
&ELECTRONS
    conv_thr = 1.0d-6
    electron_maxstep = 500
    mixing_beta = 0.10
    mixing_mode = 'local-TF'
    mixing_ndim = 12
    diagonalization = 'cg'
/
{structure(a)}
K_POINTS automatic
6 6 6 1 1 1
"""


def bands(a: float, tag: str) -> str:
    return f"""&CONTROL
    calculation = 'bands'
    prefix = 'fe_bcc128_fm_{tag}'
    pseudo_dir = '../pseudo'
    outdir = './tmp'
    disk_io = 'low'
/
&SYSTEM
    ibrav = 0
    nat = 128
    ntyp = 1
    ecutwfc = 71.0
    ecutrho = 496.0
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.02
    nspin = 2
    starting_magnetization(1) = 0.20
    nbnd = 1664
/
&ELECTRONS
    conv_thr = 1.0d-6
    electron_maxstep = 500
    startingpot = 'file'
    diagonalization = 'cg'
    diago_thr_init = 1.0d-6
    diago_cg_maxiter = 200
    diago_full_acc = .true.
/
{structure(a)}
K_POINTS crystal_b
8
0.000000 0.000000 0.000000 20
0.000000 0.500000 0.000000 20
0.500000 0.500000 0.000000 20
0.000000 0.000000 0.000000 20
0.500000 0.500000 0.500000 20
0.000000 0.500000 0.000000  0
0.500000 0.500000 0.000000 20
0.500000 0.500000 0.500000  0
"""


def nscf(a: float, tag: str) -> str:
    return f"""&CONTROL
    calculation = 'nscf'
    prefix = 'fe_bcc128_fm_{tag}'
    pseudo_dir = '../pseudo'
    outdir = './tmp'
    disk_io = 'low'
/
&SYSTEM
    ibrav = 0
    nat = 128
    ntyp = 1
    ecutwfc = 71.0
    ecutrho = 496.0
    occupations = 'tetrahedra'
    nspin = 2
    starting_magnetization(1) = 0.20
    nbnd = 1664
/
&ELECTRONS
    conv_thr = 1.0d-6
    electron_maxstep = 500
    startingpot = 'file'
    diagonalization = 'cg'
    diago_thr_init = 1.0d-6
    diago_cg_maxiter = 200
    diago_full_acc = .true.
/
{structure(a)}
K_POINTS automatic
8 8 8 0 0 0
"""


for dirname, (a, tag) in CASES.items():
    case = ROOT / dirname
    (case / "01_scf.in").write_text(scf(a, tag), encoding="utf-8")
    (case / "02_bands.in").write_text(bands(a, tag), encoding="utf-8")
    (case / "05_nscf.in").write_text(nscf(a, tag), encoding="utf-8")
