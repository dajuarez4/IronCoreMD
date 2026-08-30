# Round 17: 16 atoms / 16 magnetic species on Stampede3

This is the deliberate intermediate step between the successful 8-atom
Round 15 and the unsuccessful 32-atom attempts. It uses a cubic BCC `2x2x2`
cell with 16 Fe atoms and **16 QE species**, one species and one constrained
spin direction per atom. It does not reuse only eight species.

## What is unchanged from successful Round 15

- QE 7.5 plus the `mcons = zv * starting_magnetization` private patch.
- Target local moment 2.0 Bohr magnetons (`starting_magnetization=0.125`).
- SCF ramp: lambda `0.005 -> 0.020 -> 0.100`.
- SCF thresholds: `4e-3 -> 1e-3 -> 3e-4 Ry`.
- `mixing_mode='local-TF'`, `mixing_beta=0.010`, `mixing_ndim=20`, and CG.
- MD: lambda `0.100`, `conv_thr=4e-3 Ry`, `dt=20.67 au`, SVR at 4000 K,
  `nraise=20`, and exactly 1000 steps.
- 71/496 Ry cutoffs, Fermi-Dirac smearing 0.025334 Ry, and Gamma point.

The first eight directions are exactly Round 15. The other eight are rotated
by +47 degrees in azimuth, giving 16 distinct targets while keeping the target
vector sum numerically zero. Positions are Round 15 tiled once along z. The
upper-tile velocities use a deterministic atom/axis permutation of the Round
15 set; the combined velocities have zero center-of-mass velocity and are
rescaled exactly to 4000 K. This avoids cloning identical ionic trajectories.

## Why a private QE build is required

Stock QE 7.5 defines `ntypx=10`, so an input with `ntyp=16` cannot run in that
binary. The patch in `qe_build_stampede3/qe75_atomic_constraint_ntyp128.patch`
raises the compile-time limit to 128 and also includes the same atomic magnetic
constraint correction used for Round 15. Fe labels are `Fe01` through `Fe16`;
all are six characters or fewer.

## Build and prove the executable on Stampede3

Copy this entire Round-17 directory to Stampede3. From its
`qe_build_stampede3` directory:

```bash
mkdir -p "$SCRATCH/src" "$SCRATCH/apps"
wget -O "$SCRATCH/src/q-e-qe-7.5.tar.gz" \
  https://gitlab.com/QEF/q-e/-/archive/qe-7.5/q-e-qe-7.5.tar.gz
bash prepare_qe75_ntyp128_source.sh "$SCRATCH/src/q-e-qe-7.5.tar.gz"
sbatch -A YOUR_TACC_ALLOCATION build_qe75_ntyp128_stampede3.sbatch
```

After the build job finishes successfully, run the mandatory 16-species smoke
test:

```bash
sbatch -A YOUR_TACC_ALLOCATION test_qe75_ntyp128_stampede3.sbatch
```

Do not submit production until its log ends with:

```text
PASS: patched pw.x provenance confirms ntypx=128 and QE parsed all 16 species.
```

QE 7.5 formats the reported maximum with a two-character integer field, so
`ntypx=128` appears as `**` in `test_ntyp16.out`. The test verifies the exact
value from the installed patch-provenance file and independently verifies that
the executable parsed and ran all 16 species through `JOB DONE`.

The installed executable is:

```text
$SCRATCH/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x
```

The build and test scripts use `intel/24.0`, `impi/21.11`, `mpif90`, `ibrun`,
and an SKX node. If `module spider intel impi` shows different versions at
submission time, change both scripts and the production script consistently.

## Submit Round 17

From the Round-17 directory (not from `cases/`):

```bash
sbatch -A YOUR_TACC_ALLOCATION run_round17_stampede3.sbatch
```

The production job requests one SKX node, 16 MPI ranks, and 48 hours. It first
checks the installed executable checksum and provenance, then runs all three
SCF stages in sequence. MD starts only if every stage reports convergence. It
also refuses to mix a new attempt with existing `.out` files or `qe_tmp`.

Monitor with:

```bash
bash check_round17.sh
tail -f logs/stampede3_round17_md1000.JOBID.out
```

Success requires all of the following:

- `round17_status.tsv` says `DONE`;
- `md_1000steps.out` contains `JOB DONE`;
- exactly 1000 `Entering Dynamics` records;
- no QE fatal error, Cholesky/cdiaghg failure, or unconverged electronic step.

## Regeneration and provenance

`generate_round17.py` reads the committed Round-15 MD input and recreates all
four inputs, summaries, manifest, and velocities. It refuses to overwrite a
case after output/state exists. The pseudopotential is the same file as Round
15. Important provenance files are:

- `structure_summary.json`
- `velocity_summary.json`
- `magnetic_seed_summary.json`
- `cases/.../case_metadata.json`
- `$SCRATCH/apps/qe-7.5-atomic-fixed-ntyp128/ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt`

If the 16-atom run fails, diagnose which SCF ramp stage fails before changing
parameters. The point of this round is to change only system size/species count
relative to the validated Round-15 recipe; adding a new bootstrap scheme at the
same time would obscure the cause.
