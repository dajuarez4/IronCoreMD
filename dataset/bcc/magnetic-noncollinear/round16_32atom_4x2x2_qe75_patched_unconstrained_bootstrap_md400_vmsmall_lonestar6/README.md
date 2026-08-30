# Round 16: Round-14 32-atom workflow on Lonestar6 vm-small

## Purpose

Round 16 reproduces the complete physical workflow of Round 14 but runs it on
the Lonestar6 `vm-small` partition with eight MPI tasks instead of one `normal`
node with 64 tasks. It is intended to test whether the unconstrained bootstrap
can converge the difficult 32-atom charge/spin density with the smaller
parallel layout already used successfully for the eight-atom calculations.

This is not an independent thermodynamic replica: positions, velocities,
magnetic directions, artificial-species assignment, pseudopotential, random
seeds, SCF stages, and MD parameters are deliberately identical to Round 14.
Only the round-specific QE prefix and scheduler resources differ.

## Physical and numerical design

- BCC Fe `4x2x2` cell: 32 atoms, `9.56 x 4.78 x 4.78 angstrom`;
- thermal displacement source and selection seed `320889` from Round 14;
- Maxwell-Boltzmann velocities at 4000 K, seed `320400`;
- eight noncollinear directions in four antipodal pairs;
- `nat=32`, `ntyp=8`, with four atoms per artificial Fe species;
- two copies of every species on each BCC sublattice;
- private patched QE 7.5 with a maximum of ten species;
- Gamma point, `ecutwfc=71 Ry`, `ecutrho=496 Ry`;
- 400 MD steps, `dt=20.67`, SVR thermostat, `nraise=20`;
- final atomic constraint `lambda=0.100 Ry` and MD `conv_thr=4e-3 Ry`.

The SCF path is identical to Round 14:

1. unconstrained bootstrap: `degauss=0.08 Ry`, `conv_thr=2e-2 Ry`,
   `mixing_beta=0.05`;
2. constrained bootstrap: `lambda=0.001 Ry`, `conv_thr=1e-2 Ry`;
3. `lambda=0.005 Ry`, `conv_thr=4e-3 Ry`;
4. `lambda=0.020 Ry`, `conv_thr=1e-3 Ry`;
5. `lambda=0.100 Ry`, `conv_thr=3e-4 Ry`;
6. MD at `lambda=0.100 Ry`, `conv_thr=4e-3 Ry`.

Constrained stages use `mixing_mode='plain'`, `mixing_beta=0.01`, and a
40-vector history. The job stops at the first failed SCF.

## Lonestar6 resources

```text
partition = vm-small
nodes = 1
MPI tasks = 8
wall time = 48:00:00
```

The earlier 32-atom output estimated about `3.27 GB` total dynamical memory,
so the case should fit on `vm-small`. Eight tasks will probably be slower than
the 64-task Round-14 layout, and the 48-hour limit may still be reached. A
partial output is preserved if that happens.

The `--account=DMR26002` line from Round 14 is intentionally omitted, matching
the previously successful `vm-small` Round-11/Round-15 submission style. If
Lonestar6 requires an explicit allocation for the current account, submit with
`sbatch -A DMR26002 run_round16_lonestar6.sbatch` rather than editing the
scientific inputs.

## Copy and submit

From the Mac:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear
scp -r round16_32atom_4x2x2_qe75_patched_unconstrained_bootstrap_md400_vmsmall_lonestar6 \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/
```

On Lonestar6:

```bash
QE_FIX="$SCRATCH/apps/qe-7.5-atomic-fixed"
test -x "$QE_FIX/bin/pw.x" && echo "pw.x is executable"
sha256sum -c "$QE_FIX/pw.x.sha256"

cd /work/11381/dajuarez4/ls6/round16_32atom_4x2x2_qe75_patched_unconstrained_bootstrap_md400_vmsmall_lonestar6
mkdir -p logs
find cases \( -name '*.out' -o -name qe_tmp \) -print
sbatch run_round16_lonestar6.sbatch
```

The `find` command must produce no output before the first submission. Do not
reuse Round-12, Round-13, or Round-14 `qe_tmp` data.

## Monitor

```bash
squeue -u "$USER"
bash check_round16.sh
bash plot_round16.sh --live 30
```

For job ID `JOB_ID`:

```bash
tail -F logs/ls6_round16_vm_md400.JOB_ID.out
tail -F logs/ls6_round16_vm_md400.JOB_ID.err
sacct -j JOB_ID --format=JobID,JobName,Partition,AllocCPUS,State,Elapsed,ExitCode
```

## Acceptance criteria

1. all five SCF stages converge in order;
2. MD contains exactly 400 `Entering Dynamics` records and `JOB DONE`;
3. no QE-routine, Cholesky, `cdiaghg`, or failed-SCF error occurs;
4. every MD SCF converges below `4e-3 Ry` without approaching 500 iterations;
5. the trailing mean local moment remains preferably above `1 Bohr magneton/Fe`;
6. the time-averaged total magnetization per atom remains small;
7. temperature, energy, force, pressure, and structure stay bounded.

If the unconstrained bootstrap fails, later stages and MD must not be used. If
it converges but a constrained stage fails, the failure identifies the point
in the lambda ramp where the 32-atom constraint becomes unstable.

## Lightweight results

After a complete run:

```bash
module reset
module load python3/3.9.7
python3 collect_round16_npz.py --require-complete
ls -lh round16_results_light.npz round16_md400_dashboard.png
```

Do not copy or commit `qe_tmp`, full `.out` files, or scheduler logs. Preserve
a completed directory unchanged and use a new clean directory for any rerun.
