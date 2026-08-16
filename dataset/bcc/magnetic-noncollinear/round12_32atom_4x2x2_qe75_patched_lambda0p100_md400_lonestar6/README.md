# Round 12: 32-atom BCC 4x2x2 noncollinear MD400

## Purpose

Round 12 scales the validated patched-QE noncollinear setup from the 8-atom
`2x2x1` diagnostic to a 32-atom BCC `4x2x2` cell at `a = 2.39 angstrom`. It
runs 400 MD steps (about 200 fs) at a nominal 4000 K.

This is the first finite-size scaling test. It is substantially more useful
than the 8-atom cell for checking cancellation of disordered local moments,
but 400 steps in one 32-atom trajectory are still not production statistics.

## Model and starting state

- cell: `9.56 x 4.78 x 4.78 angstrom`, 32 Fe atoms;
- genuinely displaced starting positions, not an ideal BCC lattice;
- displacement pool taken from frame 28 (QE MD iteration 28) of
  `dataset/bcc/non-mag/2.39_4000K.npz`, where the recorded temperature is
  about 3972 K;
- selected RMS displacement: recorded in `thermalized_start_summary.json`;
- independent Maxwell-Boltzmann velocities at exactly 4000 K, seed `320400`,
  with center-of-mass motion removed;
- eight noncollinear directions forming four antipodal pairs, each repeated
  four times, so the starting magnetic-vector sum is zero;
- atomic target `2 Bohr magneton/Fe`, corrected private QE 7.5 constraint;
- `lambda: 0.005 -> 0.020 -> 0.100 Ry`;
- SCF thresholds: `4e-3 -> 1e-3 -> 3e-4 Ry`;
- MD `conv_thr = 4e-3 Ry`, `nraise = 20`, `dt = 20.67 a.u.`, 400 steps;
- Gamma point, Fermi-Dirac smearing, `degauss = 0.025334 Ry`.

QE 7.5 has a ten-species limit. Therefore this package intentionally uses
`ntyp=8`, not one artificial species per atom. Every magnetic direction is
shared by four atoms. This preserves the patched constraint and exact initial
vector cancellation without requiring another private QE modification, but it
also imposes a repeated eight-direction magnetic pattern.

Round 11 showed that `nraise=20` did not maintain a 4000 K mean in the 8-atom
cell. The larger number of ionic degrees of freedom may reduce fluctuations,
but the trailing temperature must still be checked before calling this a
4000 K trajectory.

## Regeneration

The inputs are already generated. To reproduce them locally from the full
repository:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round12_32atom_4x2x2_qe75_patched_lambda0p100_md400_lonestar6
/opt/anaconda3/bin/python generate_round12.py
```

The generator refuses to overwrite a directory containing outputs or
`qe_tmp`.

## Copy to Lonestar6

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear
scp -r round12_32atom_4x2x2_qe75_patched_lambda0p100_md400_lonestar6 \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/
```

## Verify and submit

On Lonestar6:

```bash
QE_FIX="$SCRATCH/apps/qe-7.5-atomic-fixed"
test -x "$QE_FIX/bin/pw.x" && echo "pw.x is executable"
sha256sum -c "$QE_FIX/pw.x.sha256"

cd /work/11381/dajuarez4/ls6/round12_32atom_4x2x2_qe75_patched_lambda0p100_md400_lonestar6
mkdir -p logs
find cases \( -name '*.out' -o -name qe_tmp \) -print
sbatch run_round12_lonestar6.sbatch
```

The `find` command must print nothing before the first submission. The job uses
one `normal` node, 64 MPI ranks, and a 48-hour wall-time request. The 32-atom
run is much heavier than Round 11 and could reach the wall-time limit; the
script will preserve a partial output if this happens.

## Monitor

```bash
squeue -u "$USER"
bash check_round12.sh
bash plot_round12.sh --live 30
```

Replace `JOB_ID` below with the submitted job number:

```bash
tail -F logs/ls6_round12_md400.JOB_ID.out
tail -F logs/ls6_round12_md400.JOB_ID.err
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,ExitCode
```

## Acceptance criteria

1. all three SCF stages converge;
2. 400 `Entering Dynamics` records and `JOB DONE` are present;
3. every completed MD SCF reaches `4e-3 Ry` without approaching 500 iterations;
4. no fatal QE, Cholesky, `cdiaghg`, or `efermig` error occurs;
5. the trailing mean local moment remains preferably above `1 Bohr magneton/Fe`;
6. the time-averaged total magnetization per atom stays small;
7. the final 100-step temperature mean is reasonably close to 4000 K;
8. energy, force, and pressure remain bounded without systematic instability.

## Lightweight archive

After completion:

```bash
module reset
module load python3/3.9.7
python3 collect_round12_npz.py --require-complete
ls -lh round12_results_light.npz round12_md400_dashboard.png
```

Copy back only the lightweight results:

```bash
scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round12_32atom_4x2x2_qe75_patched_lambda0p100_md400_lonestar6/round12_results_light.npz \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round12_32atom_4x2x2_qe75_patched_lambda0p100_md400_lonestar6/

scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round12_32atom_4x2x2_qe75_patched_lambda0p100_md400_lonestar6/round12_md400_dashboard.png \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round12_32atom_4x2x2_qe75_patched_lambda0p100_md400_lonestar6/
```

Do not copy or commit `qe_tmp`, full `.out` files, or scheduler logs.
