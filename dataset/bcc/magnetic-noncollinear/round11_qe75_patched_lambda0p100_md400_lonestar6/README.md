# Round 11: patched-QE lambda=0.10, 400-step MD

## Purpose

Round 10 completed 50 MD steps with stable SCF convergence and retained local
moments, but its temperature was not yet equilibrated around 4000 K. Round 11
extends the same validated 8-atom diagnostic to 400 steps (about 200 fs). It
uses the user-selected `nraise = 20`; with `dt = 20.67 a.u.` (about 0.5 fs), the
SVR thermostat time is approximately 10 fs.

This remains a small-cell convergence and stability test. It is not a
production trajectory or an independent thermodynamic sample.

## Fixed design

- BCC Fe 8-atom `2x2x1` cell, `a = 2.39 angstrom`;
- thermally displaced starting positions;
- the established 4000 K seed-`400054` velocities;
- noncollinear atomic target `2 Bohr magneton/Fe`;
- corrected private QE 7.5 atomic constraint, final `lambda = 0.10 Ry`;
- fresh SCF ramp `0.005 -> 0.020 -> 0.100 Ry`;
- SCF thresholds `4e-3 -> 1e-3 -> 3e-4 Ry`;
- MD `conv_thr = 4e-3 Ry` for exploratory dynamics;
- Gamma point, Fermi-Dirac smearing, `degauss = 0.025334 Ry`;
- 400 MD steps with `dt = 20.67 a.u.` and `nraise = 20`.

Each run builds a fresh `qe_tmp` through the three SCF stages. Do not copy a
restart directory from an earlier round.

## Copy to Lonestar6

From the Mac:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear
scp -r round11_qe75_patched_lambda0p100_md400_lonestar6 \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/
```

## Verify the private patched QE

On Lonestar6:

```bash
QE_FIX="$SCRATCH/apps/qe-7.5-atomic-fixed"
test -x "$QE_FIX/bin/pw.x" && echo "pw.x is executable"
sha256sum -c "$QE_FIX/pw.x.sha256"
cat "$QE_FIX/ATOMIC_CONSTRAINT_PATCH_INFO.txt"
```

The checksum must report `OK`.

## Submit and monitor

```bash
cd /work/11381/dajuarez4/ls6/round11_qe75_patched_lambda0p100_md400_lonestar6
mkdir -p logs
find cases \( -name '*.out' -o -name qe_tmp \) -print
sbatch run_round11_lonestar6.sbatch
```

The `find` command must print nothing before the first submission. Monitor with:

```bash
squeue -u "$USER"
bash check_round11.sh
bash plot_round11.sh --live 30
```

For the scheduler record, replace `JOB_ID` below:

```bash
tail -F logs/ls6_round11_md400.JOB_ID.out
tail -F logs/ls6_round11_md400.JOB_ID.err
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,ExitCode
```

## Interpretation and acceptance

Use roughly steps 1--100 as equilibration. Analyze steps 101--400 only if the
trailing temperature and moment curves have become bounded and reasonably
stationary. A successful technical validation must:

1. reach 400 `Entering Dynamics` records and `JOB DONE`;
2. contain no failed SCF, QE routine, Cholesky, or `cdiaghg` error;
3. keep converged MD cycles below `conv_thr = 4e-3 Ry`;
4. avoid electronic cycles approaching `electron_maxstep = 500`;
5. retain a trailing mean local moment preferably above `1 Bohr magneton/Fe`;
6. show bounded temperature, force, pressure, energy, and moment behavior.

The 8-atom cell is useful for this test but is too small for production-quality
magnetic or thermodynamic statistics.

## Lightweight results

After completion:

```bash
module reset
module load python3/3.9.7
python3 collect_round11_npz.py --require-complete
ls -lh round11_results_light.npz round11_md400_dashboard.png
```

Copy back only the lightweight products:

```bash
scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round11_qe75_patched_lambda0p100_md400_lonestar6/round11_results_light.npz \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round11_qe75_patched_lambda0p100_md400_lonestar6/

scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round11_qe75_patched_lambda0p100_md400_lonestar6/round11_md400_dashboard.png \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round11_qe75_patched_lambda0p100_md400_lonestar6/
```

Do not copy or commit `qe_tmp`, full `.out` files, or scheduler logs. For a
rerun, preserve the completed directory and copy a clean repository version to
a new directory.

## Completed result

Round 11 reached all 400 MD steps and `JOB DONE`. All 400 MD SCF cycles
converged below `4e-3 Ry`; they required a mean 6.28 electronic iterations and
a maximum of 12. No fatal QE, Cholesky, `cdiaghg`, or `efermig` error occurred.
There were 201 recoverable `c_bands` warnings.

The mean local moment remained stable at `1.595 Bohr magneton/Fe`, while the
norm of the full-trajectory mean total-magnetization vector was only
`0.0527 Bohr magneton/cell` (`0.0066` per Fe). The trajectory is therefore
consistent with a small-cell disordered-local-moment paramagnetic state.

Temperature did not equilibrate at the requested 4000 K: the trajectory mean
was about 2945 K and the final-100-step mean was about 2977 K. Thus Round 11
validates electronic and magnetic stability, but it must not be used as an
equilibrated 4000 K thermodynamic trajectory.
