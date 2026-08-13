# Exploratory BCC Fe MD at 4000 K: 10-step chunk

This folder implements the accepted exploratory setting `conv_thr=4e-3 Ry`.
It is intentionally limited to 10 steps so convergence, magnetization, and
temperature can be reviewed before any longer trajectory is attempted.

This is not production-quality AIMD. QE reported SCF force corrections of the
same order as the forces in the successful five-step pilot. See
`../EXPLORATORY_MD_DECISION.md` for the evidence and scientific limitations.

## Settings

- 16-atom BCC cell, explicit Gamma point
- noncollinear atomic-direction constraint, `lambda=0.020 Ry`
- `conv_thr=4e-3 Ry`
- local-TF mixing, beta 0.010, history 20
- CG diagonalization with `diago_thr_init=1e-3`
- atomic potential extrapolation and no wavefunction extrapolation
- 4000 K SVR thermostat, 0.5 fs timestep, 10 steps

## 1. Transfer the complete folder

Place it beside the existing Jakar `md_good` directory:

```text
bcc2_pm_4000K_convergence/
├── md_good/
└── md_exploratory_lambda0p020_conv4e-3_10steps/
```

## 2. Copy the clean preconditioned state

From inside this folder on Jakar:

```bash
cp -R ../md_good/preconditioned_lambda0p020/Fe_bcc2_pm_round4.save seed_state/
cp ../md_good/preconditioned_lambda0p020/Fe_bcc2_pm_round4.wfc* seed_state/
```

Verify the protected seed:

```bash
test -d seed_state/Fe_bcc2_pm_round4.save && echo "SAVE OK"
ls seed_state/Fe_bcc2_pm_round4.wfc* | wc -l
```

The expected wavefunction count is 40. Never copy the evolved `md_good/md_state`
into `seed_state`.

## 3. Stage the Fe pseudopotential

If this workflow is still inside the repository, copy the existing curated PAW
file without adding another copy to Git:

```bash
cp ../../../../dataset/hcp/hcp_mag/pseudo/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF pseudo/
```

If the workflow was copied to a standalone scratch directory, copy the same
pseudopotential from your trusted Jakar location into `pseudo/`.

## 4. Submit

The run directory must not already contain `md_state` or
`md_exploratory_4000K_10steps.out`.

```bash
sbatch run_exploratory10_jakar.sbatch
```

## 5. Monitor

```bash
squeue -u "$USER"
bash check_exploratory10.sh
```

For a live filtered view:

```bash
tail -f md_exploratory_4000K_10steps.out | \
grep --line-buffered -E "Entering Dynamics|iteration #|estimated scf accuracy|total magnetization|absolute magnetization|temperature|convergence has been achieved|c_bands|Cholesky|Error in routine|JOB DONE"
```

Press `Ctrl+C` to stop viewing; this does not stop QE.

## 6. Inspect after completion

```bash
bash check_exploratory10.sh
python3 extract_exploratory10.py
```

Accept the chunk only if it reaches `JOB DONE`, has no Cholesky failure, all
electronic cycles converge, temperature remains controlled, and converged net
magnetization does not run toward 5 mu_B/cell. Archive the folder before making
another attempt.
