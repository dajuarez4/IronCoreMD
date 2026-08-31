# Round 19: 48 atoms and 48 spin species on Jakar

This round scales the successful Round-17 setup to 48 atoms on Jakar. It uses
the Jakar QE executable that already contains the corrected atomic constraint
and `ntypx=128`:

```text
/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x
```

## Physical construction

- BCC Fe `6x2x2` supercell: 48 atoms, `14.34 x 4.78 x 4.78 A`.
- The Round-17 `2x2x2` initial state is tiled three times along x.
- One artificial Fe species per atom: `Fe01` through `Fe48`.
- 48 distinct noncollinear target directions in six antipodal groups.
- Each spatial tile has a zero target-spin vector sum.
- The second and third velocity tiles use deterministic row and Cartesian-axis
  permutations; center-of-mass motion is removed and the complete set is
  rescaled to exactly 4000 K.
- All 48 species use the same Fe PAW pseudopotential.

## QE recipe

- SCF ramp: `lambda=0.005, conv_thr=4e-3 Ry`;
  `lambda=0.020, conv_thr=1e-3 Ry`;
  `lambda=0.100, conv_thr=3e-4 Ry`.
- `mixing_mode='local-TF'`, `mixing_beta=0.010`, `mixing_ndim=20`.
- CG diagonalization, Gamma point, 71/496 Ry cutoffs.
- MD: 1000 steps, `dt=20.67 au`, SVR thermostat at 4000 K,
  `nraise=20`, `lambda=0.100`, `conv_thr=4e-3 Ry`.

## Jakar resources

The short submission script uses the general medium allocation and 38 MPI
tasks, without preflight pseudopotential or provenance checks:

```text
account = jakar_general
qos = jakar_medium_general
partition = medium
nodes = 1
MPI tasks = 38
```

## Copy and submit

Copy the complete folder to Jakar, then run:

```bash
cd /gpfs/scratch/dajuarez4/testing/non_collinear/round19_48atom_6x2x2_qe75_ntyp128_lambda0p100_md1000_jakar
sbatch run_round19_jakar.sbatch
```

Monitor with:

```bash
squeue -u dajuarez4
bash check_round19_jakar.sh
tail -f JOBID.log
```

If Slurm reports immediate `FAILED 0:53` and creates no log, the batch script
did not start; that reproduces the earlier Jakar scheduler/node launch problem
and is not evidence of a QE input failure. If the script begins, the first QE
output is `cases/.../stage1_4e-3.out`.

Do not mix old `*.out` files or a previous `qe_tmp` directory into this round.

## Partial result checked 2026-08-30

The returned root-level `md_1000steps.out` is incomplete. Its SHA-256 is
`573a595cdeab18203ad7b376c1123cfa5a678076a6b51d8f99a0d36cb01674ca`.
It contains 59 of the requested 1000 converged MD SCFs and 59 `Entering
Dynamics` records, but no `JOB DONE`. The file ends during the following SCF
after its first electronic iteration. There are no QE-routine, Cholesky,
`cdiaghg`, or explicit electronic-convergence failure markers in the returned
text. QE printed one recoverable `c_bands` warning.

The short segment already exposes several concerns:

- mean ionic temperature: `5387.51 K`, but the thermostat has turned the
  trajectory downward; the final value is `3677.21 K` and the final-10 mean is
  `3769.92 K`, now reasonably close to the `4000 K` target;
- mean pressure: `353.46 GPa`; the final value is `342.67 GPa`;
- mean printed local moment: only `0.10715 Bohr magneton/Fe`, with a final
  value of `0.10537`;
- more decisively, every printed constrained vector has magnitude `0.125 Bohr
  magneton`, not the requested `2.0 Bohr magneton`. Rounds 15 and 17 print
  magnitude `2.0` from the same `starting_magnetization=0.125` input because
  their corrected builds apply `mcons = Zval * starting_magnetization`;
- mean final SCF accuracy: `2.914e-3 Ry`, with a maximum of `3.9903e-3 Ry`,
  just inside the `4e-3 Ry` MD threshold;
- mean electronic iterations per completed MD SCF: `16.61`, maximum `46`;
- the printed SCF-force correction averages `0.211` of the total-force
  magnitude and triggers the large-correction warning in 57 of 59 steps;
- minimum periodic Fe--Fe separation: `1.654 A`; maximum RMS displacement
  from the starting configuration: `0.533 A`.

No quantitative phonon dispersion, phonon DOS, or TDEP force constants can be
obtained from this file. A 59-configuration, `0.059 ps` interrupted trajectory
cannot support a defensible TDEP fit or velocity-autocorrelation spectrum, and
the loose-SCF MD forces are not suitable for quantitative fitting. The HELD
section below therefore presents only the requested per-frame diagnostic
animation. More fundamentally, this trajectory used the wrong constraint
scale: the Jakar executable was not the expected corrected atomic-constraint
build, or that correction was inactive. Do not continue or use this fit
quantitatively. Verify the binary with a short test that prints a `2.0 Bohr
magneton` constrained-vector magnitude, then rerun from a clean directory. For
quantitative TDEP, recalculate selected, decorrelated snapshots with a
substantially tighter force SCF threshold.

`round19_partial_dashboard.png` summarizes every completed record and includes
an explicit phonon-readiness panel. The machine-readable audit is
`round19_analysis_summary.json`, and aligned per-step values are in
`round19_md_records.csv`. Regenerate all three with:

```bash
/opt/anaconda3/bin/python analyze_round19.py
```

## HELD animated structure and phonon dashboard

`round19_held_md_dashboard.gif` follows the visual format of the repository's
`assets/bcc_noncollinear_md_dashboard.gif`. For each available force-aligned
frame it synchronizes:

- the complete moving `6x2x2`, 48-atom structure and ideal BCC reference sites;
- atom-resolved QE local-moment arrows and the global magnetization vector;
- magnetization, temperature, and displacement histories;
- the frame-specific HELD dispersion along `Gamma-H-N-Gamma-P-H`.

The current GIF has 59 frames at 5 frames/s and is `2550x1500` pixels. The HELD
frequencies span approximately `-22.12` to `38.70 THz`; the negative branches
are shown, not suppressed. These are qualitative diagnostic branches from the
incomplete, incorrectly constrained, loose-force trajectory described above,
not validated quantitative phonons.

The complete reproducible entry point is `run_round19_held_dashboard.sh`. It
checks the Python dependencies, rebuilds the trajectory NPZ, reruns HELD,
writes the mean dispersion and heatmap, and renders both the final PNG and GIF:

```bash
bash run_round19_held_dashboard.sh
```

For a future 1000-step output, sampling every tenth frame keeps the GIF near
the length used for earlier rounds while always including the final frame:

```bash
ROUND19_STRIDE=10 bash run_round19_held_dashboard.sh
```

The Python executable, frame stride, frame rate, DPI, and QE output can also be
overridden with `ROUND19_PYTHON`, `ROUND19_STRIDE`, `ROUND19_FPS`,
`ROUND19_DPI`, and `ROUND19_QE_OUTPUT`.
