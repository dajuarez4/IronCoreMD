# Round 17: 16 species on Lonestar6 `vm-small`

This is the Lonestar6 platform variant of the same Round-17 calculation:
16 atoms, 16 QE species, 16 distinct constrained spin directions, the
successful Round-15 SCF recipe, and 1000 MD steps at 4000 K.

## Resources

The production script requests one `vm-small` VM, 16 MPI tasks, and the maximum
48-hour wall time. It uses allocation `DMR26002`, `ibrun`, and the same
`intel/19.1.1 impi/19.0.9` environment used by the working Lonestar6 rounds.

`vm-small` is intended for jobs using at most 16 cores and less than 29 GB of
memory. Confirm current queue policy with `qlimits` immediately before
submission because TACC can change limits.

## Required QE executable

The 16-species input requires the Lonestar6-local private executable:

```text
$SCRATCH/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x
```

The older `$SCRATCH/apps/qe-7.5-atomic-fixed/bin/pw.x` has the magnetic
constraint correction but should not be assumed to support more than 10
species. The Stampede3 and Jakar binaries are not used here.

If the ntypx=128 executable already exists on Lonestar6, skip directly to the
production submission. Otherwise copy this full folder to Lonestar6 and build
it once.

## One-time Lonestar6 build

From `qe_build_lonestar6`:

```bash
mkdir -p "$SCRATCH/src" "$SCRATCH/apps"
wget -O "$SCRATCH/src/q-e-qe-7.5.tar.gz" \
  https://gitlab.com/QEF/q-e/-/archive/qe-7.5/q-e-qe-7.5.tar.gz

bash prepare_qe75_ntyp128_source.sh \
  "$SCRATCH/src/q-e-qe-7.5.tar.gz"

sbatch build_qe75_ntyp128_lonestar6.sbatch
```

When the build completes, run the 16-species smoke test:

```bash
sbatch test_qe75_ntyp128_lonestar6.sbatch
```

The test log must end with:

```text
PASS: patched pw.x provenance confirms ntypx=128 and QE parsed all 16 species.
```

## Submit production

From the Round-17 Lonestar6 folder:

```bash
qlimits
sbatch run_round17_lonestar6_vmsmall.sbatch
```

The Slurm output is `JOBID.log`. QE writes the calculation outputs inside:

```text
cases/a2p39_2x2x2_m2_16species_lambda0p100_md1000/
```

Monitor with:

```bash
squeue -j JOBID
tail -f JOBID.log
bash check_round17_lonestar6.sh
```

The production script is intentionally short: it loads the established
Lonestar6 modules and runs stage1, stage2, stage3, and MD sequentially with
`ibrun`. `set -e` prevents later stages from starting when an MPI command exits
with an error.

## Calculation recipe

- SCF lambda ramp: `0.005 -> 0.020 -> 0.100`.
- SCF thresholds: `4e-3 -> 1e-3 -> 3e-4 Ry`.
- Local-TF mixing, beta 0.010, mixing dimension 20, CG diagonalization.
- MD lambda 0.100, SCF threshold 4e-3 Ry.
- SVR thermostat at 4000 K, `nraise=20`, `dt=20.67 au`.
- 1000 MD steps, Gamma point, 71/496 Ry cutoffs.
- Fe01 through Fe16 use the same pseudopotential and distinct target spins.

Provenance is recorded in `case_manifest.csv`, `structure_summary.json`,
`magnetic_seed_summary.json`, `velocity_summary.json`, and the case metadata.

## Result copied from Lonestar6 (2026-08-28)

The root-level `md_1000steps.out` is a usable but incomplete trajectory. It
contains 392 completed dynamics records (about 0.392 ps) out of the requested
1000 and has no `JOB DONE`. All 392 completed ionic SCFs report convergence;
there is no `convergence NOT achieved`, `Error in routine`, Cholesky, or
`cdiaghg` failure. The file stops during the first electronic iteration of the
next MD step, which is consistent with an external interruption or wall-time
limit rather than a QE fatal error. There are 229 nonfatal `c_bands`
eigenvalue-convergence warnings.

The last ionic state printed after dynamics step 392 has a temperature of
1570.04 K. Across the 392 dashboard configurations, the local-moment
magnitudes span 1.474--1.687 Bohr magnetons per Fe and the largest measured
spin-direction change relative to configuration 0 is only 3.09 degrees. The
small physical rotation is why the arrows appear nearly stationary even
though every GIF frame uses the corresponding measured vector.

The corrected dashboard products are:

```text
round17_vmsmall_component_dashboard.gif
round17_vmsmall_component_dashboard_final.png
```

They show all 16 atoms inside a periodically shifted 4.78-Angstrom cubic cell.
The ideal-site displacement calculation matches the thermally displaced input
atoms to their nearest periodic BCC sites before computing distances. Rebuild
the 99-frame, stride-4 GIF with:

```bash
/opt/anaconda3/bin/python3 make_round17_component_dashboard.py
```
