# Round 17: 16 atoms / 16 magnetic species on Jakar

This is the Jakar platform variant of Round 17. The physical calculation is
the same 16-atom intermediate-scale test prepared from successful Round 15:
one QE species and one constrained noncollinear spin direction per Fe atom.

## Calculation

- cubic BCC `2x2x2` cell, 16 atoms and `ntyp=16`;
- 16 distinct spin targets with numerically zero target-vector sum;
- target local moment 2.0 Bohr magnetons;
- lambda ramp `0.005 -> 0.020 -> 0.100`;
- SCF thresholds `4e-3 -> 1e-3 -> 3e-4 Ry`;
- Round-15 local-TF mixing and CG diagonalization;
- MD at lambda 0.100, 4000 K, `dt=20.67 au`, `nraise=20`;
- exactly 1000 MD steps;
- same pseudopotential and 71/496 Ry cutoffs as Round 15.

## Jakar executable and resources

This folder intentionally uses the existing modified Jakar executable:

```text
/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x
```

That executable already contains the `ntypx=128` change. No second QE build is
included here. The production script uses the requested Jakar configuration:

- account `jakar_general`;
- QoS `jakar_medium_general`;
- partition `medium`;
- one non-exclusive node allocation;
- 38 MPI ranks (`--ntasks-per-node=38`);
- Intel oneAPI from `/opt/intel/oneapi/setvars.sh --force`;
- `mpirun -np ${SLURM_NTASKS}`;
- `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `I_MPI_SHM_CMA=0`.

## Copy and submit

Copy this entire directory to Jakar, enter the copied directory, and verify the
required files:

```bash
ls cases cases.txt pseudo run_round17_jakar.sbatch
test -x /gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x
```

Submit from the Round-17 Jakar directory:

```bash
sbatch run_round17_jakar.sbatch
```

The Slurm log is written as `JOBID.log`. Monitor with:

```bash
squeue -u "$USER"
tail -f JOBID.log
bash check_round17_jakar.sh
```

The job runs the following stages sequentially in the same `qe_tmp` directory:

1. `stage1_4e-3.in`, lambda 0.005;
2. `stage2_1e-3.in`, lambda 0.020;
3. `stage3_3e-4.in`, lambda 0.100;
4. `md_1000steps.in`, lambda 0.100 and 1000 steps.

Each SCF stage must print `convergence has been achieved` before the next one
starts. Final success requires `round17_status.tsv` to say `DONE`, a `JOB DONE`
marker in `md_1000steps.out`, and exactly 1000 `Entering Dynamics` records.

The run refuses to start if the case already contains `.out` files or a
`qe_tmp` directory, preventing results from different attempts from being
mixed. Preserve a failed attempt under a separate directory before retrying.

## Provenance

- `generate_round17_jakar.py` recreates the inputs from committed Round 15.
- `structure_summary.json` records the tiled 16-atom construction.
- `magnetic_seed_summary.json` records all target directions.
- `velocity_summary.json` records the 4000 K, zero-COM velocity preparation.
- `case_manifest.csv` and `cases/.../case_metadata.json` record Jakar resources
  and the complete parameter recipe.
