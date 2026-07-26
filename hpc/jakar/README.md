# Jakar workflow

The settings below are based on the repository's current `run_jakar.sbatch`. Verify the account, QoS, partition, and task limits before each campaign because cluster policy can change.

## Current environment

The existing run uses:

- account: `jakar_partner_physj`
- QoS: `jakar_medium_physj`
- partition: `medium`
- one node with 20 MPI tasks
- Intel oneAPI initialized with `/opt/intel/oneapi/setvars.sh --force`
- Intel MPI launched with `mpirun`
- existing QE binary: `/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x`

## Install QE 7.5

On the designated login or build node:

```bash
module purge
source /opt/intel/oneapi/setvars.sh --force

which mpiifx
which mpirun

QE_ARCHIVE=/gpfs/scratch/$USER/software/qe-7.5-ReleasePack.tar.gz \
QE_INSTALL_ROOT=/gpfs/scratch/$USER/software \
BUILD_JOBS=8 \
bash /path/to/IronCoreMD/hpc/install_qe_from_source.sh
```

For a new out-of-source build, `pw.x` should be under:

```text
/gpfs/scratch/$USER/software/qe-7.5/build/bin/pw.x
```

The older in-source installation already used by this project is under `/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x`; both layouts are valid if the executable and loaded MPI environment match.

## Prepare a run

```bash
mkdir -p /gpfs/scratch/$USER/iron_runs/RUN_NAME/pseudo
cd /gpfs/scratch/$USER/iron_runs/RUN_NAME

cp /path/to/calculation.in input.in
cp /path/to/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF pseudo/
cp /path/to/IronCoreMD/hpc/jakar/qe_pw_template.sbatch run.sbatch
```

Edit the SBATCH resources and export `QE`, `INPUT`, and `OUTPUT` either in `run.sbatch` or at submission time.

## Submit and monitor

```bash
sbatch run.sbatch
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,ExitCode,MaxRSS
tail -f qe.out
```

Cancel only when necessary:

```bash
scancel JOB_ID
```

## Stage 1 to stage 2

For the noncollinear SCF-to-MD workflow:

1. Run stage 1 with `calculation='scf'` and `restart_mode='from_scratch'`.
2. Confirm that the output contains `convergence has been achieved` and `JOB DONE`.
3. Keep the complete `outdir`, including the `PREFIX.save` directory and wavefunction files.
4. Run stage 2 from the same directory with the same `prefix` and `outdir`.
5. Use `startingpot='file'` and `startingwfc='file'` only when those restart files exist and are readable.

For an interrupted MD continuation, use a separate restart input with `restart_mode='restart'`. Do not reinsert the original `ATOMIC_VELOCITIES` or replace the saved ionic state.

## Production record

Copy each actual SBATCH script into `hpc/jakar/jobs/` with a descriptive name such as:

```text
Fe_bcc_3x3x3_noncollinear_stage1_4000K.sbatch
Fe_bcc_3x3x3_noncollinear_stage2_4000K.sbatch
```
