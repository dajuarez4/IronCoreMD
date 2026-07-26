# Lonestar6 workflow

Lonestar6 uses Slurm and the TACC-specific `ibrun` launcher for MPI applications. The account, queue, task count, and currently available modules must be checked on the system rather than copied from Jakar.

## Discover the available environment

Run these commands after logging in:

```bash
module purge
module spider quantum-espresso
module spider intel
module spider gcc
module spider impi
module spider openmpi
module avail 2>&1 | less

idev -h
sinfo
```

If TACC provides a suitable QE module, record its exact name and version in the production SBATCH file. Otherwise, build QE with one compiler/MPI stack and load the same stack in every job.

## Install QE 7.5 from source

Choose and load a supported compiler/MPI stack shown by `module spider`, then run:

```bash
module purge
module load COMPILER_MODULE
module load MPI_MODULE

which mpif90
which ibrun

QE_ARCHIVE=$SCRATCH/software/qe-7.5-ReleasePack.tar.gz \
QE_INSTALL_ROOT=$SCRATCH/software \
BUILD_JOBS=8 \
bash /path/to/IronCoreMD/hpc/install_qe_from_source.sh
```

Do not copy the Jakar executable to Lonestar6. Recompile on Lonestar6 so QE links against the compiler, MPI, and numerical libraries available there.

## Prepare and submit

```bash
mkdir -p "$SCRATCH/iron_runs/RUN_NAME/pseudo"
cd "$SCRATCH/iron_runs/RUN_NAME"

cp /path/to/calculation.in input.in
cp /path/to/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF pseudo/
cp /path/to/IronCoreMD/hpc/lonestar/qe_pw_template.sbatch run.sbatch
```

Edit `YOUR_TACC_ALLOCATION`, the queue, wall time, node count, task count, module names, and `QE` path. Then:

```bash
sbatch run.sbatch
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,ExitCode,MaxRSS
tail -f qe.out
```

Use an interactive allocation for a small launch test; do not run a multi-rank QE calculation directly on a login node:

```bash
idev -N 1 -n 4 -m 30
ibrun "$QE" -in small_test.in > small_test.out
```

## Stage 1 to stage 2

The stage transition is the same as on Jakar: stage 2 must use the same `prefix` and must see the untouched stage-1 `outdir`. If the run directory is moved, move the entire `outdir` and verify permissions before submission. A true interrupted-MD continuation requires `restart_mode='restart'` and the saved ionic state; it must not be replaced by the original positions and velocities.

## Production record

Copy each actual submitted file into `hpc/lonestar/jobs/`. Record the modules, allocation, queue, QE path, source version, and build date in the script or a neighboring README.

## Official guide

See the [Lonestar6 user guide](https://docs.tacc.utexas.edu/hpc/lonestar6/) for current queues, node architecture, Slurm limits, and launcher guidance.
