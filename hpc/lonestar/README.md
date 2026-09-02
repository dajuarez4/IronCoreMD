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

## Private QE 7.5 build: atomic constraint fix and 128 species

The constrained-DLM calculations need the same two private changes validated
on Jakar:

- `Modules/parameters.f90`: `ntypx=128`;
- `PW/src/input.f90`:
  `mcons = zv(nt) * starting_magnetization(nt) * direction`.

For the Fe PAW pseudopotential used by this project, `zv=16`, so
`starting_magnetization=0.125` must print a constrained-vector norm of
`2.000000 mu_B`. A successful compilation is not sufficient without this
runtime test.

The reusable Lonestar6 files are:

```text
qe75_atomic_constraint_ntyp128.patch
prepare_qe75_atomic_fixed_ntyp128_lonestar6.sh
build_qe75_atomic_fixed_ntyp128_lonestar6.sbatch
verify_qe75_atomic_fixed_ntyp128_lonestar6.sh
test_qe75_atomic_fixed_ntyp128_lonestar6.sbatch
test_ntyp16.in
SCALING_16_TO_128.md
```

Before building, inspect the current TACC environment:

```bash
qlimits
module spider intel
module spider impi
```

Prepare a clean source tree:

```bash
cd /path/to/IronCoreMD/hpc/lonestar

bash prepare_qe75_atomic_fixed_ntyp128_lonestar6.sh \
  "$SCRATCH/software/q-e-qe-7.5.tar.gz" \
  "$SCRATCH/src"
```

The preparation script refuses to overwrite an existing tree. For a later
rebuild, choose new source and installation names explicitly:

```bash
TARGET_DIR="$SCRATCH/src/q-e-qe-7.5-atomic-fixed-ntyp128-v2" \
bash prepare_qe75_atomic_fixed_ntyp128_lonestar6.sh \
  "$SCRATCH/software/q-e-qe-7.5.tar.gz" \
  "$SCRATCH/src"

sbatch \
  --export=ALL,SOURCE_DIR="$SCRATCH/src/q-e-qe-7.5-atomic-fixed-ntyp128-v2",INSTALL_DIR="$SCRATCH/apps/qe-7.5-atomic-fixed-ntyp128-v2" \
  build_qe75_atomic_fixed_ntyp128_lonestar6.sbatch
```

Compile in `vm-small`:

```bash
sbatch build_qe75_atomic_fixed_ntyp128_lonestar6.sbatch
```

The script defaults to the unversioned `intel` and `impi` modules, records the
resolved build environment, and installs:

```text
$SCRATCH/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x
```

If the desired modules are versioned, export them explicitly:

```bash
sbatch --export=ALL,COMPILER_MODULE=intel/19.1.1,MPI_MODULE=impi/19.0.9 \
  build_qe75_atomic_fixed_ntyp128_lonestar6.sbatch
```

After the build, verify provenance and libraries using the same module stack:

```bash
COMPILER_MODULE=intel MPI_MODULE=impi \
bash verify_qe75_atomic_fixed_ntyp128_lonestar6.sh \
  "$SCRATCH/apps/qe-7.5-atomic-fixed-ntyp128"
```

Copy the trusted Fe PAW pseudopotential beside `test_ntyp16.in`, then submit
the included smoke test:

```bash
cp /path/to/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF .
sbatch test_qe75_atomic_fixed_ntyp128_lonestar6.sbatch
```

The final line of the Slurm output must be:

```text
PASS: QE parsed 16 species and printed 2.0-mu_B atomic constraints.
```

Run a short 16-species SCF through `ibrun`, then validate its output:

```bash
COMPILER_MODULE=intel MPI_MODULE=impi EXPECTED_NORM=2.0 \
bash verify_qe75_atomic_fixed_ntyp128_lonestar6.sh \
  "$SCRATCH/apps/qe-7.5-atomic-fixed-ntyp128" \
  /absolute/path/to/stage1_4e-3.out
```

The verification must report `PASS` and norms of `2.000000 mu_B`. Record the
new executable SHA256; it is expected to differ from the Jakar hash because the
compiler, MPI, and linked libraries differ.

See [`SCALING_16_TO_128.md`](SCALING_16_TO_128.md) before submitting larger
cells.

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
