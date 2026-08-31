# Jakar workflow

The settings below are based on the repository's current `run_jakar.sbatch`. Verify the account, QoS, partition, and task limits before each campaign because cluster policy can change.

## Current environment

The Round 20 constrained-DLM campaign and the private QE build use:

- account: `jakar_general`
- QoS: `jakar_medium_general`
- partition: `medium`
- one node with 16 MPI tasks for QE; one task with 16 build threads for compilation
- Intel oneAPI initialized with `/opt/intel/oneapi/setvars.sh --force`
- Intel MPI launched with `mpirun`
- existing QE binary: `/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x`

Older job templates in this repository use the `jakar_partner_physj` account
and `jakar_medium_physj` QoS. Verify current allocation policy before copying
either set of directives into a new campaign.

## Stock QE 7.5 installation

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

## Private QE 7.5 build for atomic constraints and 128 species

The noncollinear constrained-DLM calculations require two private source
changes that are not part of a stock QE 7.5 build:

1. `Modules/parameters.f90`: increase `ntypx` from 10 to 128. This permits the
   one-magnetic-species-per-atom construction used by the 16-, 32-, and
   48-atom calculations.
2. `PW/src/input.f90`: construct the atomic constraint as
   `zv(nt) * starting_magnetization(nt) * direction`. For the Fe PAW
   pseudopotential used here, `zv=16`; therefore
   `starting_magnetization=0.125` produces a `2.0 mu_B` target.

This is a project-specific patch. A source-code check alone is insufficient:
every new executable must also pass a runtime test in which QE prints the
expected constrained-moment norm.

The complete patch and scripts are stored beside this README:

```text
qe75_atomic_constraint_ntyp128.patch
prepare_qe75_atomic_fixed_ntyp128.sh
build_qe75_atomic_fixed_ntyp128_jakar.sbatch
verify_qe75_atomic_fixed_ntyp128_jakar.sh
```

### 1. Prepare a clean patched source tree

Start from an official QE 7.5 source archive. Do not prepare the new build by
copying compiled objects from another cluster.

```bash
cd /path/to/IronCoreMD/hpc/jakar

bash prepare_qe75_atomic_fixed_ntyp128.sh \
  /gpfs/scratch/$USER/software/q-e-qe-7.5.tar.gz \
  /gpfs/scratch/$USER/src
```

The command creates:

```text
/gpfs/scratch/$USER/src/q-e-qe-7.5-atomic-fixed-ntyp128
```

It refuses to overwrite an existing source directory and records the archive
and patch checksums in `ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt`.

### 2. Compile on a Jakar compute node

```bash
cd /path/to/IronCoreMD/hpc/jakar
sbatch build_qe75_atomic_fixed_ntyp128_jakar.sbatch
```

Monitor the build with:

```bash
squeue -u "$USER"
tail -f build_qe75_fix128.JOB_ID.out
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,ExitCode
```

The private installation is written to:

```text
/gpfs/scratch/$USER/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x
```

The build job refuses to overwrite an existing installation. To rebuild,
choose a new install directory by submitting with exported variables, for
example:

```bash
sbatch --export=ALL,INSTALL_DIR=/gpfs/scratch/$USER/apps/qe-7.5-atomic-fixed-ntyp128-v2 \
  build_qe75_atomic_fixed_ntyp128_jakar.sbatch
```

### 3. Verify libraries and runtime behavior

The verification script deliberately does not run `pw.x -help`: this QE build
interprets that invocation as a request to read an input file and waits on
standard input.

First verify provenance and shared libraries:

```bash
bash verify_qe75_atomic_fixed_ntyp128_jakar.sh \
  /gpfs/scratch/$USER/apps/qe-7.5-atomic-fixed-ntyp128
```

Then run a short SCF calculation with the intended pseudopotential and inspect
the resulting output:

```bash
EXPECTED_NORM=2.0 \
bash verify_qe75_atomic_fixed_ntyp128_jakar.sh \
  /gpfs/scratch/$USER/apps/qe-7.5-atomic-fixed-ntyp128 \
  /absolute/path/to/stage1_4e-3.out
```

The test passes only if every printed constrained vector has norm
`2.000000 +/- 0.0001 mu_B`. Production SBATCH scripts must point explicitly to
the verified private executable, and its SHA256 should be recorded with each
campaign.

### Repair record: 2026-08-31

The Jakar executable used initially by Round 19 and the old Round 20 replicas
had `ntypx=128`, but lacked the `zv(nt)` constraint scaling:

```text
unpatched SHA256 = dbbc2ea192a3fbddd19bff4c9d63131f57a918ceb46129047bb4364071cba842
runtime target   = 0.125000 mu_B
```

After applying the missing source change and rebuilding in place, the verified
Jakar executable was:

```text
path             = /gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x
patched SHA256   = 5aff2e0b9f36cd1dda210b5fb155b6573b6136612cf235ddc543477f41939abf
runtime target   = 2.000000 mu_B for all 16 Round 20 pilot species
```

Only outputs whose printed constrained vectors have the intended norm should
be combined into a dataset. A matching path without a matching hash/runtime
test is not sufficient provenance.

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
