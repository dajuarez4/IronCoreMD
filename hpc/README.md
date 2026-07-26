# Quantum ESPRESSO on Jakar and Lonestar6

This directory documents the reproducible path from a prepared Quantum ESPRESSO input to a submitted and archived HPC calculation.

## Directory map

```text
hpc/
├── README.md
├── install_qe_from_source.sh
├── jakar/
│   ├── README.md
│   ├── qe_pw_template.sbatch
│   └── jobs/
└── lonestar/
    ├── README.md
    ├── qe_pw_template.sbatch
    └── jobs/
```

Place copies of the SBATCH files actually used for production in the corresponding `jobs/` directory. Keep the reusable templates generic and record run-specific changes in the production copy.

## End-to-end checklist

1. Prepare the QE input and pseudopotential locally.
2. Check that `prefix`, `outdir`, `pseudo_dir`, `nat`, `ntyp`, and `K_POINTS` are consistent.
3. Copy the input, pseudopotential, and SBATCH file to the cluster.
4. Load the same compiler and MPI environment used to build QE.
5. Run a short SCF test before a long MD calculation.
6. Submit with `sbatch`, monitor with `squeue`, and inspect the QE output.
7. For a staged noncollinear run, complete stage 1 before stage 2 and preserve the entire `outdir`.
8. Copy the final input, output, scheduler logs, and provenance file back to permanent storage.

The Slurm wall time and QE `max_seconds` are separate limits. When using `max_seconds`, set it several minutes shorter than `#SBATCH --time` so QE has time to write restart data before Slurm terminates the job.

An SCF-to-MD stage transition is not the same as continuing an interrupted MD trajectory. Stage 2 may start a new MD calculation from stage-1 charge density and wavefunctions. A true continuation must retain the complete saved state, use the same `prefix` and `outdir`, and use `restart_mode='restart'` without replacing the saved trajectory state with fresh positions or velocities.

## Transfer a run directory

From the local machine:

```bash
rsync -avP prepared_qe_inputs/noncollinear_3x3x3_4000K/ \
  USER@CLUSTER:REMOTE_RUN_DIRECTORY/
```

From the cluster after completion:

```bash
rsync -avP USER@CLUSTER:REMOTE_RUN_DIRECTORY/ local_archive/
```

Replace `USER`, `CLUSTER`, and `REMOTE_RUN_DIRECTORY` with the values for the target system. Do not place large QE scratch directories in Git.

## Install QE from source

Download the QE source archive from the official release page, transfer it to the cluster, load the desired compiler/MPI environment, and run:

```bash
QE_ARCHIVE=/path/to/qe-7.5-ReleasePack.tar.gz \
QE_INSTALL_ROOT="$SCRATCH/software" \
BUILD_JOBS=8 \
bash hpc/install_qe_from_source.sh
```

The script uses the preferred out-of-source sequence:

```bash
mkdir build
cd build
../configure
make -j 8 all
```

It does not load cluster modules because the exact module stack is system-specific. Load that environment first and use the same stack when running QE.

## Record provenance

Run these commands inside every production directory before submission:

```bash
date > provenance.txt
hostname >> provenance.txt
module list 2>> provenance.txt
command -v pw.x >> provenance.txt 2>&1 || true
sha256sum *.in *.UPF >> provenance.txt 2>/dev/null || true
git -C /path/to/IronCoreMD rev-parse HEAD >> provenance.txt 2>/dev/null || true
```

On macOS, use `shasum -a 256` instead of `sha256sum`.

## References

- [Quantum ESPRESSO 7.5 user guide](https://www.quantum-espresso.org/Doc/user_guide/)
- [Quantum ESPRESSO source-build instructions](https://www.quantum-espresso.org/Doc/user_guide/node11.html)
- [Lonestar6 user guide](https://docs.tacc.utexas.edu/hpc/lonestar6/)
