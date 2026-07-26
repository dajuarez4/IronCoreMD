#!/usr/bin/env bash
#SBATCH --job-name=qe_pw
#SBATCH --account=YOUR_TACC_ALLOCATION
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

module purge
module load COMPILER_MODULE
module load MPI_MODULE

ulimit -s unlimited
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

QE=${QE:-"$SCRATCH/software/qe-7.5/build/bin/pw.x"}
INPUT=${INPUT:-input.in}
OUTPUT=${OUTPUT:-qe.out}

if [[ ! -x "$QE" ]]; then
    echo "QE executable is missing or not executable: $QE" >&2
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "QE input is missing: $INPUT" >&2
    exit 1
fi

echo "Start: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Directory: $PWD"
echo "Tasks: $SLURM_NTASKS"
echo "QE: $QE"
echo "Input: $INPUT"
module list 2>&1
command -v ibrun

ibrun "$QE" -in "$INPUT" > "$OUTPUT"

echo "End: $(date --iso-8601=seconds)"

