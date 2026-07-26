#!/usr/bin/env bash
#SBATCH --job-name=qe_pw
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

module purge
source /opt/intel/oneapi/setvars.sh --force

ulimit -s unlimited
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

QE=${QE:-/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x}
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
command -v mpirun
command -v mpiifx

mpirun -np "$SLURM_NTASKS" "$QE" -in "$INPUT" > "$OUTPUT"

echo "End: $(date --iso-8601=seconds)"

