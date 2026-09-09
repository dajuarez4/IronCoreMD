#!/bin/bash
#SBATCH --job-name=r21_scf_test
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --array=0-2%2
#SBATCH --output=logs/%A_%a.log

set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_SHM_CMA=0

cd "${SLURM_SUBMIT_DIR}"
ulimit -s unlimited
module purge
# Intel's setvars.sh references optional variables before defining them and is
# not compatible with `set -u`. Disable nounset only while sourcing it.
set +u
source /opt/intel/oneapi/setvars.sh --force
set -u

QE=/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x
CASE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" cases.txt)
cd "cases/$CASE"

for STEM in stage1_4e-3 stage2_1e-3 stage3_3e-4 md_100steps; do
    mpirun -np "${SLURM_NTASKS}" "$QE" -in "$STEM.in" > "$STEM.out"
done
