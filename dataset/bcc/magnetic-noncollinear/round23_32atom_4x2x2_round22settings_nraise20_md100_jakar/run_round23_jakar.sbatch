#!/bin/bash
#SBATCH --job-name=r23_bcc32
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH --output=logs/%j.log

set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_SHM_CMA=0

cd "${SLURM_SUBMIT_DIR}"
ulimit -s unlimited
module purge
set +u
source /opt/intel/oneapi/setvars.sh --force
set -u

QE=/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x
CASE=$(sed -n '1p' cases.txt)
cd "cases/$CASE"

for STEM in stage1_4e-3 stage2_1e-3 stage3_1e-4 md_100steps; do
    mpirun -np "${SLURM_NTASKS}" "$QE" -in "$STEM.in" > "$STEM.out"
done
