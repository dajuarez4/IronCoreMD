#!/bin/bash
#SBATCH --job-name=r20_dlm_pilot
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=%j.log

set -e
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_SHM_CMA=0

cd "${SLURM_SUBMIT_DIR}"
ulimit -s unlimited
module purge
source /opt/intel/oneapi/setvars.sh --force

QE=/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x
CASE=$(head -n 1 cases.txt)
cd "cases/$CASE"

mpirun -np "${SLURM_NTASKS}" "$QE" -in stage1_4e-3.in > stage1_4e-3.out
mpirun -np "${SLURM_NTASKS}" "$QE" -in stage2_1e-3.in > stage2_1e-3.out
mpirun -np "${SLURM_NTASKS}" "$QE" -in stage3_3e-4.in > stage3_3e-4.out
mpirun -np "${SLURM_NTASKS}" "$QE" -in md_100steps.in > md_100steps.out
