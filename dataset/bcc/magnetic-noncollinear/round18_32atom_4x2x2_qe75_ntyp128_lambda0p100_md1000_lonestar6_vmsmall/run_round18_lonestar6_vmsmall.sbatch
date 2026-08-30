#!/bin/bash
#SBATCH --job-name=bcc32_r18
#SBATCH --account=DMR26002
#SBATCH --partition=vm-small
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --output=%j.log

set -e

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}"
module purge
module load intel/19.1.1 impi/19.0.9

QE=${SCRATCH}/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x
CASE=cases/a2p39_4x2x2_m2_32species_lambda0p100_md1000

cd "$CASE"

ibrun "$QE" -in stage1_4e-3.in > stage1_4e-3.out
ibrun "$QE" -in stage2_1e-3.in > stage2_1e-3.out
ibrun "$QE" -in stage3_3e-4.in > stage3_3e-4.out
ibrun "$QE" -in md_1000steps.in > md_1000steps.out
