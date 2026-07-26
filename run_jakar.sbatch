#!/bin/bash
#SBATCH --job-name=2.76_1700K
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --output=%j.log
#SBATCH --error=%j.err

pwd
hostname
date
ulimit -s unlimited

module purge
source /opt/intel/oneapi/setvars.sh --force

which mpirun
which mpiifx

QE=/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x
INPUT=fe_bcc_2x_2.75_non_col_0.5.in
OUTPUT=fe1.out

mpirun -np ${SLURM_NTASKS} ${QE} -in ${INPUT} > ${OUTPUT}