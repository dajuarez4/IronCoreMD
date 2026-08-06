#!/bin/bash
#SBATCH --job-name=fe128_fm_a278
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --output=%j.log

# Intel oneAPI's setvars.sh references optional variables that may be unset.
set -eo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p tmp

pwd
hostname
date
ulimit -s unlimited

source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \
  -in 01_scf.in > 01_scf.out

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \
  -in 02_bands.in > 02_bands.out

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/bands.x \
  -in 03_bands_up.in > 03_bands_up.out

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/bands.x \
  -in 04_bands_down.in > 04_bands_down.out

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \
  -in 05_nscf.in > 05_nscf.out

mpirun --mca pml ob1 --mca btl ^openib,uct \
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/projwfc.x \
  -in 06_projwfc.in > 06_projwfc.out

date
