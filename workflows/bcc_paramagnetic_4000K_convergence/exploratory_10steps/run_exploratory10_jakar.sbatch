#!/bin/bash
#SBATCH --job-name=bcc_exp10
#SBATCH --account=jakar_partner_physj
#SBATCH --qos=jakar_medium_physj
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
#SBATCH --output=%j.log

set -eo pipefail

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_SHM_CMA=0

cd "${SLURM_SUBMIT_DIR}"

prefix=Fe_bcc2_pm_round4
input=md_exploratory_4000K_10steps.in
output=md_exploratory_4000K_10steps.out

if [[ ! -d "seed_state/${prefix}.save" ]]; then
    echo "Missing seed_state/${prefix}.save" >&2
    echo "Copy the untouched lambda=0.020 preconditioned state into seed_state." >&2
    exit 2
fi

shopt -s nullglob
seed_wfc=(seed_state/${prefix}.wfc*)
if (( ${#seed_wfc[@]} != 40 )); then
    echo "Expected 40 seed wavefunction files; found ${#seed_wfc[@]}." >&2
    exit 3
fi

if [[ -e md_state || -e "${output}" ]]; then
    echo "md_state or ${output} already exists; refusing to mix attempts." >&2
    exit 4
fi

cp -R seed_state md_state

module purge
source /opt/intel/oneapi/setvars.sh --force
set -u

QE=/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x

echo "Starting exploratory BCC 4000 K ten-step MD at $(date)"
mpirun -np "${SLURM_NTASKS}" "${QE}" -in "${input}" > "${output}"

if ! grep -q "JOB DONE\." "${output}"; then
    echo "Exploratory MD stopped before JOB DONE." >&2
    exit 5
fi

echo "Exploratory BCC MD completed at $(date)"
