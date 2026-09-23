#!/bin/bash
#SBATCH --job-name=bands_r25
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH --output=bands-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
[[ -f checkpoint.json && -f 01_snapshot_scf.in && -f 02_bands.in ]]
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_SHM_CMA=0
ulimit -s unlimited
module purge
set +u
source /opt/intel/oneapi/setvars.sh --force
set -u
QE=/gpfs/scratch/dajuarez4/qe-7.5/bin
[[ -x "$QE/pw.x" ]] || { echo "Missing executable: $QE/pw.x" >&2; exit 1; }
# One stage per allocation; resubmit this same script after a completed stage.
exec 9>.bands.lock
flock -n 9 || exit 2
if [[ ! -f scf.complete ]]; then
    [[ ! -e 01_snapshot_scf.out ]]
    mpirun -np "$SLURM_NTASKS" "$QE/pw.x" -in 01_snapshot_scf.in > 01_snapshot_scf.out
    grep -q 'JOB DONE' 01_snapshot_scf.out
    grep -q 'convergence has been achieved' 01_snapshot_scf.out
    ! grep -q 'convergence NOT achieved' 01_snapshot_scf.out
    cp qe_tmp/Fe_bcc32_r25_conv1e4_beta1e2.save/data-file-schema.xml snapshot_scf.xml
    touch scf.complete
    echo 'Snapshot SCF complete. Submit this script again for the band path.'
    exit 0
fi
if [[ ! -f bands.complete ]]; then
    [[ ! -e 02_bands.out ]]
    mpirun -np "$SLURM_NTASKS" "$QE/pw.x" -in 02_bands.in > 02_bands.out
    grep -q 'JOB DONE' 02_bands.out
    grep -q 'End of band structure calculation' 02_bands.out
    touch bands.complete
fi
if [[ ! -f post.complete ]]; then
    [[ -x "$QE/bands.x" ]] || { echo "SCF and band path finished; postprocessing pending: $QE/bands.x is unavailable. Preserve qe_tmp." >&2; exit 3; }
    [[ ! -e 03_bands_post.out ]]
    mpirun -np "$SLURM_NTASKS" "$QE/bands.x" -in 03_bands_post.in > 03_bands_post.out
    grep -q 'JOB DONE' 03_bands_post.out
    touch post.complete
fi
echo 'Bands complete. Run python3 ../plot_bands.py from this directory.'
