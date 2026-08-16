#!/bin/bash
#SBATCH --job-name=build_qe75_fix
#SBATCH --partition=vm-small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=logs/build_qe75_fix.%j.out
#SBATCH --error=logs/build_qe75_fix.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch}"
module purge
module load intel/19.1.1 impi/19.0.9
export OMP_NUM_THREADS=1
SOURCE_DIR=${SCRATCH:?}/src/q-e-qe-7.5-atomic-fixed
INSTALL_DIR=${SCRATCH}/apps/qe-7.5-atomic-fixed
[[ -f "$SOURCE_DIR/ATOMIC_CONSTRAINT_PATCH_INFO.txt" ]] || { echo "Patched source not prepared." >&2; exit 2; }
grep -q 'mcons(1,nt) = zv(nt) \* starting_magnetization' "$SOURCE_DIR/PW/src/input.f90" || { echo "Source patch missing." >&2; exit 3; }
cd "$SOURCE_DIR"
if [[ ! -f make.inc ]]; then
    ./configure MPIF90=mpif90 F90=ifort CC=icc
fi
make -j "${SLURM_CPUS_PER_TASK:-16}" pw
[[ -x bin/pw.x ]] || { echo "Build finished without bin/pw.x" >&2; exit 4; }
mkdir -p "$INSTALL_DIR/bin"
cp -p bin/pw.x "$INSTALL_DIR/bin/pw.x"
cp -p ATOMIC_CONSTRAINT_PATCH_INFO.txt "$INSTALL_DIR/"
sha256sum "$INSTALL_DIR/bin/pw.x" > "$INSTALL_DIR/pw.x.sha256"
echo "Installed patched executable: $INSTALL_DIR/bin/pw.x"
