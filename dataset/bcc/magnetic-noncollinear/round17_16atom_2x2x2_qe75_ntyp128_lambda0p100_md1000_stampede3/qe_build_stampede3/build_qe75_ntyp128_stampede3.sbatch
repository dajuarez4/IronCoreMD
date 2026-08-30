#!/bin/bash
#SBATCH --job-name=build_qe75_n128
#SBATCH --partition=skx-dev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=logs/build_qe75_n128.%j.out
#SBATCH --error=logs/build_qe75_n128.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"
module purge
module load intel/24.0 impi/21.11
export OMP_NUM_THREADS=1

SOURCE_DIR=${SCRATCH:?}/src/q-e-qe-7.5-atomic-fixed-ntyp128
INSTALL_DIR=${SCRATCH}/apps/qe-7.5-atomic-fixed-ntyp128
INFO=ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt

[[ -f "$SOURCE_DIR/$INFO" ]] || { echo "Patched source not prepared." >&2; exit 2; }
[[ ! -e "$INSTALL_DIR" ]] || {
    echo "$INSTALL_DIR already exists; refusing to overwrite a private installation." >&2
    exit 2
}
grep -Eq 'INTEGER,[[:space:]]*PARAMETER[[:space:]]*::[[:space:]]*ntypx[[:space:]]*=[[:space:]]*128' \
    "$SOURCE_DIR/Modules/parameters.f90" || { echo "ntypx=128 source patch missing." >&2; exit 3; }
grep -q 'mcons(1,nt) = zv(nt) \* starting_magnetization' "$SOURCE_DIR/PW/src/input.f90" || {
    echo "Atomic-constraint source patch missing." >&2; exit 3;
}
command -v mpif90
mpif90 -show

cd "$SOURCE_DIR"
./configure MPIF90=mpif90 CC=mpicc
grep -q '__MPI' make.inc || { echo "configure did not enable MPI." >&2; exit 4; }
make -j "${SLURM_CPUS_PER_TASK:-16}" pw
[[ -x bin/pw.x ]] || { echo "Build finished without bin/pw.x." >&2; exit 5; }

mkdir -p "$INSTALL_DIR/bin"
cp -p bin/pw.x "$INSTALL_DIR/bin/pw.x"
cp -p "$INFO" "$INSTALL_DIR/$INFO"
cp -p make.inc "$INSTALL_DIR/make.inc"
sha256sum "$INSTALL_DIR/bin/pw.x" > "$INSTALL_DIR/pw.x.sha256"
ldd "$INSTALL_DIR/bin/pw.x" > "$INSTALL_DIR/pw.x.ldd.txt"

echo "Installed private QE: $INSTALL_DIR/bin/pw.x"
cat "$INSTALL_DIR/pw.x.sha256"
