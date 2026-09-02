#!/bin/bash
#SBATCH --job-name=build_qe75_fix128
#SBATCH --account=DMR26002
#SBATCH --partition=vm-small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=build_qe75_fix128.%j.out
#SBATCH --error=build_qe75_fix128.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit this file with sbatch}"
module purge

# Override at submission time if TACC has changed the defaults, for example:
# sbatch --export=ALL,COMPILER_MODULE=intel/19.1.1,MPI_MODULE=impi/19.0.9 ...
COMPILER_MODULE=${COMPILER_MODULE:-intel}
MPI_MODULE=${MPI_MODULE:-impi}
module load "$COMPILER_MODULE" "$MPI_MODULE"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

SOURCE_DIR=${SOURCE_DIR:-${SCRATCH:?}/src/q-e-qe-7.5-atomic-fixed-ntyp128}
INSTALL_DIR=${INSTALL_DIR:-${SCRATCH}/apps/qe-7.5-atomic-fixed-ntyp128}
INFO_FILE=ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt

[[ -f "$SOURCE_DIR/$INFO_FILE" ]] || {
    echo "Prepared-source provenance is missing: $SOURCE_DIR/$INFO_FILE" >&2
    exit 2
}
[[ ! -e "$INSTALL_DIR" ]] || {
    echo "$INSTALL_DIR already exists; refusing to overwrite it." >&2
    exit 2
}
grep -Eq 'ntypx[[:space:]]*=[[:space:]]*128' "$SOURCE_DIR/Modules/parameters.f90" || {
    echo "ntypx=128 is missing from the source." >&2
    exit 3
}
grep -q 'mcons(1,nt) = zv(nt) \* starting_magnetization' "$SOURCE_DIR/PW/src/input.f90" || {
    echo "The zv(nt) atomic-constraint patch is missing." >&2
    exit 3
}

command -v ibrun
command -v mpif90
mpif90 -show

if command -v ifort >/dev/null 2>&1; then
    SERIAL_F90=ifort
elif command -v ifx >/dev/null 2>&1; then
    SERIAL_F90=ifx
else
    echo "Neither ifort nor ifx is available after loading $COMPILER_MODULE." >&2
    exit 4
fi

if command -v icc >/dev/null 2>&1; then
    SERIAL_CC=icc
elif command -v icx >/dev/null 2>&1; then
    SERIAL_CC=icx
else
    echo "Neither icc nor icx is available after loading $COMPILER_MODULE." >&2
    exit 4
fi

cd "$SOURCE_DIR"
./configure MPIF90=mpif90 F90="$SERIAL_F90" CC="$SERIAL_CC"
grep -q '__MPI' make.inc || {
    echo "QE configure did not enable MPI." >&2
    exit 5
}
make -j "${SLURM_CPUS_PER_TASK:-16}" pw
[[ -x bin/pw.x ]] || { echo "Build completed without bin/pw.x." >&2; exit 6; }

mkdir -p "$INSTALL_DIR/bin"
cp -Lp bin/pw.x "$INSTALL_DIR/bin/pw.x"
cp -p "$INFO_FILE" "$INSTALL_DIR/$INFO_FILE"
cp -p make.inc "$INSTALL_DIR/make.inc"

{
    echo "Built UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "Build host: $(hostname)"
    echo "Compiler module requested: $COMPILER_MODULE"
    echo "MPI module requested: $MPI_MODULE"
    echo "Fortran compiler: $SERIAL_F90"
    echo "C compiler: $SERIAL_CC"
    echo "MPI wrapper: $(command -v mpif90)"
    echo "MPI wrapper expansion:"
    mpif90 -show
    echo "Loaded modules:"
    module list 2>&1
    echo "Executable SHA256:"
    sha256sum "$INSTALL_DIR/bin/pw.x"
} > "$INSTALL_DIR/BUILD_ENVIRONMENT.txt"

cat "$INSTALL_DIR/BUILD_ENVIRONMENT.txt" >> "$INSTALL_DIR/$INFO_FILE"
sha256sum "$INSTALL_DIR/bin/pw.x" > "$INSTALL_DIR/pw.x.sha256"
ldd "$INSTALL_DIR/bin/pw.x" > "$INSTALL_DIR/pw.x.ldd.txt"
if grep -q 'not found' "$INSTALL_DIR/pw.x.ldd.txt"; then
    echo "The executable has unresolved shared libraries." >&2
    cat "$INSTALL_DIR/pw.x.ldd.txt" >&2
    exit 7
fi

echo "Installed private QE: $INSTALL_DIR/bin/pw.x"
cat "$INSTALL_DIR/pw.x.sha256"
echo "Run a short SCF and verify the printed 2.0-mu_B constraint before production."
