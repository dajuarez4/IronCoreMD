#!/bin/bash
#SBATCH --job-name=build_qe75_fix128
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=build_qe75_fix128.%j.out

set -euo pipefail
module purge
set +u
source /opt/intel/oneapi/setvars.sh --force
set -u

SOURCE_DIR=${SOURCE_DIR:-/gpfs/scratch/${USER}/src/q-e-qe-7.5-atomic-fixed-ntyp128}
INSTALL_DIR=${INSTALL_DIR:-/gpfs/scratch/${USER}/apps/qe-7.5-atomic-fixed-ntyp128}
INFO_FILE=ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt

[[ -f "$SOURCE_DIR/$INFO_FILE" ]] || {
    echo "Prepared source provenance is missing: $SOURCE_DIR/$INFO_FILE" >&2
    exit 2
}
[[ ! -e "$INSTALL_DIR" ]] || {
    echo "$INSTALL_DIR already exists; refusing to overwrite it." >&2
    exit 2
}
grep -Eq 'ntypx[[:space:]]*=[[:space:]]*128' "$SOURCE_DIR/Modules/parameters.f90" || {
    echo "ntypx=128 is missing from the prepared source." >&2
    exit 3
}
grep -q 'mcons(1,nt) = zv(nt) \* starting_magnetization' "$SOURCE_DIR/PW/src/input.f90" || {
    echo "The zv(nt) atomic-constraint patch is missing." >&2
    exit 3
}

command -v mpirun
command -v mpif90 || command -v mpiifx || command -v mpiifort

cd "$SOURCE_DIR"
./configure
grep -q '__MPI' make.inc || {
    echo "QE configure did not enable MPI; inspect make.inc and the oneAPI environment." >&2
    exit 4
}
make -j "${SLURM_CPUS_PER_TASK:-16}" pw
[[ -x bin/pw.x ]] || { echo "Build finished without bin/pw.x." >&2; exit 5; }

mkdir -p "$INSTALL_DIR/bin"
cp -Lp bin/pw.x "$INSTALL_DIR/bin/pw.x"
cp -p "$INFO_FILE" "$INSTALL_DIR/$INFO_FILE"
cp -p make.inc "$INSTALL_DIR/make.inc"

{
    echo "Built UTC: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "Build host: $(hostname)"
    echo "Source directory: $SOURCE_DIR"
    echo "Install directory: $INSTALL_DIR"
    echo "Executable SHA256:"
    sha256sum "$INSTALL_DIR/bin/pw.x"
} >> "$INSTALL_DIR/$INFO_FILE"

sha256sum "$INSTALL_DIR/bin/pw.x" > "$INSTALL_DIR/pw.x.sha256"
ldd "$INSTALL_DIR/bin/pw.x" > "$INSTALL_DIR/pw.x.ldd.txt"
if grep -q 'not found' "$INSTALL_DIR/pw.x.ldd.txt"; then
    echo "The executable has unresolved shared libraries." >&2
    cat "$INSTALL_DIR/pw.x.ldd.txt" >&2
    exit 6
fi

echo "Installed private QE executable: $INSTALL_DIR/bin/pw.x"
cat "$INSTALL_DIR/pw.x.sha256"
echo "A calculation output containing constrained moments is still required for runtime validation."
