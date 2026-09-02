#!/bin/bash
#SBATCH --job-name=test_qe75_fix128
#SBATCH --account=DMR26002
#SBATCH --partition=vm-small
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
#SBATCH --output=test_qe75_fix128.%j.out
#SBATCH --error=test_qe75_fix128.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit this file with sbatch}"
module purge

COMPILER_MODULE=${COMPILER_MODULE:-intel}
MPI_MODULE=${MPI_MODULE:-impi}
module load "$COMPILER_MODULE" "$MPI_MODULE"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

INSTALL_DIR=${INSTALL_DIR:-${SCRATCH:?}/apps/qe-7.5-atomic-fixed-ntyp128}
PW="$INSTALL_DIR/bin/pw.x"
INFO="$INSTALL_DIR/ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt"
INPUT=${INPUT:-test_ntyp16.in}
OUTPUT=${OUTPUT:-test_ntyp16.out}
TMPDIR_QE=test_ntyp16_tmp
PSEUDO=Fe.pbe-spn-kjpaw_psl.1.0.0.UPF

[[ -x "$PW" ]] || { echo "Missing executable: $PW" >&2; exit 2; }
[[ -s "$INFO" ]] || { echo "Missing provenance: $INFO" >&2; exit 2; }
[[ -s "$INPUT" ]] || { echo "Missing input: $INPUT" >&2; exit 2; }
[[ -s "$PSEUDO" ]] || { echo "Copy $PSEUDO into $SLURM_SUBMIT_DIR first." >&2; exit 2; }
[[ ! -e "$OUTPUT" && ! -e "$TMPDIR_QE" ]] || {
    echo "$OUTPUT or $TMPDIR_QE already exists; archive it before rerunning." >&2
    exit 2
}

grep -q 'ntypx increased from 10 to 128' "$INFO" || {
    echo "Installed provenance does not confirm ntypx=128." >&2
    exit 3
}
grep -q 'atomic mcons target multiplied by zv(nt)' "$INFO" || {
    echo "Installed provenance does not confirm the atomic-constraint patch." >&2
    exit 3
}
if ldd "$PW" | grep 'not found'; then
    echo "Executable has unresolved shared libraries." >&2
    exit 3
fi

echo "QE executable: $(readlink -f "$PW")"
sha256sum "$PW"
module list 2>&1
command -v ibrun

ibrun "$PW" -in "$INPUT" > "$OUTPUT"

grep -Eq 'number of atomic types[[:space:]]*=[[:space:]]*16' "$OUTPUT" || {
    echo "QE did not parse all 16 species." >&2
    exit 4
}
grep -q 'JOB DONE' "$OUTPUT" || { echo "Smoke test did not finish." >&2; exit 4; }
if grep -Eq 'Error in routine|Cholesky|cdiaghg' "$OUTPUT"; then
    echo "Smoke test contains a QE error." >&2
    exit 4
fi

awk '
/constrained moment/ {
    x=$(NF-2); y=$(NF-1); z=$NF
    gsub(/[(),]/,"",x); gsub(/[(),]/,"",y); gsub(/[(),]/,"",z)
    norm=sqrt(x*x+y*y+z*z)
    if (n == 0 || norm < minimum) minimum=norm
    if (n == 0 || norm > maximum) maximum=norm
    total += norm
    n++
}
END {
    if (n < 16) {
        print "FAIL: fewer than 16 constrained moments were printed." > "/dev/stderr"
        exit 5
    }
    mean=total/n
    printf "constraint_count=%d min=%.9f mean=%.9f max=%.9f\n", \
           n, minimum, mean, maximum
    if (minimum < 1.9999 || maximum > 2.0001) {
        print "FAIL: runtime target is not 2.0 mu_B." > "/dev/stderr"
        exit 5
    }
}' "$OUTPUT"

echo "PASS: QE parsed 16 species and printed 2.0-mu_B atomic constraints."
