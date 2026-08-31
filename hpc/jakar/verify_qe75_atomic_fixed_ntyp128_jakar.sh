#!/usr/bin/env bash
# Verify build provenance/libraries and, optionally, constrained moments in a QE output.
# Usage:
#   bash verify_qe75_atomic_fixed_ntyp128_jakar.sh [install_dir] [qe_output]
# Environment:
#   EXPECTED_NORM=2.0 TOLERANCE=0.0001 bash verify_... [install_dir] output.out
set -euo pipefail

INSTALL_DIR=${1:-/gpfs/scratch/${USER}/apps/qe-7.5-atomic-fixed-ntyp128}
QE_OUTPUT=${2:-}
EXPECTED_NORM=${EXPECTED_NORM:-2.0}
TOLERANCE=${TOLERANCE:-0.0001}
PW="$INSTALL_DIR/bin/pw.x"
INFO="$INSTALL_DIR/ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt"

[[ -x "$PW" ]] || { echo "Missing executable: $PW" >&2; exit 2; }
[[ -s "$INFO" ]] || { echo "Missing provenance: $INFO" >&2; exit 2; }
grep -q 'ntypx increased from 10 to 128' "$INFO" || {
    echo "Provenance does not confirm ntypx=128." >&2
    exit 3
}
grep -q 'atomic mcons target multiplied by zv(nt)' "$INFO" || {
    echo "Provenance does not confirm the atomic-constraint patch." >&2
    exit 3
}

module purge
set +u
source /opt/intel/oneapi/setvars.sh --force
set -u

echo "===== EXECUTABLE ====="
readlink -f "$PW"
sha256sum "$PW"
if ldd "$PW" | grep 'not found'; then
    echo "FAIL: unresolved shared libraries." >&2
    exit 4
else
    echo "PASS: all shared libraries found."
fi

if [[ -z "$QE_OUTPUT" ]]; then
    echo "Build checks passed. Provide a QE output to validate the runtime target norm."
    exit 0
fi
[[ -s "$QE_OUTPUT" ]] || { echo "Missing QE output: $QE_OUTPUT" >&2; exit 2; }

echo "===== RUNTIME CONSTRAINED-MOMENT NORMS ====="
awk -v expected="$EXPECTED_NORM" -v tolerance="$TOLERANCE" '
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
    if (n == 0) {
        print "FAIL: no constrained moment lines found." > "/dev/stderr"
        exit 5
    }
    mean=total/n
    printf "count=%d min=%.9f mean=%.9f max=%.9f expected=%.9f\n", \
           n, minimum, mean, maximum, expected
    if (minimum < expected-tolerance || maximum > expected+tolerance) {
        print "FAIL: constrained-moment norm is outside tolerance." > "/dev/stderr"
        exit 6
    }
    print "PASS: runtime constrained-moment norm matches the requested target."
}' "$QE_OUTPUT"
