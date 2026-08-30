#!/bin/bash
#SBATCH --job-name=test_qe75_n128
#SBATCH --account=DMR26002
#SBATCH --partition=vm-small
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
#SBATCH --output=logs/test_qe75_n128.%j.out
#SBATCH --error=logs/test_qe75_n128.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"
module purge
module load intel/19.1.1 impi/19.0.9
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

INSTALL=${SCRATCH:?}/apps/qe-7.5-atomic-fixed-ntyp128
PWX=$INSTALL/bin/pw.x
INFO=$INSTALL/ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt
[[ -x "$PWX" && -s "$INFO" ]] || { echo "Missing executable/provenance: $INSTALL" >&2; exit 2; }
grep -q 'ntypx increased from 10 to 128' "$INFO" || {
    echo "Installed provenance does not confirm ntypx=128." >&2; exit 2;
}
rm -f test_ntyp16.out
rm -rf test_ntyp16_tmp
ibrun "$PWX" -in test_ntyp16.in > test_ntyp16.out

# QE 7.5 prints ntypx with an I2 field. Values above 99 therefore appear as
# "**"; exact 128 is established by the installed patch provenance above.
grep -Eq 'Max number of different atomic species.*=[[:space:]]*(128|\*\*)' test_ntyp16.out || {
    echo "Executable did not report the expected ntypx field." >&2; exit 3;
}
grep -Eq 'number of atomic types[[:space:]]*=[[:space:]]*16' test_ntyp16.out || {
    echo "Executable did not parse 16 atomic species." >&2; exit 3;
}
grep -q 'JOB DONE' test_ntyp16.out || { echo "Smoke test did not finish." >&2; exit 3; }
grep -Eq 'Error in routine|Cholesky|cdiaghg' test_ntyp16.out && {
    echo "Smoke test contains a QE error." >&2; exit 3;
}

echo "PASS: patched pw.x provenance confirms ntypx=128 and QE parsed all 16 species."
