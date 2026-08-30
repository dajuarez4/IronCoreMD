#!/bin/bash
#SBATCH --job-name=bcc16_r17_md1000
#SBATCH --partition=skx
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --output=logs/stampede3_round17_md1000.%j.out
#SBATCH --error=logs/stampede3_round17_md1000.%j.err

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit this script from the Round-17 directory}"
module purge
module load intel/24.0 impi/21.11
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

INSTALL=${SCRATCH:?}/apps/qe-7.5-atomic-fixed-ntyp128
PWX=$INSTALL/bin/pw.x
INFO=$INSTALL/ATOMIC_CONSTRAINT_NTYP128_PATCH_INFO.txt
[[ -x "$PWX" && -s "$INFO" && -s "$INSTALL/pw.x.sha256" ]] || {
    echo "Missing verified private QE 7.5 ntypx=128 installation: $INSTALL" >&2
    exit 2
}
grep -q 'ntypx increased from 10 to 128' "$INFO" || {
    echo "Installed QE provenance does not confirm ntypx=128." >&2
    exit 2
}
grep -q 'atomic mcons target scaled' "$INFO" || {
    echo "Installed QE provenance does not confirm the atomic-constraint fix." >&2
    exit 2
}
(cd "$INSTALL" && sha256sum -c pw.x.sha256) || exit 2
command -v ibrun >/dev/null || { echo "ibrun unavailable" >&2; exit 2; }

mapfile -t cases < cases.txt
stages=(stage1_4e-3 stage2_1e-3 stage3_3e-4)
printf 'case\tstatus\tfinished_at\n' > round17_status.tsv
failed=0

for case_name in "${cases[@]}"; do
    folder="cases/$case_name"
    echo "=== $case_name starting $(date -Is) ==="
    if [[ -e "$folder/qe_tmp" ]] || compgen -G "$folder/*.out" >/dev/null; then
        echo "$case_name: existing output/state; refusing to mix attempts." >&2
        printf '%s\tSKIPPED_EXISTING\t%s\n' "$case_name" "$(date -Is)" >> round17_status.tsv
        failed=$((failed + 1))
        continue
    fi

    bad=0
    for stage in "${stages[@]}"; do
        echo "$case_name: $stage starting $(date -Is)"
        (cd "$folder" && ibrun "$PWX" -in "$stage.in" > "$stage.out")
        rc=$?
        if [[ "$rc" -ne 0 ]] || ! grep -q 'convergence has been achieved' "$folder/$stage.out"; then
            printf '%s\tFAILED_%s\t%s\n' "$case_name" "$stage" "$(date -Is)" >> round17_status.tsv
            failed=$((failed + 1)); bad=1; break
        fi
    done
    [[ "$bad" -eq 0 ]] || continue

    echo "$case_name: 1000-step MD starting $(date -Is)"
    (cd "$folder" && ibrun "$PWX" -in md_1000steps.in > md_1000steps.out)
    rc=$?
    steps=$(grep -c 'Entering Dynamics' "$folder/md_1000steps.out" || true)
    if [[ "$rc" -ne 0 ]] || ! grep -q 'JOB DONE' "$folder/md_1000steps.out" || [[ "$steps" -ne 1000 ]]; then
        printf '%s\tFAILED_OR_PARTIAL_MD_%s_STEPS\t%s\n' "$case_name" "$steps" "$(date -Is)" >> round17_status.tsv
        failed=$((failed + 1))
        continue
    fi
    printf '%s\tDONE\t%s\n' "$case_name" "$(date -Is)" >> round17_status.tsv
done

column -t -s $'\t' round17_status.tsv 2>/dev/null || cat round17_status.tsv
[[ "$failed" -eq 0 ]] || exit 1
