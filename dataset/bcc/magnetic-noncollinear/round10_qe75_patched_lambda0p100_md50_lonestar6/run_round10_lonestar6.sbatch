#!/bin/bash
#SBATCH --job-name=bcc8_r10_md50
#SBATCH --partition=vm-small
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=48:00:00
#SBATCH --output=logs/ls6_round10_md50.%j.out
#SBATCH --error=logs/ls6_round10_md50.%j.err

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"
module purge
module load intel/19.1.1 impi/19.0.9
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

PWX=${SCRATCH:?}/apps/qe-7.5-atomic-fixed/bin/pw.x
INFO=${SCRATCH}/apps/qe-7.5-atomic-fixed/ATOMIC_CONSTRAINT_PATCH_INFO.txt
[[ -x "$PWX" && -s "$INFO" ]] || {
    echo "Missing verified private patched QE 7.5 executable." >&2
    exit 2
}
command -v ibrun >/dev/null || { echo "ibrun unavailable" >&2; exit 2; }
echo "Using private patched QE: $PWX"
sha256sum "$PWX"

mapfile -t cases < cases.txt
stages=(stage1_4e-3 stage2_1e-3 stage3_3e-4)
printf 'case\tstatus\tfinished_at\n' > round10_status.tsv
failed=0

for case_name in "${cases[@]}"; do
    folder="cases/$case_name"
    echo "=== $case_name starting $(date) ==="
    if [[ -e "$folder/qe_tmp" ]] || compgen -G "$folder/*.out" >/dev/null; then
        echo "$case_name: existing state/output; skipped to avoid mixing attempts." >&2
        printf '%s\tSKIPPED_EXISTING\t%s\n' "$case_name" "$(date -Is)" >> round10_status.tsv
        failed=$((failed + 1))
        continue
    fi

    bad=0
    for stage in "${stages[@]}"; do
        echo "$case_name: $stage starting $(date)"
        (cd "$folder" && ibrun "$PWX" -in "$stage.in" > "$stage.out")
        rc=$?
        if [[ "$rc" -ne 0 ]] || ! grep -q 'convergence has been achieved' "$folder/$stage.out"; then
            printf '%s\tFAILED_%s\t%s\n' "$case_name" "$stage" "$(date -Is)" >> round10_status.tsv
            failed=$((failed + 1)); bad=1; break
        fi
    done
    [[ "$bad" -eq 0 ]] || continue

    echo "$case_name: 50-step MD starting $(date)"
    (cd "$folder" && ibrun "$PWX" -in md_50steps.in > md_50steps.out)
    rc=$?
    steps=$(grep -c 'Entering Dynamics' "$folder/md_50steps.out" || true)
    if [[ "$rc" -ne 0 ]] || ! grep -q 'JOB DONE' "$folder/md_50steps.out" || [[ "$steps" -ne 50 ]]; then
        printf '%s\tFAILED_OR_PARTIAL_MD\t%s\n' "$case_name" "$(date -Is)" >> round10_status.tsv
        failed=$((failed + 1))
        continue
    fi
    printf '%s\tDONE\t%s\n' "$case_name" "$(date -Is)" >> round10_status.tsv
done

column -t -s $'\t' round10_status.tsv 2>/dev/null || cat round10_status.tsv
[[ "$failed" -eq 0 ]] || exit 1
