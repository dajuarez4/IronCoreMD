#!/usr/bin/env bash
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"

vector_mean() {
    awk -v key="$2" '$0 ~ "^[[:space:]]*" key "[[:space:]]*:" {
        for(i=1;i<=NF;i++) if($i==":") {
            n++; slot=(n-1)%8; value[slot]=sqrt($(i+1)^2+$(i+2)^2+$(i+3)^2); break
        }
    } END {
        if(n>=8){sum=0; for(i=0;i<8;i++)sum+=value[i]; printf "%.6f",sum/8}
    }' "$1"
}

printf '%-30s %-7s %-18s %-8s %-8s %-10s %-15s %-12s %-12s %-11s\n' \
    case lambda status md_steps scf_iter max_cycle last_accuracy temp_K target_print mean_local
while IFS=, read -r index case_name lattice temperature target seed l1 l2 l3 requested_steps velocity_seed requirement; do
    [[ "$index" == serial_index ]] && continue
    folder="cases/$case_name"; md="$folder/md_10steps.out"
    status=NOT_STARTED; steps=-; total_iter=-; max_cycle=-; accuracy=-; temp=-; printed=-; moment=-
    latest=""
    for stem in stage1_4e-3 stage2_1e-3 stage3_3e-4; do
        [[ -f "$folder/$stem.out" ]] && latest="$folder/$stem.out"
    done
    if [[ -f "$md" ]]; then
        latest="$md"; steps=$(grep -c 'Entering Dynamics' "$md" || true)
        total_iter=$(grep -c 'iteration #' "$md" || true)
        max_cycle=$(awk '/iteration #/{n++} /convergence has been achieved|convergence NOT achieved/{if(n>m)m=n;n=0} END{print m+0}' "$md")
        if grep -q 'JOB DONE' "$md" && [[ "$steps" -eq "$requested_steps" ]]; then status=DONE
        elif grep -Eq 'Error in routine|convergence NOT achieved|Cholesky|cdiaghg' "$md"; then status=FAILED_MD
        else status=RUNNING_MD; fi
        value=$(grep 'temperature[[:space:]]*=' "$md" | tail -1 | awk '{print $(NF-1)}'); [[ -n "$value" ]] && temp=$value
    elif [[ -n "$latest" ]]; then
        if grep -Eq 'Error in routine|convergence NOT achieved|Cholesky|cdiaghg' "$latest"; then status=FAILED_SCF
        else status=PRECONDITIONING; fi
    fi
    if [[ -n "$latest" ]]; then
        value=$(grep 'estimated scf accuracy' "$latest" | tail -1 | awk '{print $(NF-1)}'); [[ -n "$value" ]] && accuracy=$value
        value=$(vector_mean "$latest" 'constrained moment'); [[ -n "$value" ]] && printed=$value
        value=$(vector_mean "$latest" 'magnetization'); [[ -n "$value" ]] && moment=$value
    fi
    printf '%-30s %-7s %-18s %-8s %-8s %-10s %-15s %-12s %-12s %-11s\n' \
        "$case_name" "$l3" "$status" "$steps" "$total_iter" "$max_cycle" "$accuracy" "$temp" "$printed" "$moment"
done < case_manifest.csv

if [[ -f round9_serial_status.tsv ]]; then
    echo
    column -t -s $'\t' round9_serial_status.tsv 2>/dev/null || cat round9_serial_status.tsv
fi
