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
        if(n>=32){sum=0; for(i=0;i<32;i++)sum+=value[i]; printf "%.6f",sum/32}
    }' "$1"
}

vector_tail_mean() {
    awk -v key="$2" -v window="$3" '$0 ~ "^[[:space:]]*" key "[[:space:]]*:" {
        for(i=1;i<=NF;i++) if($i==":") {
            slot=n%32; value[slot]=sqrt($(i+1)^2+$(i+2)^2+$(i+3)^2); n++;
            if(n%32==0){sum=0; for(j=0;j<32;j++)sum+=value[j]; report++; average[report]=sum/32}
            break
        }
    } END {
        first=(report-window+1>1 ? report-window+1 : 1); sum=0; count=0;
        for(i=first;i<=report;i++){sum+=average[i];count++}
        if(count) printf "%.6f",sum/count
    }' "$1"
}

printf '%-30s %-7s %-18s %-8s %-8s %-10s %-15s %-12s %-12s %-12s %-11s %-12s\n' \
    case lambda status md_steps scf_iter max_cycle last_accuracy temp_K temp_mean100 target_print mean_local moment_mean100
while IFS=, read -r index case_name lattice temperature target seed l1 l2 l3 requested_steps velocity_seed requirement; do
    [[ "$index" == serial_index ]] && continue
    folder="cases/$case_name"; md="$folder/md_400steps.out"
    status=NOT_STARTED; steps=-; total_iter=-; max_cycle=-; accuracy=-; temp=-; temp100=-; printed=-; moment=-; moment100=-
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
        value=$(grep 'temperature[[:space:]]*=' "$md" | awk '{v[NR]=$(NF-1)} END{first=(NR>100?NR-99:1);for(i=first;i<=NR;i++){s+=v[i];n++}if(n)printf "%.2f",s/n}'); [[ -n "$value" ]] && temp100=$value
    elif [[ -n "$latest" ]]; then
        if grep -Eq 'Error in routine|convergence NOT achieved|Cholesky|cdiaghg' "$latest"; then status=FAILED_SCF
        else status=PRECONDITIONING; fi
    fi
    if [[ -n "$latest" ]]; then
        value=$(grep 'estimated scf accuracy' "$latest" | tail -1 | awk '{print $(NF-1)}'); [[ -n "$value" ]] && accuracy=$value
        value=$(vector_mean "$latest" 'constrained moment'); [[ -n "$value" ]] && printed=$value
        value=$(vector_mean "$latest" 'magnetization'); [[ -n "$value" ]] && moment=$value
        value=$(vector_tail_mean "$latest" 'magnetization' 100); [[ -n "$value" ]] && moment100=$value
    fi
    printf '%-30s %-7s %-18s %-8s %-8s %-10s %-15s %-12s %-12s %-12s %-11s %-12s\n' \
        "$case_name" "$l3" "$status" "$steps" "$total_iter" "$max_cycle" "$accuracy" "$temp" "$temp100" "$printed" "$moment" "$moment100"
done < case_manifest.csv

if [[ -f round12_status.tsv ]]; then
    echo
    column -t -s $'\t' round12_status.tsv 2>/dev/null || cat round12_status.tsv
fi
