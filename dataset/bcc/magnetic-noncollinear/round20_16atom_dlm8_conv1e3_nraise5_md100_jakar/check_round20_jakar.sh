#!/usr/bin/env bash
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"

printf '%-33s %-22s %7s %9s %10s %10s %8s\n' \
    case status steps last_T_K mean50_K mean_muB warnings

while IFS= read -r case_name; do
    folder="cases/$case_name"
    md="$folder/md_100steps.out"
    latest=""
    status=NOT_STARTED
    steps=0
    last_temperature=-
    mean_temperature=-
    mean_moment=-
    warnings=0

    for stem in stage1_4e-3 stage2_1e-3 stage3_3e-4; do
        [[ -f "$folder/$stem.out" ]] && latest="$folder/$stem.out"
    done
    if [[ -f "$md" ]]; then
        latest="$md"
        steps=$(grep -c 'Entering Dynamics' "$md" || true)
        warnings=$(grep -c 'SCF correction compared to forces is large' "$md" || true)
        if grep -q 'JOB DONE' "$md" && [[ "$steps" -eq 100 ]]; then
            status=DONE
        elif grep -Eq 'Error in routine|convergence NOT achieved|Cholesky|cdiaghg' "$md"; then
            status=FAILED_MD
        else
            status=RUNNING_OR_PARTIAL_MD
        fi
        value=$(grep 'temperature[[:space:]]*=' "$md" | tail -1 | awk '{print $(NF-1)}')
        [[ -n "$value" ]] && last_temperature=$(printf '%.0f' "$value")
        value=$(grep 'temperature[[:space:]]*=' "$md" | awk '{v[NR]=$(NF-1)} END{first=(NR>50?NR-49:1);for(i=first;i<=NR;i++){s+=v[i];n++}if(n)printf "%.0f",s/n}')
        [[ -n "$value" ]] && mean_temperature=$value
    elif [[ -n "$latest" ]]; then
        if grep -Eq 'Error in routine|convergence NOT achieved|Cholesky|cdiaghg' "$latest"; then
            status=FAILED_SCF
        else
            status=PRECONDITIONING
        fi
    fi

    if [[ -n "$latest" ]]; then
        value=$(awk '/^[[:space:]]*magnetization[[:space:]]*:/ {
            slot=(n%16); v[slot]=sqrt($3*$3+$4*$4+$5*$5); n++
        } END {if(n>=16){s=0;for(i=0;i<16;i++)s+=v[i];printf "%.3f",s/16}}' "$latest")
        [[ -n "$value" ]] && mean_moment=$value
    fi
    printf '%-33s %-22s %3s/100 %9s %10s %10s %8s\n' \
        "$case_name" "$status" "$steps" "$last_temperature" \
        "$mean_temperature" "$mean_moment" "$warnings"
done < cases.txt
