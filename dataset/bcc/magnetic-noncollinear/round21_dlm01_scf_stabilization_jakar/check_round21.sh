#!/usr/bin/env bash
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"

printf '%-24s %-12s %7s %8s %10s %10s %10s\n' case status steps warnings scf_iters last_T_K mean50_K
while IFS= read -r case_name; do
    md="cases/$case_name/md_100steps.out"
    status=NOT_STARTED; steps=0; warnings=0; iterations=0; last=-; mean50=-
    if [[ -f "$md" ]]; then
        steps=$(grep -c 'Entering Dynamics' "$md" || true)
        warnings=$(grep -c 'SCF correction compared to forces is large' "$md" || true)
        iterations=$(grep -c 'iteration #' "$md" || true)
        if grep -q 'convergence NOT achieved' "$md"; then status=FAILED
        elif grep -q 'JOB DONE' "$md" && [[ "$steps" -eq 100 ]]; then status=DONE
        else status=RUNNING; fi
        value=$(awk '/temperature[[:space:]]*=/{v[++n]=$(NF-1)} END{if(n)printf "%.0f",v[n]}' "$md")
        [[ -n "$value" ]] && last=$value
        value=$(awk '/temperature[[:space:]]*=/{v[++n]=$(NF-1)} END{f=(n>50?n-49:1);for(i=f;i<=n;i++){s+=v[i];m++}if(m)printf "%.0f",s/m}' "$md")
        [[ -n "$value" ]] && mean50=$value
    fi
    printf '%-24s %-12s %7s %8s %10s %10s %10s\n' "$case_name" "$status" "$steps" "$warnings" "$iterations" "$last" "$mean50"
done < cases.txt
