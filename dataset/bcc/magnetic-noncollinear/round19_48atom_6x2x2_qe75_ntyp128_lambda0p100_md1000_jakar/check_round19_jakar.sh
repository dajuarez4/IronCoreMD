#!/usr/bin/env bash
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"

case_name=$(tail -n +2 case_manifest.csv | cut -d, -f2)
requested_steps=$(tail -n +2 case_manifest.csv | cut -d, -f10)
folder="cases/$case_name"
md="$folder/md_1000steps.out"
status=NOT_STARTED
steps=0
iterations=0
max_cycle=0
accuracy=-
temperature=-
temperature_mean100=-
moment=-
moment_mean100=-

latest=""
for stem in stage1_4e-3 stage2_1e-3 stage3_3e-4; do
    [[ -f "$folder/$stem.out" ]] && latest="$folder/$stem.out"
done

if [[ -f "$md" ]]; then
    latest="$md"
    steps=$(grep -c 'Entering Dynamics' "$md" || true)
    iterations=$(grep -c 'iteration #' "$md" || true)
    max_cycle=$(awk '/iteration #/{n++} /convergence has been achieved|convergence NOT achieved/{if(n>m)m=n;n=0} END{print m+0}' "$md")
    if grep -q 'JOB DONE' "$md" && [[ "$steps" -eq "$requested_steps" ]]; then
        status=DONE
    elif grep -Eq 'Error in routine|convergence NOT achieved|Cholesky|cdiaghg' "$md"; then
        status=FAILED_MD
    else
        status=RUNNING_OR_PARTIAL_MD
    fi
    value=$(grep 'temperature[[:space:]]*=' "$md" | tail -1 | awk '{print $(NF-1)}')
    [[ -n "$value" ]] && temperature=$value
    value=$(grep 'temperature[[:space:]]*=' "$md" | awk '{v[NR]=$(NF-1)} END{first=(NR>100?NR-99:1);for(i=first;i<=NR;i++){s+=v[i];n++}if(n)printf "%.2f",s/n}')
    [[ -n "$value" ]] && temperature_mean100=$value
elif [[ -n "$latest" ]]; then
    if grep -Eq 'Error in routine|convergence NOT achieved|Cholesky|cdiaghg' "$latest"; then status=FAILED_SCF
    else status=PRECONDITIONING; fi
fi

if [[ -n "$latest" ]]; then
    value=$(grep 'estimated scf accuracy' "$latest" | tail -1 | awk '{print $(NF-1)}')
    [[ -n "$value" ]] && accuracy=$value
    value=$(awk '/^[[:space:]]*magnetization[[:space:]]*:/ {
        slot=(n%48); v[slot]=sqrt($3*$3+$4*$4+$5*$5); n++
    } END {if(n>=48){s=0;for(i=0;i<48;i++)s+=v[i];printf "%.6f",s/48}}' "$latest")
    [[ -n "$value" ]] && moment=$value
    value=$(awk '/^[[:space:]]*magnetization[[:space:]]*:/ {
        slot=(n%48); v[slot]=sqrt($3*$3+$4*$4+$5*$5); n++;
        if(n%48==0){s=0;for(i=0;i<48;i++)s+=v[i]; frame++; avg[frame]=s/48}
    } END {first=(frame>100?frame-99:1);s=0;c=0;for(i=first;i<=frame;i++){s+=avg[i];c++}if(c)printf "%.6f",s/c}' "$latest")
    [[ -n "$value" ]] && moment_mean100=$value
fi

printf '%-22s %s\n' status "$status"
printf '%-22s %s/%s\n' md_steps "$steps" "$requested_steps"
printf '%-22s %s\n' scf_iterations "$iterations"
printf '%-22s %s\n' max_scf_cycle "$max_cycle"
printf '%-22s %s\n' last_accuracy_Ry "$accuracy"
printf '%-22s %s\n' last_temperature_K "$temperature"
printf '%-22s %s\n' temperature_mean100_K "$temperature_mean100"
printf '%-22s %s\n' last_mean_moment "$moment"
printf '%-22s %s\n' moment_mean100 "$moment_mean100"
