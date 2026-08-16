#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
vector_mean() {
 awk -v key="$2" '$0 ~ "^[[:space:]]*" key "[[:space:]]*:" {for(i=1;i<=NF;i++)if($i==":"){n++;slot=(n-1)%8;v[slot]=sqrt($(i+1)^2+$(i+2)^2+$(i+3)^2);break}} END{if(n>=8){s=0;for(i=0;i<8;i++)s+=v[i];printf "%.6f",s/8}}' "$1"
}
printf "%-27s %-8s %-17s %-10s %-15s %-12s %-12s %-10s\n" case final_lam status iterations accuracy target_print measured retention
while IFS=, read -r index case_name lattice target start l1 l2 l3 requirement; do
 [[ "$index" == serial_index ]] && continue
 folder="cases/$case_name"; latest=""; complete=0; state=NOT_STARTED
 for stem in stage1_4e-3 stage2_1e-3 stage3_3e-4; do
  out="$folder/$stem.out"; [[ -f "$out" ]] || continue; latest="$out"
  grep -q 'convergence has been achieved' "$out" && complete=$((complete+1))
  grep -Eq 'Error in routine|convergence NOT achieved' "$out" && state=FAILED_$stem
 done
 if [[ "$state" == NOT_STARTED && -n "$latest" ]]; then [[ "$complete" -eq 3 ]] && state=DONE || state=SCF_$((complete+1))_ACTIVE; fi
 iterations=-; accuracy=-; printed=-; measured=-; retention=-
 if [[ -n "$latest" ]]; then
  iterations=$(grep -c 'iteration #' "$latest" || true)
  value=$(grep 'estimated scf accuracy' "$latest" | tail -1 | awk '{print $(NF-1)}'); [[ -n "$value" ]] && accuracy=$value
  value=$(vector_mean "$latest" 'constrained moment'); [[ -n "$value" ]] && printed=$value
  value=$(vector_mean "$latest" 'magnetization'); [[ -n "$value" ]] && measured=$value
  [[ "$measured" != - ]] && retention=$(awk -v m="$measured" -v t="$target" 'BEGIN{printf "%.4f",m/t}')
 fi
 printf "%-27s %-8s %-17s %-10s %-15s %-12s %-12s %-10s\n" "$case_name" "$l3" "$state" "$iterations" "$accuracy" "$printed" "$measured" "$retention"
done < case_manifest.csv
if [[ -s round8_serial_status.tsv ]]; then echo; column -t -s $'\t' round8_serial_status.tsv 2>/dev/null || cat round8_serial_status.tsv; fi
