#!/usr/bin/env bash
set -u
export LC_ALL=C
ROOT=$(cd "$(dirname "$0")" && pwd); cd "$ROOT"
MODE=${1:---once}; INTERVAL=${2:-30}; DATA=.round8_monitor; PNG=round8_dashboard.png; TMP=.round8_dashboard.$$.png
[[ "$MODE" == --once || "$MODE" == --live ]] || { echo "Usage: bash $0 [--once|--live seconds]" >&2; exit 2; }
command -v gnuplot >/dev/null || { echo 'gnuplot is required' >&2; exit 3; }
mkdir -p "$DATA"; trap 'rm -f "$TMP"' EXIT; trap 'exit 0' INT TERM
vmean(){ awk -v key="$2" '$0~"^[[:space:]]*"key"[[:space:]]*:"{for(i=1;i<=NF;i++)if($i==":"){n++;q=(n-1)%8;v[q]=sqrt($(i+1)^2+$(i+2)^2+$(i+3)^2);break}}END{if(n>=8){s=0;for(i=0;i<8;i++)s+=v[i];print s/8}}' "$1"; }
render(){
 : > "$DATA/hist.dat"; : > "$DATA/summary.dat"; done_n=0; active=0; failed=0
 while IFS=, read -r idx case_name lattice target start l1 l2 l3 requirement; do
  [[ "$idx" == serial_index ]] && continue; folder="cases/$case_name"; off=0; latest=""; complete=0; bad=0; ns=()
  for stem in stage1_4e-3 stage2_1e-3 stage3_3e-4; do
   out="$folder/$stem.out"; n=0
   if [[ -f "$out" ]]; then latest="$out"; n=$(grep -c 'iteration #' "$out"||true); awk -v o="$off" '/estimated scf accuracy/{k++;print o+k,$(NF-1)}' "$out" >> "$DATA/hist.dat"; grep -q 'convergence has been achieved' "$out"&&complete=$((complete+1)); grep -Eq 'Error in routine|convergence NOT achieved' "$out"&&bad=1; fi
   ns+=("$n"); off=$((off+n))
  done
  printf '\n\n' >> "$DATA/hist.dat"; [[ "$bad" -eq 1 ]]&&failed=$((failed+1))||{ [[ "$complete" -eq 3 ]]&&done_n=$((done_n+1))||{ [[ -n "$latest" ]]&&active=$((active+1)); }; }
  acc=NaN; printed=NaN; measured=NaN
  if [[ -n "$latest" ]]; then value=$(grep 'estimated scf accuracy' "$latest"|tail -1|awk '{print $(NF-1)}');[[ -n "$value" ]]&&acc=$value; value=$(vmean "$latest" 'constrained moment');[[ -n "$value" ]]&&printed=$value; value=$(vmean "$latest" 'magnetization');[[ -n "$value" ]]&&measured=$value; fi
  printf '%s %s %s %s %s %s %s %s\n' "$idx" "$l3" "${ns[0]}" "${ns[1]}" "${ns[2]}" "$acc" "$printed" "$measured" >> "$DATA/summary.dat"
 done < case_manifest.csv
 maxx=$(awk '{n=$3+$4+$5;if(n>m)m=n}END{print(m?m:10)}' "$DATA/summary.dat")
 gnuplot <<GP
set term pngcairo size 2200,1600 font 'Arial,13';set output '$TMP';set multiplot layout 2,2 title 'Round 8 - patched QE 7.5 atomic constraint validation' font ',20';set grid
set title 'SCF histories';set xlabel 'Cumulative iteration';set ylabel 'Accuracy (Ry)';set logscale y;set xrange [0:$maxx]
plot for[i=0:3]'$DATA/hist.dat' index i u 1:2 w l lw 2 title word('0.01 0.02 0.05 0.10',i+1)
unset logscale y;set autoscale x
set title 'Printed constrained target';set xlabel 'Final lambda';set ylabel 'Bohr magneton/Fe';set logscale x;set xrange[0.008:0.13];set yrange[0:2.3]
plot '$DATA/summary.dat' u 2:7 w lp pt 7 lw 2 title 'printed',2 w l dt 2 title 'intended 2'
set title 'Measured final local moment - DONE: $done_n ACTIVE: $active FAILED: $failed';set ylabel 'Bohr magneton/Fe'
plot '$DATA/summary.dat' u 2:8 w lp pt 7 lw 2 title 'measured',2 w l dt 2 title 'target'
unset logscale x;set autoscale;set title 'Iterations per stage';set xlabel 'Final lambda';set ylabel 'Iterations';set logscale x;set xrange[0.008:0.13]
plot '$DATA/summary.dat' u 2:3 w lp title '4e-3','$DATA/summary.dat' u 2:4 w lp title '1e-3','$DATA/summary.dat' u 2:5 w lp title '3e-4'
unset multiplot
GP
 [[ -f "$PNG" ]]&&{ cp "$TMP" "$PNG";rm -f "$TMP"; }||mv "$TMP" "$PNG";echo "Updated: $ROOT/$PNG";bash check_round8.sh
}
while true;do render;[[ "$MODE" == --once ]]&&break;sleep "$INTERVAL";done
