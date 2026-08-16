#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
command -v gnuplot >/dev/null || { echo "gnuplot is required" >&2; exit 2; }
mkdir -p .round9_plot

prepare_data() {
    local index=0 case_name folder output
    : > .round9_plot/files.txt
    while IFS=, read -r serial case_name rest; do
        [[ "$serial" == serial_index ]] && continue
        folder="cases/$case_name"; output="$folder/md_10steps.out"
        index=$((index + 1))
        printf '%s\t%s\n' "$index" "$case_name" >> .round9_plot/files.txt
        : > ".round9_plot/accuracy_${index}.dat"
        : > ".round9_plot/cycles_${index}.dat"
        : > ".round9_plot/temp_${index}.dat"
        : > ".round9_plot/moment_${index}.dat"
        [[ -f "$output" ]] || continue
        awk '/estimated scf accuracy/{n++; print n,$(NF-1)}' "$output" > ".round9_plot/accuracy_${index}.dat"
        awk '/iteration #/{n++} /convergence has been achieved|convergence NOT achieved/{c++; print c,n; n=0}' "$output" > ".round9_plot/cycles_${index}.dat"
        awk '/temperature[[:space:]]*=/{n++; print n,$(NF-1)}' "$output" > ".round9_plot/temp_${index}.dat"
        awk '/^[[:space:]]*magnetization[[:space:]]*:/ {
            slot=(n%8); value[slot]=sqrt($3*$3+$4*$4+$5*$5); n++;
            if(n%8==0){sum=0; for(i=0;i<8;i++)sum+=value[i]; record++; print record,sum/8}
        }' "$output" > ".round9_plot/moment_${index}.dat"
    done < case_manifest.csv
}

make_plot() {
    prepare_data
    gnuplot <<'GNUPLOT'
set terminal pngcairo size 1800,1200 enhanced font "Arial,13"
set output "round9_md_dashboard.png"
set multiplot layout 2,2 title "BCC Fe 2x2x1 - Round 9 patched-QE noncollinear MD pilots" font ",20"
set grid ytics lc rgb "#dddddd"
set key top right
set logscale y
set xlabel "Electronic accuracy record"
set ylabel "Estimated SCF accuracy (Ry)"
set xrange [0:100]
set yrange [1e-4:100]
plot ".round9_plot/accuracy_1.dat" u 1:2 w l lw 2 lc rgb "#00a878" title "lambda=0.05", \
     ".round9_plot/accuracy_2.dat" u 1:2 w l lw 2 lc rgb "#9b00ff" title "lambda=0.10", \
     4e-3 w l dt 2 lc rgb "#555555" title "MD conv_thr"

unset logscale y
set autoscale y
set xlabel "SCF cycle (initial plus ionic steps)"
set ylabel "Electronic iterations"
set xrange [0:12]
set yrange [0:25]
plot ".round9_plot/cycles_1.dat" u 1:2 w lp lw 2 pt 7 lc rgb "#00a878" title "lambda=0.05", \
     ".round9_plot/cycles_2.dat" u 1:2 w lp lw 2 pt 7 lc rgb "#9b00ff" title "lambda=0.10"

set xlabel "Printed ionic-temperature record"
set ylabel "Temperature (K)"
set xrange [0:12]
set yrange [1500:6000]
plot ".round9_plot/temp_1.dat" u 1:2 w lp lw 2 pt 7 lc rgb "#00a878" title "lambda=0.05", \
     ".round9_plot/temp_2.dat" u 1:2 w lp lw 2 pt 7 lc rgb "#9b00ff" title "lambda=0.10", \
     4000 w l dt 2 lc rgb "#555555" title "target"

set autoscale y
set xlabel "Printed 8-atom magnetic record"
set ylabel "Mean local moment (Bohr magneton/Fe)"
set xrange [0:12]
set yrange [0:2.2]
plot ".round9_plot/moment_1.dat" u 1:2 w l lw 2 lc rgb "#00a878" title "lambda=0.05", \
     ".round9_plot/moment_2.dat" u 1:2 w l lw 2 lc rgb "#9b00ff" title "lambda=0.10", \
     2.0 w l dt 2 lc rgb "#555555" title "constraint target", \
     1.0 w l dt 3 lc rgb "#999999" title "acceptance guide"
unset multiplot
GNUPLOT
    echo "Updated: $ROOT/round9_md_dashboard.png"
    bash check_round9.sh
}

if [[ "${1:-}" == "--live" ]]; then
    interval=${2:-30}
    while true; do make_plot; sleep "$interval"; done
else
    make_plot
fi
