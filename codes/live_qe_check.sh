#!/bin/bash

# Usage:
#   bash plot_qe_live.sh fe1.out

if [ $# -lt 1 ]; then
    echo "Usage: bash plot_qe_live.sh QE_output_file"
    exit 1
fi

OUTFILE="$1"

if [ ! -f "$OUTFILE" ]; then
    echo "Error: file '$OUTFILE' not found."
    exit 1
fi

if ! command -v gnuplot >/dev/null 2>&1; then
    echo "Error: gnuplot is not installed."
    exit 1
fi

SCF_DAT="scf_convergence.dat"
ENER_DAT="energy_vs_step.dat"
TMAG_DAT="total_magnetization.dat"
AMAG_DAT="absolute_magnetization.dat"
PRES_DAT="pressure_vs_step.dat"
FORCE_DAT="total_force_vs_step.dat"
PNG_OUT="qe_live_dashboard.png"
TEMP_DAT="temperature_vs_step.dat"

UPDATE_SEC=10
RY_TO_EV=13.605693009
KBAR_TO_GPA=0.1

while true; do
    rm -f "$SCF_DAT" "$ENER_DAT" "$TMAG_DAT" "$AMAG_DAT" "$PRES_DAT" "$FORCE_DAT" "$TEMP_DAT"

    awk -v ry2ev="$RY_TO_EV" -v kbar2gpa="$KBAR_TO_GPA" '
    BEGIN{
        scf_iter = ""
        ener_step = 0
        mag_step = 0
        pres_step = 0
        force_step = 0
        temp_step = 0
        nat = 0
    }

    /number of atoms\/cell[[:space:]]*=/ {
        if (match($0, /=[[:space:]]*[0-9]+/)) {
            tmp = substr($0, RSTART, RLENGTH)
            sub(/=/, "", tmp)
            gsub(/[[:space:]]/, "", tmp)
            nat = tmp + 0
        }
    }

    /iteration #[[:space:]]*[0-9]+/ {
        if (match($0, /iteration #[[:space:]]*[0-9]+/)) {
            scf_iter = substr($0, RSTART, RLENGTH)
            sub(/iteration #[[:space:]]*/, "", scf_iter)
        }
    }

    /estimated scf accuracy/ {
        if (match($0, /<[[:space:]]*[-+]?[0-9]*\.?[0-9]+([Ee][-+]?[0-9]+)?/)) {
            acc = substr($0, RSTART, RLENGTH)
            sub(/</, "", acc)
            gsub(/[[:space:]]/, "", acc)
            if (scf_iter != "" && acc != "") {
                print scf_iter + 0, acc + 0 >> "'"$SCF_DAT"'"
            }
        }
    }

    /^![[:space:]]+total energy[[:space:]]*=/ {
        if (match($0, /=[[:space:]]*[-+]?[0-9]*\.?[0-9]+([Ee][-+]?[0-9]+)?/)) {
            Eval = substr($0, RSTART, RLENGTH)
            sub(/=/, "", Eval)
            gsub(/[[:space:]]/, "", Eval)
            ener_step++
            if (nat > 0) {
                EevAtom = (Eval + 0) * ry2ev / nat
                print ener_step, EevAtom >> "'"$ENER_DAT"'"
            }
        }
    }

    /total magnetization[[:space:]]*=/ {
        line = $0
        sub(/^.*=/, "", line)
        gsub(/Bohr mag\/cell/, "", line)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)

        n = split(line, vals, /[[:space:]]+/)
        if (n >= 3) {
            mag_step++
            print mag_step, vals[1]+0, vals[2]+0, vals[3]+0 >> "'"$TMAG_DAT"'"
        }
    }

    /absolute magnetization[[:space:]]*=/ {
        if (match($0, /=[[:space:]]*[-+]?[0-9]*\.?[0-9]+([Ee][-+]?[0-9]+)?/)) {
            amag = substr($0, RSTART, RLENGTH)
            sub(/=/, "", amag)
            gsub(/[[:space:]]/, "", amag)
            if (mag_step > 0) {
                print mag_step, amag + 0 >> "'"$AMAG_DAT"'"
            }
        }
    }

    /\(kbar\)[[:space:]]*P=/ {
        if (match($0, /P=[[:space:]]*[-+]?[0-9]*\.?[0-9]+([Ee][-+]?[0-9]+)?/)) {
            pval = substr($0, RSTART, RLENGTH)
            sub(/P=/, "", pval)
            gsub(/[[:space:]]/, "", pval)
            pres_step++
            print pres_step, (pval + 0) * kbar2gpa >> "'"$PRES_DAT"'"
        }
    }

    /Total force =/ {
        if (match($0, /Total force =[[:space:]]*[-+]?[0-9]*\.?[0-9]+([Ee][-+]?[0-9]+)?/)) {
            fval = substr($0, RSTART, RLENGTH)
            sub(/Total force =/, "", fval)
            gsub(/[[:space:]]/, "", fval)
            force_step++
            print force_step, fval + 0 >> "'"$FORCE_DAT"'"
        }
    }

    /(temperature|temp)[[:space:]]*=/ {
    if (match($0, /(temperature|temp)[[:space:]]*=[[:space:]]*[-+]?[0-9]*\.?[0-9]+([Ee][-+]?[0-9]+)?/)) {
        tval = substr($0, RSTART, RLENGTH)
        sub(/.*=/, "", tval)
        gsub(/[[:space:]]/, "", tval)
        temp_step++
        print temp_step, tval + 0 >> "'"$TEMP_DAT"'"
        }
    }
    ' "$OUTFILE"

    echo "----------------------------------------"
    echo "Updated at: $(date)"
    [ -f "$SCF_DAT" ]   && echo "SCF points:  $(wc -l < "$SCF_DAT")"
    [ -f "$ENER_DAT" ]  && echo "ENER points: $(wc -l < "$ENER_DAT")"
    [ -f "$TMAG_DAT" ]  && echo "TMAG points: $(wc -l < "$TMAG_DAT")"
    [ -f "$AMAG_DAT" ]  && echo "AMAG points: $(wc -l < "$AMAG_DAT")"
    [ -f "$PRES_DAT" ]  && echo "PRES points: $(wc -l < "$PRES_DAT")"
    [ -f "$FORCE_DAT" ] && echo "FORCE points: $(wc -l < "$FORCE_DAT")"
    [ -f "$TEMP_DAT" ]  && echo "TEMP points:  $(wc -l < "$TEMP_DAT")"

    echo "Last SCF points:"
    [ -s "$SCF_DAT" ] && tail -3 "$SCF_DAT"

    echo "Last energy points:"
    [ -s "$ENER_DAT" ] && tail -3 "$ENER_DAT"

    echo "Last total magnetization points:"
    [ -s "$TMAG_DAT" ] && tail -3 "$TMAG_DAT"

    echo "Last absolute magnetization points:"
    [ -s "$AMAG_DAT" ] && tail -3 "$AMAG_DAT"

    echo "Last pressure points:"
    [ -s "$PRES_DAT" ] && tail -3 "$PRES_DAT"

    echo "Last force points:"
    [ -s "$FORCE_DAT" ] && tail -3 "$FORCE_DAT"

    echo "Last temperature points:"
    [ -s "$TEMP_DAT" ] && tail -3 "$TEMP_DAT"

    if [ -s "$SCF_DAT" ]; then
        SCF_PLOT="plot '$SCF_DAT' using 1:2 with linespoints lw 2 pt 7"
    else
        SCF_PLOT="plot 1/0 notitle"
    fi

    if [ -s "$ENER_DAT" ]; then
        ENER_PLOT="plot '$ENER_DAT' using 1:2 with linespoints lw 2 pt 7"
    else
        ENER_PLOT="plot 1/0 notitle"
    fi

    if [ -s "$TMAG_DAT" ]; then
        TMAG_PLOT="plot '$TMAG_DAT' using 1:2 with linespoints lw 2 pt 7 title 'Mx', \
                         '$TMAG_DAT' using 1:3 with linespoints lw 2 pt 7 title 'My', \
                         '$TMAG_DAT' using 1:4 with linespoints lw 2 pt 7 title 'Mz'"
    else
        TMAG_PLOT="plot 1/0 notitle"
    fi

    if [ -s "$AMAG_DAT" ]; then
        AMAG_PLOT="plot '$AMAG_DAT' using 1:2 with linespoints lw 2 pt 7"
    else
        AMAG_PLOT="plot 1/0 notitle"
    fi

    if [ -s "$PRES_DAT" ]; then
        PRES_PLOT="plot '$PRES_DAT' using 1:2 with linespoints lw 2 pt 7"
    else
        PRES_PLOT="plot 1/0 notitle"
    fi

    if [ -s "$FORCE_DAT" ]; then
        FORCE_PLOT="plot '$FORCE_DAT' using 1:2 with linespoints lw 2 pt 7"
    else
        FORCE_PLOT="plot 1/0 notitle"
    fi

    if [ -s "$TEMP_DAT" ]; then
        TEMP_PLOT="plot '$TEMP_DAT' using 1:2 with linespoints lw 2 pt 7"
    else
        TEMP_PLOT="plot 1/0 notitle"
    fi
    gnuplot <<EOF
set terminal pngcairo size 1600,900
set output "$PNG_OUT"
set multiplot layout 4,2 title sprintf("Live QE Monitoring: %s", "$OUTFILE")
set grid

unset key
set title "SCF Convergence vs Iteration"
set xlabel "SCF Iteration"
set ylabel "Estimated SCF Accuracy (Ry)"
set logscale y
$SCF_PLOT

unset logscale y
unset key
set title "Converged Total Energy vs Step"
set xlabel "Step"
set ylabel "Total Energy (eV/atom)"
$ENER_PLOT

set key top right
set title "Total Magnetization vs Step"
set xlabel "Step"
set ylabel "Total Magnetization (Bohr mag/cell)"
$TMAG_PLOT

unset key
set title "Absolute Magnetization vs Step"
set xlabel "Step"
set ylabel "Absolute Magnetization (Bohr mag/cell)"
$AMAG_PLOT

unset key
set title "Pressure vs Step"
set xlabel "Step"
set ylabel "Pressure (GPa)"
$PRES_PLOT

unset key
set title "Total Force vs Step"
set xlabel "Step"
set ylabel "Total Force (Ry/au)"
$FORCE_PLOT

unset key
set title "Temperature vs Step"
set xlabel "Step"
set ylabel "Temperature (K)"
$TEMP_PLOT

unset multiplot
EOF

    ls -lh "$PNG_OUT" 2>/dev/null
    sleep "$UPDATE_SEC"
done
