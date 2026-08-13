#!/usr/bin/env bash

# ================================================================
# Live Quantum ESPRESSO monitoring dashboard for VS Code image viewers
# Compatible with mawk, gawk, and standard awk implementations.
#
# Usage:
#   bash live_qe_en_conv_vscode.sh fe1.out
#   bash live_qe_en_conv_vscode.sh fe1.out --once
#
# This variant renders each dashboard to a temporary PNG and then copies its
# contents over the existing visible PNG. Keeping the destination inode and
# explicitly touching it helps GeriYoco Image Gallery detect remote changes.
#
# Magnetization support:
#   Collinear total magnetization:       step M
#   Noncollinear total magnetization:    step Mx My Mz
#   Absolute magnetization:              step |M|
# ================================================================

set -u
export LC_ALL=C

if [ "$#" -lt 1 ]; then
    echo "Usage: bash $0 QE_output_file [--once]"
    exit 1
fi

OUTFILE="$1"
RUN_ONCE=0

if [ "${2:-}" = "--once" ]; then
    RUN_ONCE=1
fi

if [ ! -f "$OUTFILE" ]; then
    echo "Error: file '$OUTFILE' not found."
    exit 1
fi

if ! command -v awk >/dev/null 2>&1; then
    echo "Error: awk is not installed."
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
TEMP_DAT="temperature_vs_step.dat"
EMPTY_DAT=".qe_empty_plot.dat"
PNG_OUT="qe_live_dashboard_vscode.png"
PNG_TMP=".qe_live_dashboard_vscode.$$.png"

UPDATE_SEC=10
RY_TO_EV=13.605693009
KBAR_TO_GPA=0.1

cleanup() {
    rm -f "$PNG_TMP"
}

stop_monitor() {
    echo
    echo "Live monitoring stopped."
    exit 0
}

trap cleanup EXIT
trap stop_monitor INT TERM

while true; do
    rm -f "$SCF_DAT" "$ENER_DAT" "$TMAG_DAT" "$AMAG_DAT" \
          "$PRES_DAT" "$FORCE_DAT" "$TEMP_DAT"

    # A valid invisible point prevents gnuplot from printing
    # "all points y value undefined" when a dataset is empty.
    printf "1 1\n" > "$EMPTY_DAT"

    awk \
        -v ry2ev="$RY_TO_EV" \
        -v kbar2gpa="$KBAR_TO_GPA" \
        -v scf_file="$SCF_DAT" \
        -v energy_file="$ENER_DAT" \
        -v tmag_file="$TMAG_DAT" \
        -v amag_file="$AMAG_DAT" \
        -v pressure_file="$PRES_DAT" \
        -v force_file="$FORCE_DAT" \
        -v temperature_file="$TEMP_DAT" '

    function numeric_value(value) {
        gsub(/[dD]/, "E", value)
        return value + 0
    }

    BEGIN {
        scf_iter = ""
        scf_record = 0
        energy_step = 0
        total_mag_step = 0
        absolute_mag_step = 0
        pressure_step = 0
        force_step = 0
        temperature_step = 0
        nat = 0
    }

    /number of atoms\/cell[[:space:]]*=/ {
        if (match($0, /=[[:space:]]*[0-9]+/)) {
            value = substr($0, RSTART, RLENGTH)
            sub(/=/, "", value)
            gsub(/[[:space:]]/, "", value)
            nat = value + 0
        }
    }

    /iteration #[[:space:]]*[0-9]+/ {
        if (match($0, /iteration #[[:space:]]*[0-9]+/)) {
            scf_iter = substr($0, RSTART, RLENGTH)
            sub(/iteration #[[:space:]]*/, "", scf_iter)
        }
    }

    /estimated scf accuracy/ {
        if (match($0, /<[[:space:]]*[-+]?[0-9]*\.?[0-9]+([EeDd][-+]?[0-9]+)?/)) {
            accuracy = substr($0, RSTART, RLENGTH)
            sub(/</, "", accuracy)
            gsub(/[[:space:]]/, "", accuracy)
            if (scf_iter != "" && accuracy != "") {
                scf_record++
                print scf_record, numeric_value(accuracy), scf_iter + 0 >> scf_file
            }
        }
    }

    /^[[:space:]]*![[:space:]]+total energy[[:space:]]*=/ {
        if (match($0, /=[[:space:]]*[-+]?[0-9]*\.?[0-9]+([EeDd][-+]?[0-9]+)?/)) {
            energy = substr($0, RSTART, RLENGTH)
            sub(/=/, "", energy)
            gsub(/[[:space:]]/, "", energy)
            energy_step++
            if (nat > 0) {
                energy_ev_atom = numeric_value(energy) * ry2ev / nat
                print energy_step, energy_ev_atom >> energy_file
            }
        }
    }

    # Total magnetization:
    #   scalar collinear format:    total magnetization = M
    #   vector noncollinear format: total magnetization = Mx My Mz
    /^[[:space:]]*total magnetization[[:space:]]*=/ {
        line = $0
        sub(/^.*=/, "", line)
        gsub(/[Bb]ohr[[:space:]]+mag\/cell/, "", line)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
        number_values = split(line, values, /[[:space:]]+/)

        if (number_values >= 3) {
            total_mag_step++
            print total_mag_step, numeric_value(values[1]), numeric_value(values[2]), numeric_value(values[3]) >> tmag_file
        } else if (number_values >= 1 && values[1] != "") {
            total_mag_step++
            print total_mag_step, numeric_value(values[1]) >> tmag_file
        }
    }

    # Absolute magnetization is stored independently from total M.
    /^[[:space:]]*absolute magnetization[[:space:]]*=/ {
        line = $0
        sub(/^.*=/, "", line)
        gsub(/[Bb]ohr[[:space:]]+mag\/cell/, "", line)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", line)
        number_values = split(line, values, /[[:space:]]+/)

        if (number_values >= 1 && values[1] != "") {
            absolute_mag_step++
            print absolute_mag_step, numeric_value(values[1]) >> amag_file
        }
    }

    /\(kbar\)[[:space:]]*P=/ {
        if (match($0, /P=[[:space:]]*[-+]?[0-9]*\.?[0-9]+([EeDd][-+]?[0-9]+)?/)) {
            pressure = substr($0, RSTART, RLENGTH)
            sub(/P=/, "", pressure)
            gsub(/[[:space:]]/, "", pressure)
            pressure_step++
            pressure_gpa = numeric_value(pressure) * kbar2gpa
            print pressure_step, pressure_gpa >> pressure_file
        }
    }

    /Total force[[:space:]]*=/ {
        if (match($0, /Total force[[:space:]]*=[[:space:]]*[-+]?[0-9]*\.?[0-9]+([EeDd][-+]?[0-9]+)?/)) {
            force = substr($0, RSTART, RLENGTH)
            sub(/Total force[[:space:]]*=/, "", force)
            gsub(/[[:space:]]/, "", force)
            force_step++
            print force_step, numeric_value(force) >> force_file
        }
    }

    /(temperature|temp)[[:space:]]*=/ {
        if (match($0, /(temperature|temp)[[:space:]]*=[[:space:]]*[-+]?[0-9]*\.?[0-9]+([EeDd][-+]?[0-9]+)?/)) {
            temperature = substr($0, RSTART, RLENGTH)
            sub(/.*=/, "", temperature)
            gsub(/[[:space:]]/, "", temperature)
            temperature_step++
            print temperature_step, numeric_value(temperature) >> temperature_file
        }
    }

    ' "$OUTFILE"

    echo
    echo "================================================"
    echo "Quantum ESPRESSO live monitor"
    echo "Output file: $OUTFILE"
    echo "Updated:     $(date)"
    echo "================================================"

    if [ -s "$SCF_DAT" ]; then
        echo "SCF records:                  $(wc -l < "$SCF_DAT")"
    else
        echo "SCF records:                  0"
    fi

    if [ -s "$ENER_DAT" ]; then
        echo "Energy records:               $(wc -l < "$ENER_DAT")"
    else
        echo "Energy records:               0"
    fi

    if [ -s "$TMAG_DAT" ]; then
        TMAG_NCOL=$(awk 'NF {print NF; exit}' "$TMAG_DAT")
        echo "Total magnetization records:  $(wc -l < "$TMAG_DAT")"
        if [ "$TMAG_NCOL" -ge 4 ]; then
            echo "Total magnetization format:   noncollinear Mx My Mz"
        else
            echo "Total magnetization format:   collinear scalar M"
        fi
    else
        TMAG_NCOL=0
        echo "Total magnetization records:  0"
        echo "Total magnetization format:   not detected yet"
    fi

    if [ -s "$AMAG_DAT" ]; then
        echo "Absolute magnetization:       $(wc -l < "$AMAG_DAT")"
    else
        echo "Absolute magnetization:       0"
    fi

    if [ -s "$PRES_DAT" ]; then
        echo "Pressure records:             $(wc -l < "$PRES_DAT")"
    else
        echo "Pressure records:             0"
    fi

    if [ -s "$FORCE_DAT" ]; then
        echo "Force records:                $(wc -l < "$FORCE_DAT")"
    else
        echo "Force records:                0"
    fi

    if [ -s "$TEMP_DAT" ]; then
        echo "Temperature records:          $(wc -l < "$TEMP_DAT")"
    else
        echo "Temperature records:          0"
    fi

    echo "------------------------------------------------"

    echo "Last SCF records:"
    if [ -s "$SCF_DAT" ]; then tail -3 "$SCF_DAT"; else echo "No SCF values found."; fi

    echo
    echo "Last energy records:"
    if [ -s "$ENER_DAT" ]; then tail -3 "$ENER_DAT"; else echo "No converged energies found yet."; fi

    echo
    echo "Last total magnetization records:"
    if [ -s "$TMAG_DAT" ]; then tail -3 "$TMAG_DAT"; else echo "No total magnetization values found."; fi

    echo
    echo "Last absolute magnetization records:"
    if [ -s "$AMAG_DAT" ]; then tail -3 "$AMAG_DAT"; else echo "No absolute magnetization values found."; fi

    echo
    echo "Last pressure records:"
    if [ -s "$PRES_DAT" ]; then tail -3 "$PRES_DAT"; else echo "No pressure values found yet."; fi

    echo
    echo "Last force records:"
    if [ -s "$FORCE_DAT" ]; then tail -3 "$FORCE_DAT"; else echo "No force values found yet."; fi

    echo
    echo "Last temperature records:"
    if [ -s "$TEMP_DAT" ]; then tail -3 "$TEMP_DAT"; else echo "No temperature values found yet."; fi

    NO_DATA_PLOT="set label 99 'No data yet' at graph 0.5,0.5 center; plot '$EMPTY_DAT' using 1:2 with points ps 0 notitle; unset label 99"

    if [ -s "$SCF_DAT" ]; then
        SCF_PLOT="plot '$SCF_DAT' using 1:2 with linespoints lw 2 pt 7 title 'SCF accuracy'"
    else
        SCF_PLOT="$NO_DATA_PLOT"
    fi

    if [ -s "$ENER_DAT" ]; then
        ENER_PLOT="plot '$ENER_DAT' using 1:2 with linespoints lw 2 pt 7 title 'Energy'"
    else
        ENER_PLOT="$NO_DATA_PLOT"
    fi

    if [ -s "$TMAG_DAT" ]; then
        if [ "$TMAG_NCOL" -ge 4 ]; then
            TMAG_PLOT="plot '$TMAG_DAT' using 1:2 with linespoints lw 2 pt 7 title 'Mx', '$TMAG_DAT' using 1:3 with linespoints lw 2 pt 7 title 'My', '$TMAG_DAT' using 1:4 with linespoints lw 2 pt 7 title 'Mz'"
        else
            TMAG_PLOT="plot '$TMAG_DAT' using 1:2 with linespoints lw 2 pt 7 title 'Total M'"
        fi
    else
        TMAG_PLOT="$NO_DATA_PLOT"
    fi

    if [ -s "$AMAG_DAT" ]; then
        AMAG_PLOT="plot '$AMAG_DAT' using 1:2 with linespoints lw 2 pt 7 title 'Absolute M'"
    else
        AMAG_PLOT="$NO_DATA_PLOT"
    fi

    if [ -s "$PRES_DAT" ]; then
        PRES_PLOT="plot '$PRES_DAT' using 1:2 with linespoints lw 2 pt 7 title 'Pressure'"
    else
        PRES_PLOT="$NO_DATA_PLOT"
    fi

    if [ -s "$FORCE_DAT" ]; then
        FORCE_PLOT="plot '$FORCE_DAT' using 1:2 with linespoints lw 2 pt 7 title 'Total force'"
    else
        FORCE_PLOT="$NO_DATA_PLOT"
    fi

    if [ -s "$TEMP_DAT" ]; then
        TEMP_PLOT="plot '$TEMP_DAT' using 1:2 with linespoints lw 2 pt 7 title 'Temperature'"
    else
        TEMP_PLOT="$NO_DATA_PLOT"
    fi

    gnuplot <<EOF
set terminal pngcairo size 1600,1000 enhanced
set output "$PNG_TMP"
set multiplot layout 4,2 rowsfirst title sprintf("Live Quantum ESPRESSO Monitoring: %s", "$OUTFILE")
set grid
set border linewidth 1.2
set tics out
set datafile separator whitespace

unset key
set title "SCF Convergence"
set xlabel "SCF Record"
set ylabel "Estimated SCF Accuracy (Ry)"
set logscale y
$SCF_PLOT

unset logscale y
unset key
set title "Converged Total Energy"
set xlabel "Energy Record"
set ylabel "Total Energy (eV/atom)"
$ENER_PLOT

set key top right
set title "Total Magnetization"
set xlabel "Magnetization Record"
set ylabel "Magnetization (Bohr mag/cell)"
$TMAG_PLOT

set key top right
set title "Absolute Magnetization"
set xlabel "Magnetization Record"
set ylabel "Absolute Magnetization (Bohr mag/cell)"
$AMAG_PLOT

unset key
set title "Pressure"
set xlabel "Pressure Record"
set ylabel "Pressure (GPa)"
$PRES_PLOT

unset key
set title "Total Force"
set xlabel "Force Record"
set ylabel "Total Force (Ry/au)"
$FORCE_PLOT

unset key
set title "Temperature"
set xlabel "Temperature Record"
set ylabel "Temperature (K)"
$TEMP_PLOT

unset key
unset grid
unset border
unset tics
unset xlabel
unset ylabel
set title ""
set label 100 "Waiting for additional QE data" at graph 0.5,0.5 center
plot "$EMPTY_DAT" using 1:2 with points ps 0 notitle
unset label 100

unset multiplot
EOF

    if [ ! -s "$PNG_TMP" ]; then
        echo "Error: gnuplot did not create a valid temporary dashboard." >&2
        exit 2
    fi

    # Preserve the visible file path/inode while updating its contents. Some
    # custom VS Code image editors watch "modified" but miss atomic replacement
    # events over Remote SSH.
    cp -f "$PNG_TMP" "$PNG_OUT"
    rm -f "$PNG_TMP"
    touch "$PNG_OUT"

    echo
    ls -lh "$PNG_OUT" 2>/dev/null || true

    if [ "$RUN_ONCE" -eq 1 ]; then
        break
    fi

    sleep "$UPDATE_SEC"
done
