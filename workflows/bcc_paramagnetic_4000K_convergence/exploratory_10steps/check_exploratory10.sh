#!/bin/bash
set -u

output=md_exploratory_4000K_10steps.out

if [[ ! -e "${output}" ]]; then
    echo "status=NOT_STARTED"
    exit 0
fi

if grep -q "JOB DONE\." "${output}"; then
    status=COMPLETE
elif grep -q "Error in routine\|convergence NOT achieved" "${output}"; then
    status=FAILED
else
    status=RUNNING_OR_INCOMPLETE
fi

echo "status=${status}"
echo "ionic_positions=$(grep -c 'Entering Dynamics:' "${output}" || true)"
echo "converged_scf_cycles=$(grep -c 'convergence has been achieved' "${output}" || true)"
echo "c_bands_warnings=$(grep -c 'c_bands:' "${output}" || true)"
echo "cholesky_errors=$(grep -ci 'Cholesky' "${output}" || true)"
echo "large_force_correction_warnings=$(grep -c 'SCF correction compared to forces is large' "${output}" || true)"

awk '
/Entering Dynamics:/ {ionic=$NF}
/temperature[[:space:]]*=/ {temperature=$(NF-1)}
/estimated scf accuracy/ {accuracy=$(NF-1)}
/total magnetization/ {
    mx=$(NF-4); my=$(NF-3); mz=$(NF-2)
    magnitude=sqrt(mx*mx+my*my+mz*mz)
}
/absolute magnetization/ {absolute=$(NF-2)}
/convergence has been achieved/ {
    printf "ionic=%s temperature_K=%s scf_iterations=%s accuracy_Ry=%s M=(%.2f,%.2f,%.2f) Mnorm=%.3f Mabs=%s\n", ionic, temperature, $6, accuracy, mx, my, mz, magnitude, absolute
}
/Total force =/ {
    ratio=($4 != 0 ? $9/$4 : 0)
    printf "force_Ry_per_bohr=%s scf_correction_Ry_per_bohr=%s correction_ratio=%.3f\n", $4, $9, ratio
}
' "${output}"

grep -E "maximum number of steps|End of molecular dynamics|JOB DONE|Error in routine|convergence NOT achieved" "${output}" | tail -12 || true

last_mnorm=$(awk '
/total magnetization/ {
    mx=$(NF-4); my=$(NF-3); mz=$(NF-2)
    magnitude=sqrt(mx*mx+my*my+mz*mz)
}
/convergence has been achieved/ {last=magnitude}
END {if (last != "") printf "%.6f", last}
' "${output}")

if [[ -n "${last_mnorm}" ]] && awk -v value="${last_mnorm}" 'BEGIN {exit !(value >= 5.0)}'; then
    echo "magnetic_review=REQUIRED: converged Mnorm reached ${last_mnorm} mu_B/cell"
else
    echo "magnetic_review=monitor trend; suggested rejection level is 5 mu_B/cell"
fi
