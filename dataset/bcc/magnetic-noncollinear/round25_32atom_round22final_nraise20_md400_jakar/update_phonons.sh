#!/usr/bin/env bash
# Run on this Mac after refreshing the selected MD outputs.
set -euo pipefail
export ROUND25_SEGMENTS="${ROUND25_SEGMENTS:-md_400steps.run000.out md_400steps.run002.out}"
python_bin="${PYTHON:-python3}"
work='/Users/dajuarez4/Documents/Fe/dataset/bcc/magnetic-noncollinear/round25'
cd "$work"
mkdir .update_phonons.lock 2>/dev/null || { echo "Another update is running (or remove a stale .update_phonons.lock)."; exit 1; }
trap 'rmdir .update_phonons.lock' EXIT
"$python_bin" update_dashboard.py
"$python_bin" update_held_dashboard.py
"$python_bin" plot_tdep_comparison.py
