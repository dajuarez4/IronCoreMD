#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

for case_dir in a_2.50 a_2.78; do
    echo "Submitting ${case_dir}"
    (
        cd "${root_dir}/${case_dir}"
        sbatch run_jakar.sbatch
    )
done
