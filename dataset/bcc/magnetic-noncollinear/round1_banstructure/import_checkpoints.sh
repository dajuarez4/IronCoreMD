#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
python3 prepare_checkpoint.py --round 24 --source-case /scratch/dajuarez4/testing/non_collinear/round24_48atom_6x2x2_round23settings_nraise20_md400_jakar/cases/conv1e-4_beta1e-2_48atom
python3 prepare_checkpoint.py --round 25 --source-case /scratch/dajuarez4/testing/non_collinear/round25_32atom_round22final_nraise20_md400_jakar/cases/conv1e-4_beta1e-2_32atom
