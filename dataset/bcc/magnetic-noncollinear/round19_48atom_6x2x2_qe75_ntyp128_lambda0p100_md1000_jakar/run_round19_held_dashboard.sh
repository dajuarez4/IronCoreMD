#!/usr/bin/env bash
set -euo pipefail

# Complete reproducible Round-19 postprocessing entry point.
# Optional environment controls:
#   ROUND19_PYTHON=/path/to/python
#   ROUND19_STRIDE=1       # use 10 for a future 1000-step output
#   ROUND19_FPS=5
#   ROUND19_DPI=150
#   ROUND19_QE_OUTPUT=/path/to/md_1000steps.out

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PYTHON_BIN=${ROUND19_PYTHON:-/opt/anaconda3/bin/python}
STRIDE=${ROUND19_STRIDE:-1}
FPS=${ROUND19_FPS:-5}
DPI=${ROUND19_DPI:-150}
QE_OUTPUT=${ROUND19_QE_OUTPUT:-$SCRIPT_DIR/md_1000steps.out}
QE_INPUT=${ROUND19_QE_INPUT:-$SCRIPT_DIR/cases/a2p39_6x2x2_m2_48species_lambda0p100_md1000/md_1000steps.in}

if [[ ! -x "$PYTHON_BIN" ]]; then
    echo "Python executable not found: $PYTHON_BIN" >&2
    exit 2
fi
if [[ ! -f "$QE_OUTPUT" ]]; then
    echo "QE output not found: $QE_OUTPUT" >&2
    exit 2
fi
if [[ ! -f "$QE_INPUT" ]]; then
    echo "QE input not found: $QE_INPUT" >&2
    exit 2
fi

"$PYTHON_BIN" - <<'PY'
import importlib.util
missing = [name for name in ("numpy", "scipy", "matplotlib", "PIL", "ase", "spglib")
           if importlib.util.find_spec(name) is None]
if missing:
    raise SystemExit("Missing Python dependencies: " + ", ".join(missing))
PY

cd "$SCRIPT_DIR"
"$PYTHON_BIN" analyze_round19.py --input "$QE_INPUT" --output "$QE_OUTPUT" >/dev/null
"$PYTHON_BIN" make_round19_held_dashboard.py \
    --input "$QE_INPUT" \
    --output "$QE_OUTPUT" \
    --stride "$STRIDE" \
    --fps "$FPS" \
    --dpi "$DPI"

echo
echo "Round-19 HELD dashboard complete"
echo "GIF:  $SCRIPT_DIR/round19_held_md_dashboard.gif"
echo "PNG:  $SCRIPT_DIR/round19_held_md_dashboard_final.png"
echo "HELD: $SCRIPT_DIR/round19_held_summary.json"
