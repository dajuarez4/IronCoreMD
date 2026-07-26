# Noncollinear BCC Fe AIMD dataset

This directory contains the processed 54-atom, `3 x 3 x 3` noncollinear BCC
iron trajectory at a nominal ionic temperature of `4000 K`.

## Current trajectory

- Raw Quantum ESPRESSO output: `fe1.out`
- Parsed position records: `65`
- Complete force frames: `64`
- Complete atom-resolved local-spin frames: `64`
- Synchronized MD iterations: `3--66`
- Synchronized trajectory span: `0.063 ps`
- Cell: cubic, `7.65 A`

MD iteration `67` has positions but no completed force or local-magnetization
evaluation. Quantum ESPRESSO stopped during its SCF with a `cdiaghg` Cholesky
error, so that frame is excluded by the validity masks in `simulation.npz`.

## Archive contents

`simulation.npz` includes positions, forces, energies, temperature, pressure,
global magnetization, absolute magnetization, and the atom-resolved
noncollinear vectors in `local_magnetization_Bohr`. The synchronized dashboard
uses

```text
frame_valid & local_magnetization_frame_valid
```

to prevent incomplete ionic, force, spin, or phonon records from being mixed.
`md_processing_summary.json` records the frame counts and termination state.

## HELD diagnostic

`held_heatmap_steps_all_complete.npz` stores the step-resolved HELD phonons for
all 64 complete frames. The cache has shape `(64, 301, 3)`: 64 MD frames, 301
wave-vector samples, and three BCC branches. The arithmetic-mean dispersion
extends to approximately `15.08 THz` and has no negative-frequency points on
the sampled `Gamma-H-N-Gamma-P-H` path.

These frames are consecutive and strongly correlated, and the `3 x 3 x 3`
supercell is smaller than required for a definitive fifth-shell cutoff test.
The result is therefore a trajectory diagnostic, not a converged
finite-temperature phonon spectrum.

## Visualization

The dashboard generator is `../../../codes/make_noncollinear_dashboard.py`.
It renders the atomic trajectory, transparent ideal BCC sites, evolving QE
local-spin arrows, global magnetization, temperature, displacement histories,
and synchronized HELD phonons:

```bash
python codes/make_noncollinear_dashboard.py
```

The generated PNG and GIF are stored in `assets/`.
