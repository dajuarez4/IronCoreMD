# Round 22 per-frame HELD phonons

Open `round22_held_md_dashboard.gif` for the animation or
`round22_held_md_dashboard_final.png` for the final dashboard.
All 62 complete force frames in the current local output are fitted and animated.
The unfinished next SCF cycle is excluded. No MD or DFT calculation is launched.

The dashboard shows phonon dispersion, atom positions and local spin directions,
temperature, pressure, cell magnetization components, SCF iteration counts,
relative SCF force correction, displacement from ideal BCC, electronic and
constraint energies, HELD force residuals, and the fraction of sampled path
frequencies that are imaginary. Atom color indicates local spin z direction
(blue negative, red positive); arrows indicate spin direction, not magnitude.
The header lists the principal fixed calculation settings. The complete DFT
parameters remain in `cases/conv1e-4_beta1e-2_32atom/md_100steps.in`.

## Model and alignment

This uses the repository's cubic-symmetry HELD model with one neighbor shell
(two independent coefficients), fitted independently to every force frame.
Its 2.0698 Å cutoff is less than half the shortest cell dimension (2.39 Å).
The existing HELD minimum-image pair mapping cannot unambiguously represent
all five shells in this short 4×2×2 cell, so this calculation deliberately
uses 1NN rather than the nonmagnetic notebook's 5NN setting.
Longer-range interactions are omitted. The per-component force RMSE is
1.55–4.06 eV/Å, so these are limited-model diagnostics, not converged phonons.
The dashed reference dispersion is calculated from the mean coefficients,
not from an arithmetic average of frequencies. Negative frequencies represent
imaginary modes; the plotted percentage is over this sampled symmetry path,
not a Brillouin-zone volume fraction.

QE's first force block belongs to the input coordinates. Force frame k uses
the positions printed after propagation k−1, and its temperature is the
pre-propagation temperature (4000 K for frame 1). Thus the 62 force frames
cover approximately 0–61 fs. The output also prints positions and temperature
after propagation 62 at 62 fs, but these lack the next completed force block
and are excluded from fitting. Initial ideal BCC sites are assigned to the
input atom order with a periodic minimum-distance assignment.

Validation checks require sequential step IDs, converged SCF segments,
complete observable records, full design-matrix rank for every fit, finite
frequencies, a cutoff below the minimum-image bound, acoustic Γ modes below
1e-4 THz, and a GIF frame count equal to the fitted frame count.

## Files and regeneration

- `round22_held_phonon_frames.npz`: all step dispersions, q-points, path distances,
  labels, mean-coefficient dispersion, per-frame HELD coefficients and step IDs.
- `round22_held_trajectory.npz`: aligned coordinates, forces, local spins,
  thermodynamic quantities and fit diagnostics.
- `round22_dashboard_parameters.csv`: scalar observables and fit diagnostics.
- `round22_held_coefficients.csv`: HELD basis coefficients.
- `round22_held_mean_dispersion.png` / `.dat`: mean-coefficient dispersion.
- `round22_held_summary.json`: provenance, source hash, model and validation summary.

Run from this directory after adding more complete MD steps:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/anaconda3/bin/python3 make_round22_held_dashboard.py --dpi 110
```

The script regenerates derived outputs from the local `md_100steps.out`.
Default animation speed is 5 frames/second. Full histories remain visible
while red cursors track the current frame.
