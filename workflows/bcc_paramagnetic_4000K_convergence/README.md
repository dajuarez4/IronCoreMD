# BCC constrained-paramagnetic Fe at 4000 K

This curated workflow records the convergence study for a 16-atom BCC Fe
cell with noncollinear atomic-direction constraints. It contains only
repository-ready inputs, scripts, and decisions; raw Quantum ESPRESSO output,
Slurm logs, wavefunctions, and restart directories remain outside Git.

The accepted numerical choice is `conv_thr=4e-3 Ry`, `lambda=0.020 Ry`,
local-TF mixing, CG diagonalization, atomic potential extrapolation, and no
wavefunction extrapolation. It is accepted for **exploratory MD only** because
QE's SCF force correction remained comparable to the total force.

- [`EXPLORATORY_MD_DECISION.md`](EXPLORATORY_MD_DECISION.md) contains the
  measured evidence, limitations, and rejection rules.
- [`exploratory_10steps/`](exploratory_10steps/) is a guarded ten-step Jakar
  workflow that starts from an external clean preconditioned state.

Do not use this workflow for quantitative TDEP, phonon, diffusion,
free-energy, or BCC-versus-HCP results without tighter force validation.
