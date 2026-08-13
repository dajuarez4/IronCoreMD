# Lessons learned: constrained-paramagnetic HCP Fe at 5000 K

## System and intent

- Structure: HCP Fe, 3x3x3 supercell, 54 atoms.
- Primitive parameters: a = 2.186731 A, c = 3.553438 A, c/a = 1.625.
- Electronic temperature: 5000 K using Fermi-Dirac smearing and
  degauss = 0.031668 Ry.
- Magnetic model: noncollinear, one Fe type per site, random zero-net starting
  directions, and `constrained_magnetization='atomic direction'`.
- Current k-point sampling is an explicit Gamma point. It is suitable for
  workflow screening, not final Brillouin-zone convergence.

## What the calculations demonstrated

1. A direct weak-to-strong penalty jump is unstable. The original
   `lambda=0.02 -> 0.10 Ry` transition failed after 500 iterations even though
   the charge density and wavefunctions were read correctly.
2. A gentler ramp works. Round 2 converged at lambda 0.03, 0.05, and 0.075 Ry.
3. The remaining `0.075 -> 0.10 Ry` jump was still too large. It approached
   1.8e-3 Ry, then diverged and ended at 2.78 Ry with a large magnetic spike.
4. A QE `JOB DONE` marker means the executable ended normally; it does not mean
   that the SCF converged. Always check for `convergence has been achieved`.
5. A single shared `outdir` is unsafe for a sequential convergence search. A
   failed later stage overwrites the last good charge density and wavefunctions.
6. No HCP stage showed the Davidson Cholesky failure seen in the BCC tests, but
   the unstable lambda=0.10 stages produced many unconverged-band warnings.
7. The SCF threshold is extensive. The present 1e-3 Ry target applies to the
   complete 54-atom cell. Tightening to 1e-4 Ry should happen only after the
   lambda ramp is stable.
8. The QE `atomic direction` penalty primarily constrains the polar direction.
   It should not be interpreted as a full fixed three-component local moment.

## Workflow rules going forward

- Change only one hard control at a time.
- Increase lambda in small increments and stop immediately after a failed SCF.
- Store each converged lambda in an independent state directory.
- Preserve charge density and distributed wavefunctions together.
- Do not start MD until the final SCF stage is reproducibly converged.
- For MD, avoid second-order wavefunction extrapolation until it has been shown
  stable; the BCC reference failed in diagonalization immediately after using it.
- Validate the selected workflow with a denser k-point mesh before comparing
  final energies, pressures, or phase stability.

## Current next experiment

Round 3 first recovers lambda=0.075 Ry from the failed Round 2 active state,
then attempts 0.080, 0.085, 0.090, 0.095, and 0.100 Ry. Each stage has its own
checkpoint. MD remains disabled.

Separately, `md_pilot_lambda0p075/` provides a conservative 50-step MD test at
the last demonstrated stable penalty. It reconstructs the proven 0.03 -> 0.05
-> 0.075 Ry SCF path, copies the final state before MD, uses `conv_thr=1e-3 Ry`
and `wfc_extrapolation='none'`, and includes a CSV result extractor. This pilot
must succeed before creating a long trajectory.
