# Noncollinear BCC Fe AIMD

This directory contains the processed 3x3x3 noncollinear BCC Fe trajectory at a nominal 4000 K.

## Processed data

- `fe1.out`: Quantum ESPRESSO output used as the source.
- `simulation.npz`: parsed trajectory with positions, forces, thermodynamic observables, converged global magnetization, and atom-resolved local magnetization vectors for each complete MD frame.
- `md_processing_summary.json`: frame counts and sampled temperature, pressure, and time ranges.

The current output contains 23 position records and 22 complete force frames (MD iterations 3--24). The final position record is incomplete because `fe1.out` was still running when it was processed, so it is excluded from the plots, GIFs, and HELD fit.

## Figures and animation

- `initial_SQS_spin_directions_3d.png`: atom-resolved starting-spin directions read from the QE input.
- `initial_SQS_spin_vectors.csv`: Cartesian unit vectors for those starting directions.
- `noncollinear_md_spin_displacement_results.png`: global magnetization, lattice displacement, temperature, and pressure evolution.
- `bcc_noncollinear_3x3x3_4000K_no_vito.gif`: 22-frame trajectory animation generated directly with Matplotlib, without OVITO.

Quantum ESPRESSO reports both the converged total magnetization vector and 54 atom-resolved local magnetization vectors integrated inside atomic spheres for every complete MD step. The initial SQS figure remains useful as the prescribed starting texture, while the HELD dashboard animates the evolved local vectors.

## HELD analysis

The diagnostic HELD case is in `HELD/examples/iron/bcc_noncollinear_4000K_3x3x3`. It uses an ideal 3x3x3 BCC reference lattice rather than the initially displaced coordinates.

The 22 complete frames span only 0.021 ps and are consecutive, strongly correlated samples. In addition, the fifth-shell cutoff exceeds half the 7.65-Angstrom box length. The resulting phonons are useful as an early trajectory diagnostic, not as a statistically converged finite-temperature spectrum.
