# SCF and magnetization analysis

Both output files stop during the **first electronic SCF cycle**. No ionic MD step was completed, so every plot uses SCF iteration—not MD time step.

| Quantity | try1 | try2 |
|---|---:|---:|
| Complete printed SCF iterations | 115 | 86 |
| Last estimated accuracy (Ry) | 1.641 | 4.227 |
| Minimum estimated accuracy (Ry) | 1.18 | 1.824 |
| Last absolute M (muB/cell) | 45.560 | 48.670 |
| Last absolute M (muB/atom) | 0.356 | 0.380 |
| Late absolute-M mean (muB/atom) | 0.568 | 0.421 |
| Late absolute-M std (muB/atom) | 0.240 | 0.067 |
| Late absolute-M lag-1 correlation | 0.561 | 0.874 |
| Apparent absolute-M period (iterations) | 95.00 | 33.00 |

## Interpretation

- Neither run is converged: the requested threshold is 1e-4 Ry, while the final printed residuals remain many orders of magnitude larger. The ionic integrator therefore never starts.
- The huge early absolute magnetization is an initialization transient. It collapses quickly, then fluctuates because charge/spin mixing is cycling among substantially different electronic states.
- The component sign changes are not physical 5000 K spin dynamics: they occur inside one unconverged electronic minimization. The `atomic direction` penalty constrains directions but not fixed moment magnitudes, and the 128 independently typed atoms make the SCF problem unusually stiff.
- try1's larger beta updates the density ten times more aggressively and shows stronger overshoot. try2 damps the update, but its residual stalls instead of approaching 1e-4 Ry. Thus the smaller beta reduces aggression but does not cure the underlying instability.
- Any Fourier-derived period is only descriptive, not evidence of a stable physical oscillation: the traces are nonstationary, short, and unconverged.

## Practical next checks

Converge a static SCF at the initial geometry before MD; restart from its saved density/wavefunctions; test a staged constraint lambda and starting moment; and compare `plain`/`local-TF` mixing with smaller beta and larger history. Do not use forces, energies, or magnetization from these two outputs as MD data yet.
