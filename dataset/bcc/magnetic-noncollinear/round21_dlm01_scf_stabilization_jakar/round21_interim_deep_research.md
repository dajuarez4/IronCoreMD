# Round 21 interim decision

Date: 2026-09-03

## Recommendation

Let `conv1e-4_beta1e-2` finish, but do not extend or select any Round-21 case yet.
The two completed failures reject `conv_thr=3e-4` for this difficult seed, and
`mixing_beta=0.005` is substantially worse than 0.010. If the active case is
not both reliable and affordable, the next same-seed pilot should use:

```text
conv_thr            = 1.0e-4 Ry
mixing_beta         = 0.010
mixing_mode         = local-TF
lambda              = 0.05 Ry
scf_must_converge   = .true.
```

QE 7.5 emits the force-correction warning when the SCF-correction norm exceeds
10% of the total-force norm, and its source explicitly recommends reducing
`conv_thr`. This explains why the active `1e-4` case has zero warnings so far.
[QE 7.5 force source](https://gitlab.com/QEF/q-e/-/raw/qe-7.5/PW/src/forces.f90)

Tighter convergence does not remove the underlying constraint stiffness. QE's
official input reference says to reduce `lambda` when constrained-magnetization
SCF does not converge. The `atomic` option penalizes the full local-moment vector,
so it constrains both magnitude and direction.
[QE 7.5 input reference](https://www.quantum-espresso.org/Doc/INPUT_PW.html)

Round 20 maintained exceptionally tight directions at `lambda=0.1`: mean angular
errors were 0.66–0.85°, 95th percentiles 1.38–1.78°, and the observed maximum
was 2.62°. A lambda=0.05 test therefore has measurable room to improve SCF
stability without automatically sacrificing the DLM texture.

Do not switch to QE's `atomic direction` mode: its implementation constrains
only the polar angle relative to z, not each arbitrary three-dimensional target
including azimuth. [QE input source](https://gitlab.com/QEF/q-e/blob/master/PW/src/input.f90)

## Acceptance gate for the lambda=0.05 pilot

- 100/100 MD steps and no unconverged SCF cycle.
- Force-correction warnings on no more than 5% of steps.
- Mean angular error no more than 3° and 95th percentile no more than 7°.
- No systematic local-moment collapse.
- SCF iterations per step low enough to support the intended production length.

The angular and warning gates are project-specific engineering criteria inferred
from Round 20, not universal literature thresholds. DLM-MD studies establish the
importance of sampling coupled magnetic and vibrational disorder, while downstream
TDEP fitting depends directly on reliable matched forces and displacements.
[Alling et al., PRB 93, 224411 (2016)](https://d-nb.info/1229919929/34),
[Mozafari et al., PRB 94, 054111 (2016)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.054111)

## Current evidence limitation

This is an interim decision: the `1e-4` trajectory was only at step 15. Its high
temperature is not yet interpretable as a thermostat failure, and its final SCF
cost and stability remain unknown.
