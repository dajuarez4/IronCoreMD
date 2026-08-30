# Round 14: unconstrained-bootstrap 32-atom BCC noncollinear MD400

Round 14 tests a different solution to the Round-12 SCF failure. It first
converges the noncollinear electronic state without an atomic magnetic
constraint, then enables the patched atomic constraint and ramps its strength.

The geometry, thermal displacements, velocities, eight magnetic directions,
and spatially mixed eight-species assignment are identical to Round 13. This
keeps the Round-13/14 comparison focused on the SCF initialization route.

## Why Rounds 13 and 14 were created

Round 12 stopped in its first SCF stage after all 500 allowed electronic
iterations. Its target was `conv_thr=4e-3 Ry`, but the best estimated accuracy
was about `4.835 Ry` and the final value was about `4.993 Ry`. The residual
settled near `5--6 Ry` for hundreds of iterations while the energy and magnetic
state continued to oscillate. There was no eigensolver, Cholesky, MPI, or
patched-QE runtime error; the failure was in charge/spin-density convergence.

The Round-12 negative-density diagnostic grew to approximately `3.83`, whereas
the successfully converged eight-atom Round 11 calculation ended near
`1e-6`. Round 12 printed `JOB DONE` because QE exited normally after the failed
SCF; this does not mean the SCF converged. Its results archive correctly marks
stage 1 as `FAILED` and all later stages and MD as `NOT_STARTED`.

Round 12 also assigned Fe01--Fe04 to one BCC sublattice and their antipodal
directions Fe05--Fe08 to the other. The full target-vector sum was zero, but
the spatial pattern was ordered and sublattice polarized rather than a mixed
disordered-local-moment seed.

Round 13 corrected that spatial assignment and introduced a weak constrained
bootstrap. Round 14 retains the corrected Round-13 structure and magnetic
layout but adds an unconstrained first stage. The intended comparison is:

- Round 13: weak atomic constraint from the first SCF;
- Round 14: converge charge and spin density without the constraint, then
  activate and ramp the patched atomic constraint.

If Round 13 fails while Round 14 converges, the initial constrained SCF is the
likely source of the instability. If both fail similarly, the next suspects
are the 32-atom initialization, density mixing, or the repeated-species
constraint representation rather than merely the value of the first lambda.

## QE species limit

Every input uses `nat=32` and `ntyp=8`. The artificial species Fe01--Fe08 are
each assigned to four atoms, with two copies on each BCC sublattice. This stays
below the private QE 7.5 executable's ten-species limit.

The magnetic correction does not create one species per atom. Every generated
Round-14 input has:

- `nat=32`;
- `ntyp=8`;
- eight `ATOMIC_SPECIES` cards, Fe01 through Fe08;
- four atoms per artificial species;
- two copies of every species on each 16-site BCC sublattice.

Consequently, each sublattice and the full 32-atom cell independently have a
target magnetic-vector sum equal to zero to floating-point precision.

## SCF path

1. `bootstrap_unconstrained_2e-2`: `constrained_magnetization='none'`,
   `degauss=0.08 Ry`, `conv_thr=2e-2 Ry`, and `mixing_beta=0.05`;
2. `bootstrap_l0p001_1e-2`: atomic constraint at lambda 0.001 Ry,
   `degauss=0.05 Ry`, and `conv_thr=1e-2 Ry`;
3. `stage1_l0p005_4e-3`: lambda 0.005 Ry and `conv_thr=4e-3 Ry`;
4. `stage2_l0p020_1e-3`: lambda 0.020 Ry, production smearing, and
   `conv_thr=1e-3 Ry`;
5. `stage3_l0p100_3e-4`: lambda 0.100 Ry and `conv_thr=3e-4 Ry`;
6. 400 MD steps at lambda 0.100 Ry.

All stages use `mixing_mode='plain'` and a 40-vector mixing history. The first
stage uses `mixing_beta=0.05`; constrained stages and MD use 0.01.

The job stops at the first failed SCF and starts in its own `qe_tmp`. Do not
copy or reuse Round 12's failed electronic state.

## Local validation completed

The generated package was checked locally before submission:

- the generator and collector compile with Python;
- all shell and Slurm scripts pass shell syntax validation;
- all six QE inputs contain `nat=32`, `ntyp=8`, and exactly eight species;
- the unconstrained input contains no lambda and uses
  `constrained_magnetization='none'`;
- every subsequent SCF and MD input uses the atomic constraint at the intended
  lambda;
- Round-13 and Round-14 atomic labels and coordinates compare exactly equal;
- both sublattices contain exactly two copies of Fe01--Fe08;
- the lightweight collector correctly reports `NOT_STARTED` before outputs
  exist and refuses `--require-complete` on an incomplete trajectory.

These checks validate package construction only. The actual electronic
convergence must be tested with the private patched QE executable on
Lonestar6.

## Round 11 phonon feasibility

Round 11 successfully completed 400 MD steps and is useful as an electronic
and magnetic stability demonstration. A phonon calculation can technically be
attempted from it, but the result should be labeled an exploratory trajectory
diagnostic, not a converged 4000 K TDEP phonon spectrum.

Important limitations are:

- the cell is only eight atoms and `2x2x1`, so it is small and especially
  short along one lattice direction for fitting spatial force constants;
- 400 steps at about 0.5 fs cover only about 200 fs, and consecutive frames
  are strongly correlated;
- the trajectory mean temperature is about `2945 K`, with about `2977 K` over
  the final 100 steps, rather than the requested 4000 K;
- all 400 steps print `SCF correction compared to forces is large`;
- the median printed SCF-correction/total-force ratio is about `1.01`, and the
  ratio exceeds one in 203 of 400 steps;
- the existing lightweight archive contains 400 position blocks but 2,800
  force blocks because its parser retained the total-force table plus six
  component tables for every MD step. It is therefore not directly suitable
  for the standard NPZ-to-TDEP converter.

The full QE output does contain 400 final total-force tables. A corrected
parser could construct synchronized samples using the input position with the
first force table and output positions 1--399 with force tables 2--400. Even
then, the loose MD force convergence remains the dominant accuracy problem.

The preferred phonon route is to select decorrelated, post-equilibration
snapshots and recompute their forces using much tighter static SCF thresholds,
then fit only the synchronized total forces. This would still be a small-cell
diagnostic. Quantitative finite-temperature phonons should instead use a
successful larger, three-dimensionally extended cell and a longer equilibrated
trajectory with force-quality criteria enforced during or after MD.

## Copy and submit

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear
scp -r round14_32atom_4x2x2_qe75_patched_unconstrained_bootstrap_md400_lonestar6 \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/
```

On Lonestar6:

```bash
cd /work/11381/dajuarez4/ls6/round14_32atom_4x2x2_qe75_patched_unconstrained_bootstrap_md400_lonestar6
mkdir -p logs
sbatch run_round14_lonestar6.sbatch
```

Monitor and collect with:

```bash
bash check_round14.sh
bash plot_round14.sh --live 30
python3 collect_round14_npz.py --require-complete
```
