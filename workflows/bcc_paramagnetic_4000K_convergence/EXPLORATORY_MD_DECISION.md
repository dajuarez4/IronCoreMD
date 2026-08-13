# Accepted exploratory-MD settings for 16-atom BCC Fe at 4000 K

## Decision

Use `conv_thr=4e-3 Ry` for exploratory constrained-magnetic MD. Do not label
these trajectories as production-quality AIMD, and do not use their forces for
quantitative thermodynamic, vibrational, or transport results.

Accepted electronic and magnetic settings:

```text
constrained_magnetization='atomic direction'
lambda=0.020000

startingpot='file'
startingwfc='file'
conv_thr=4.0d-3
electron_maxstep=500
scf_must_converge=.true.
mixing_mode='local-TF'
mixing_beta=0.010d0
mixing_ndim=20
diagonalization='cg'
diago_thr_init=1.0d-3
diago_cg_maxiter=200

pot_extrapolation='atomic'
wfc_extrapolation='none'
```

The MD must start from an untouched copy of the converged `lambda=0.020` SCF
state. Do not seed it from an MD state that has already developed a large net
magnetization.

## Evidence

The five-step `conv_thr=4e-3 Ry` pilot completed with `JOB DONE`, no Cholesky
failure, and one recoverable `c_bands` warning. Its five converged electronic
cycles required 7, 26, 8, 20, and 10 SCF iterations. Temperature ended at
3979 K.

The converged total-magnetization magnitude increased as follows:

```text
0.783 -> 1.171 -> 1.120 -> 1.677 -> 2.185 mu_B/cell
```

Absolute magnetization stayed between 5.45 and 5.79 mu_B/cell. This is much
better than the rejected weak-constraint trajectory, whose net moment grew to
about 13.3 mu_B/cell, but the upward drift must still be monitored.

The ratio of QE's total SCF force correction to total force was:

```text
32.75, 2.42, 2.00, 1.60, 1.04
```

QE warned after every force evaluation that the SCF correction was large.
Consequently, numerical completion does not imply quantitatively converged
forces.

The tighter `conv_thr=1e-3 Ry` test improved the later force-correction ratios
to approximately 0.96 and 0.93, but failed with
`cdiaghg: problems computing cholesky` at MD iteration 4. It is therefore not
the accepted stable workflow.

## Permitted use

The accepted `4e-3 Ry` trajectory may be used for:

- numerical and workflow testing;
- qualitative inspection of finite-temperature structures;
- diagnosing magnetic drift and MD stability;
- generating candidates that will later be recomputed with tighter settings.

It must not be used without further validation for:

- quantitative forces or force constants;
- phonons or TDEP fitting;
- diffusion coefficients or transport properties;
- equations of state, phase boundaries, or free-energy differences;
- final BCC-versus-HCP energetic comparisons.

## Monitoring and rejection rules

Run short chunks first and inspect only magnetization values from converged SCF
cycles. Temporary values inside an unconverged SCF loop are not trajectory
observables.

Stop and review the trajectory if any of the following occurs:

- Cholesky or unrecovered `c_bands` failure;
- an SCF cycle reaches `electron_maxstep`;
- persistent monotonic growth of the converged net magnetization;
- net magnetization approaches the previously rejected partially ferromagnetic
  range (roughly 5 mu_B/cell or larger for this 16-atom cell);
- temperature control becomes unstable;
- the scientific task requires quantitatively reliable forces.

Archive every attempt in its own directory. Always keep the untouched
preconditioned state separate from the evolving MD `outdir`.
