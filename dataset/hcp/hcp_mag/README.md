# Robust noncollinear HCP Fe MD input

`fe.in` is based on the supplied 128-atom HCP `a=2.14`, `c=3.42`, 5000 K
noncollinear input with velocities. It keeps all atomic directions and the
target constraint strength `lambda=0.2`.

The previous `fe1.out` is strongly unstable: at `mixing_beta=0.10`, the SCF
accuracy remains of order 1--13 Ry after more than 100 iterations and the
absolute magnetization oscillates strongly. The revised input uses:

- `mixing_beta=0.01`
- `mixing_mode='local-TF'`
- `mixing_ndim=20`
- conjugate-gradient diagonalization with `diago_thr_init=1e-6`
- `diago_cg_maxiter=200`
- initial normalized moments reduced from `0.15` to `0.05`

If the first ionic step still oscillates, preconverge the same structure with
`calculation='scf'` and `lambda=0.05`, then restart the production MD at
`lambda=0.2`. Lowering the production constraint itself changes the intended
magnetic model and should not be done silently.
