# Recommended HCP non-collinear workflow

The existing HCP `try1` and `try2` outputs never finish their first SCF cycle. The successful BCC reference differs in one crucial respect: it is non-collinear but does not apply `constrained_magnetization='atomic direction'`. It uses `local-TF`, `mixing_beta=0.01`, `mixing_ndim=20`, and a production threshold of `4.0d-3 Ry`; 65 electronic cycles converge with a median of 11 iterations, producing 64 complete configurations (MD markers 3--66) before the next configuration fails in `cdiaghg`.

## Run order

1. Run `01_hcp_unconstrained_preconverge.in`. It performs a tight (`1.0d-4 Ry`) unconstrained SCF and writes the charge density and wavefunctions to `./hcp_preconverge`.
2. Only if step 1 reports `convergence has been achieved`, run `02_hcp_constrained_md_from_scf.in` from the same directory. It reads the saved density/wavefunctions, uses the BCC-tested production threshold (`4.0d-3 Ry`), and starts with a gentler direction penalty `lambda=0.05 Ry`.
3. Do not raise `lambda` to `0.2 Ry` until several ionic steps converge cleanly. If the first constrained cycle stalls, insert a constrained SCF stage at `lambda=0.02 Ry` before MD.

Example commands (adjust the executable/launcher for the cluster):

```bash
pw.x -in 01_hcp_unconstrained_preconverge.in > 01_hcp_unconstrained_preconverge.out
pw.x -in 02_hcp_constrained_md_from_scf.in > 02_hcp_constrained_md_from_scf.out
```

Both runs must use the same working directory, `prefix`, `outdir`, pseudopotential, MPI layout, and Quantum ESPRESSO build. The `hcp_preconverge` directory must remain intact between stages.

## Why these changes

| Setting | HCP try2 | Recommended production | Successful BCC reference |
|---|---:|---:|---:|
| Preconverged density | no | yes | subsequent steps reuse prior state |
| `mixing_mode` | local-TF | local-TF | local-TF |
| `mixing_beta` | 0.01 | 0.01 | 0.01 |
| `mixing_ndim` | 20 | 20 | 20 |
| `conv_thr` (Ry) | 1e-4 | 4e-3 | 4e-3 |
| Direction constraint | 0.2 Ry immediately | 0.05 Ry after SCF | none |
| Initial eigensolver threshold | 1e-6 | 1e-4 | automatically relaxed |

The staged workflow preserves the intended frozen-direction model in production, but it is still more difficult than the unconstrained BCC calculation. If physical frozen spin directions are not required, removing the constraint entirely is the most BCC-like and most robust option.
