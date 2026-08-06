# BCC Fe noncollinear AIMD at 4000 K

This directory contains the 54-atom, 3x3x3 noncollinear BCC Quantum ESPRESSO trajectory in `fe1.out`.

## Processed trajectory

`simulation.npz` contains 65 position records. Of these, 64 configurations have complete forces and atom-resolved local magnetization vectors:

- complete MD iterations: 3--66;
- trajectory interval: 0.003--0.066 ps;
- analyzed time span: 0.063 ps;
- atoms per frame: 54.

MD iteration 67 has positions but is excluded because its SCF calculation terminated with `cdiaghg: problems computing cholesky` before a complete force and local-spin evaluation.

The exact parsing status is recorded in `md_processing_summary.json`.

## HELD phonons

The synchronized 64-frame HELD diagnostic is stored in:

`../../../HELD/examples/iron/bcc_noncollinear_4000K_3x3x3`

The mean BCC dispersion spans approximately 0--15.08 THz along `Gamma-H-N-Gamma-P-H`, with no negative points on the sampled mean path. These consecutive configurations are strongly correlated, so the result remains a trajectory diagnostic rather than a statistically converged finite-temperature phonon spectrum.

## Dashboard

The atom motion, evolved local-spin arrows, global magnetization, temperature, displacement history, and step-resolved HELD phonons are synchronized in:

`../../../HELD/local_visualizations/bcc_noncollinear_4000K_dashboard/bcc_noncollinear_md_dashboard.gif`

The corresponding final-frame image is:

`../../../HELD/local_visualizations/bcc_noncollinear_4000K_dashboard/bcc_noncollinear_md_dashboard_final.png`

## Rebuild

From the workspace root:

```bash
/opt/anaconda3/bin/python IronCoreMD/codes/data_compress.py \
  dataset/bcc/magnetic-noncollinear \
  --output-dir dataset/bcc/magnetic-noncollinear
```

From the `HELD` directory:

```bash
/opt/anaconda3/bin/python examples/iron/bcc_noncollinear_4000K_3x3x3/run_case.py
```

From the workspace root:

```bash
/opt/anaconda3/bin/python \
  HELD/local_visualizations/bcc_noncollinear_4000K_dashboard/make_dashboard.py
```
