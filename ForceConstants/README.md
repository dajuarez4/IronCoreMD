# BCC force constants through fifth-nearest neighbors

`bcc_held_5nn.ipynb` discovers the non-magnetic BCC Fe trajectory NPZ files
under `IronCoreMD/dataset/bcc/non-mag`, fits each compatible trajectory with
the repository's unified HELD implementation using five neighbor shells, and
writes compressed NPZ results under `results/`.

The notebook uses `held_bcc_5nn.py` so the calculation is reproducible outside
Jupyter and testable without executing every expensive fit interactively.

## Stored quantities

Each case NPZ contains:

- `fc_labels = [alpha_0, alpha_1, beta_1, ..., alpha_5, beta_5]`;
- `fc_mean` and `fc_per_md_step` in `eV/angstrom^2`;
- the five HELD shell distances;
- complete HELD basis labels, mean coefficients, uncertainties, and per-frame
  coefficients;
- mean and per-frame onsite/offsite force-constant tensors;
- source path, selected frame indices, MD step IDs, and JSON metadata.

These are the 14 traditional monoatomic-BCC Born-von Karman components through
5NN: one onsite component and 13 interatomic components. Only 13 are
independent because HELD derives `alpha_0` from the interatomic tensors using
the acoustic sum rule. Symmetry-equivalent tensors are rotated to the canonical
BCC bond orientation before averaging. The full HELD tensors are retained.

The workflow is intentionally restricted to elemental Fe files in `non-mag`.
Magnetic BCC trajectories, B2 alloys, prior HELD caches, and unrelated NPZ
files are outside this notebook's input scope.

`results/bcc_held_5nn_index.npz` provides one row per attempted trajectory,
including status/error fields and the 14 mean BvK values. Files that
are not valid HELD trajectories are recorded separately in the same index.

## Environment

Use the scientific Python installation that provides NumPy, ASE, spglib, and
the dependencies listed in `../HELD/requirements.txt`:

```bash
/opt/anaconda3/bin/jupyter notebook bcc_held_5nn.ipynb
```

The full dataset contains many 128-atom, 400-frame trajectories and can take a
long time. Run the notebook's smoke test first. Existing results are cached;
set `overwrite=True` only when a deliberate refit is required.
