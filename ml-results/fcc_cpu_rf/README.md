# FCC CPU Regression Baseline

This folder contains a CPU-only machine-learning baseline for the FCC non-magnetic Fe dataset at `5000 K`.

The model was trained with [ml_cpu_regression.py](/Users/dajuarez4/Documents/Fe/IronCoreMD/codes/ml_cpu_regression.py), which was designed as a CUDA-free alternative to the original GraphDot workflow in [ml_gpr.py](/Users/dajuarez4/Documents/Fe/IronCoreMD/codes/ml_gpr.py).

## Result

For this FCC-only run, the test-set error is:

- `RMSE = 0.0160045 eV/atom = 16.0 meV/atom`
- `MAE = 0.0073561 eV/atom = 7.36 meV/atom`
- `N_train = 2128`
- `N_test = 1064`

This is a strong baseline for a simple CPU-side model using only geometric and cell-based descriptors.

## Dataset

The run used all eight FCC NPZ archives:

- `2.85_5000K.npz`
- `2.90_5000K.npz`
- `2.95_5000K.npz`
- `3.00_5000K.npz`
- `3.05_5000K.npz`
- `3.10_5000K.npz`
- `3.15_5000K.npz`
- `3.20_5000K.npz`

Each archive contributed `399` valid frames, for a total of `3192` structures.

The train/test split is random with seed `0`, and uses:

\[
N_{\mathrm{test}} = \left\lfloor \frac{N}{3} \right\rfloor,
\qquad
N_{\mathrm{train}} = N - N_{\mathrm{test}}
\]

with `N = 3192`.

## Model

The fitted regressor is a **random forest**:

\[
\hat{y}(\mathbf{x}) = \frac{1}{T}\sum_{t=1}^{T} f_t(\mathbf{x})
\]

where:

- `\mathbf{x}` is the feature vector for one structure,
- `f_t` is the prediction of tree `t`,
- `T = 400` trees in this run.

Each tree is trained on a bootstrap sample of the training set, and split decisions are chosen to reduce the variance of the target energy per atom. In practice, the forest learns a nonlinear map:

\[
\mathbf{x} \mapsto E_{\mathrm{per\,atom}}
\]

without requiring CUDA, graph kernels, or message-passing descriptors.

## Structural Descriptor

For each structure, the feature vector is:

\[
\mathbf{x} =
\Big[
N_{\mathrm{atoms}},
V/N,
N/V,
a,b,c,\alpha,\beta,\gamma,
\{\mu_k\}_{k=1}^{12},
\{\sigma_k\}_{k=1}^{12},
\text{phase one-hot}
\Big]
\]

where:

- `V/N` is the volume per atom,
- `N/V` is the atomic density,
- `a,b,c,\alpha,\beta,\gamma` are the unit-cell parameters,
- `\mu_k` and `\sigma_k` are the mean and standard deviation of the `k`-th nearest-neighbor distance over all atoms in the frame.

More explicitly, if `d_i^{(k)}` is the distance from atom `i` to its `k`-th nearest neighbor, then:

\[
\mu_k = \frac{1}{N_{\mathrm{atoms}}}\sum_{i=1}^{N_{\mathrm{atoms}}} d_i^{(k)}
\]

and

\[
\sigma_k =
\sqrt{
\frac{1}{N_{\mathrm{atoms}}}
\sum_{i=1}^{N_{\mathrm{atoms}}}
\left(d_i^{(k)} - \mu_k\right)^2
}
\]

For this run, `12` neighbor shells were used.

Because this dataset is FCC-only, the phase one-hot encoding does not provide the same benefit it would in a mixed `bcc + fcc + hcp` training set, but it is kept for consistency with the general pipeline.

## Target

The regression target is the **energy per atom**:

\[
E_{\mathrm{per\,atom}} = \frac{E_{\mathrm{total}}}{N_{\mathrm{atoms}}}
\]

where the total energy is read from the QE-derived NPZ archives and converted from Ry to eV.

## Error Metrics

The reported metrics are:

\[
\mathrm{MAE} =
\frac{1}{M}\sum_{j=1}^{M}
\left| y_j - \hat{y}_j \right|
\]

\[
\mathrm{RMSE} =
\sqrt{
\frac{1}{M}\sum_{j=1}^{M}
\left(y_j - \hat{y}_j\right)^2
}
\]

where `M` is the number of test structures.

For this run:

\[
\mathrm{RMSE} = 1.60045\times10^{-2}\ \mathrm{eV/atom}
\]

which is `16.0 meV/atom`.

## Most Important Features

The learned importance ranking is dominated by local geometric features, especially the mid-range nearest-neighbor shells:

1. `nn5_mean`
2. `nn4_mean`
3. `nn6_mean`
4. `nn7_mean`
5. `nn3_mean`
6. `a_A`
7. `c_A`
8. `b_A`
9. `volume_A3_atom`
10. `density_atoms_A3`

This is physically reasonable: the energy is strongly controlled by the local environment and effective lattice spacing.

## How To Reproduce

From the `IronCoreMD` repository root:

```bash
python3 codes/ml_cpu_regression.py \
  --phase fcc \
  --run-name fcc_cpu_rf \
  --output-root ml-results
```

The exact hyperparameters used here were:

- model: `random_forest`
- `n_estimators = 400`
- `n_neighbors = 12`
- `seed = 0`
- `min_samples_leaf = 1`

## Files In This Folder

- [fcc_cpu_rf_dataset_preview.csv](/Users/dajuarez4/Documents/Fe/IronCoreMD/ml-results/fcc_cpu_rf/fcc_cpu_rf_dataset_preview.csv)
- [fcc_cpu_rf_train_split.csv](/Users/dajuarez4/Documents/Fe/IronCoreMD/ml-results/fcc_cpu_rf/fcc_cpu_rf_train_split.csv)
- [fcc_cpu_rf_test_split.csv](/Users/dajuarez4/Documents/Fe/IronCoreMD/ml-results/fcc_cpu_rf/fcc_cpu_rf_test_split.csv)
- [fcc_cpu_rf_dataset_preview.png](/Users/dajuarez4/Documents/Fe/IronCoreMD/ml-results/fcc_cpu_rf/fcc_cpu_rf_dataset_preview.png)
- [fcc_cpu_rf_parity_plot.png](/Users/dajuarez4/Documents/Fe/IronCoreMD/ml-results/fcc_cpu_rf/fcc_cpu_rf_parity_plot.png)
- [fcc_cpu_rf_feature_importance.csv](/Users/dajuarez4/Documents/Fe/IronCoreMD/ml-results/fcc_cpu_rf/fcc_cpu_rf_feature_importance.csv)
- [fcc_cpu_rf_feature_importance.png](/Users/dajuarez4/Documents/Fe/IronCoreMD/ml-results/fcc_cpu_rf/fcc_cpu_rf_feature_importance.png)
- [fcc_cpu_rf_test_predictions.csv](/Users/dajuarez4/Documents/Fe/IronCoreMD/ml-results/fcc_cpu_rf/fcc_cpu_rf_test_predictions.csv)
- [fcc_cpu_rf_random_forest_EperAtom_2128.pkl](/Users/dajuarez4/Documents/Fe/IronCoreMD/ml-results/fcc_cpu_rf/fcc_cpu_rf_random_forest_EperAtom_2128.pkl)

## Interpretation

This model is not meant to replace the original graph-kernel GPR on physical grounds. Instead, it provides:

- a practical CPU-only baseline,
- a fast sanity check for the dataset,
- a useful regression benchmark when CUDA is unavailable,
- and an interpretable connection between energy and local structural descriptors.

The `16 meV/atom` RMSE shows that a relatively simple geometric descriptor plus ensemble regression already captures much of the variation in the FCC dataset.
