# Corrected HCP TDEP references (2026-09-28)

All 21 existing nonmagnetic AIMD archives were refitted with corrected HCP reference cells. **No MD simulation was rerun, and no original NPZ positions or forces were changed.** An incorrect ideal starting/reference geometry does not establish that an equilibrated MD trajectory is unusable.

The updated phonon overlay shows every refit, with mean-position reference diagnostics separated into two panels. All TDEP calculations, inputs, outputs, logs, and tooling live outside the repository in `/Users/dajuarez4/Documents/Fe/dataset/hcp`. Corrected results are in `tdep_corrected/`, workflow scripts in `codes/tdep_workflow/`, and previous repository images in `archived_repository_figures/hcp_pre_reference_correction/`. The repository retains only current figures and this report.

## What changed

The old TDEP cell paired 60-degree basal vectors with the fractional basis `(0,0,0), (2/3,1/3,1/2)`, giving Cmcm rather than HCP symmetry. Corrected cells retain the lattice dimensions and use either `(1/3,1/3,1/2)` or `(2/3,2/3,1/2)` for the second atom, choosing the HCP stacking variant closest to the selected trajectory's mean positions. All 21 corrected unit cells identify as **P6_3/mmc (194)** with spglib.

`infile.ssposcar` is rebuilt as the matching 4×4×4, 128-atom ideal supercell. The original atom ordering and forces are retained. Each frame is translated as a whole to remove center-of-mass drift; no atom is independently moved to force an HCP arrangement. The two slightly rounded input lattices are made symmetry-exact at the same nominal a and c, as in the previous converter; this is a small metric approximation rather than a new DFT calculation.

The first 200 frame indices are excluded as an equilibration window. Finite frames in indices 200–399 are used (199 or 200 frames depending on archive). Force conversion remains Ry/bohr to eV/Å. The second-order fit uses the same requested `-rc2 5` cutoff as before. This short-window refit is not a supercell/cutoff/sampling convergence study.

The path is now explicitly read with **`--readpath`**. Merely writing `infile.qpoints_dispersion` did not activate it in the old commands. The path is Γ–M–K–Γ–A–L–H–A | L–M | K–H, with 120 points per segment. For the 60-degree real-space basis:

| Point | Reciprocal fractional coordinates |
| --- | --- |
| Γ | (0, 0, 0) |
| M | (1/2, 0, 0) |
| K | (2/3, 1/3, 0) |
| A | (0, 0, 1/2) |
| L | (1/2, 0, 1/2) |
| H | (2/3, 1/3, 1/2) |

Explicit coordinates avoid the legacy automatic HEX2 corner/edge naming issue. Plots break lines at discontinuous path segments. DOS is recomputed on a 32×32×32 grid. Free-energy output is retained for provenance, but the old HCP thermodynamic plots have not been revalidated by this phonon correction.

## Reference and sampling diagnostics

Three trajectories have small mean offsets from a single ordered HCP reference over the retained window:

| a, c (Å) | Force-fit R² | Γ optical frequencies (THz) | RMS difference between half-window phonons (THz) |
| --- | --- | --- | --- |
| 2.14, 3.40 | 0.895363 | 9.08167, 9.08167, 20.09103 | 0.499 |
| 2.16, 3.42 | 0.882452 | 8.87791, 8.87791, 19.47492 | 0.451 |
| 2.18, 3.40 | 0.893764 | 9.32430, 9.32430, 19.44481 | 0.355 |

All three have nonnegative frequencies on the sampled path, including independent fits of frames 200–299 and 300–399. Maximum between-half differences are 1.491, 1.243, and 0.822 THz, respectively. These differences quantify finite sampling sensitivity; they are not statistical confidence intervals. The shorter halves do not individually pass every mean-position screen because their site averages fluctuate more strongly.

Ten of the 21 refits retain negative frequencies below −0.001 THz somewhere on the sampled path; these are retained rather than clipped or presented as stable.

The other 18 refits are retained in the lower overlay panel. Their larger mean-site offsets mean that a single ideal HCP reference with this atom mapping is less well supported. This does **not** establish that the simulations are wrong, nor does it identify their phases. Positive phonons alone also do not validate the reference or convergence. These runs need further structural/mapping analysis before quantitative HCP interpretation.

The transparent, heuristic screen uses: mean-site RMS ≤0.10 Å, maximum mean-site offset ≤0.20 Å, either half-window mean-site RMS ≤0.15 Å, and RMS change of site means between halves ≤0.20 Å. All per-case values are in `reference_screening.json` and `hcp_reference_validation.json`. This is a conservative mapping diagnostic, not an HCP phase classifier or thermodynamic stability criterion.

## Reproduce

From the external HCP workspace, with NumPy, Matplotlib, and the local TDEP executables available:

```bash
cd /Users/dajuarez4/Documents/Fe/dataset/hcp
python codes/tdep_workflow/recompute_hcp_reference.py \
  --dataset-dir . --tdep-root ../../tdep/build/src
python codes/tdep_workflow/plot_corrected_hcp.py \
  --results tdep_corrected \
  --output ../../IronCoreMD/assets/hcp_phonon_dispersion_overlay.png
```

The recompute command deliberately retains all cases as diagnostic refits. The general NPZ converter rejects a failed reference screen by default; `--allow-unvalidated-hcp` explicitly permits diagnostic runs. `--prepare-only` writes inputs without executing TDEP, and `--run-existing` runs prepared inputs. Results use a separate directory to preserve old calculations.

For each sampling half, run the same recompute command with a separate `--output-root`, `--skip 200` or `--skip 300`, `--max-frames 100`, and `--targets a_2.14_c_3.40_5000K.npz a_2.16_c_3.42_5000K.npz a_2.18_c_3.40_5000K.npz`.

All generated inputs (including positions and forces), numeric results, HDF5 files, and fit logs remain in the external workspace. Per-case SHA-256 provenance identifies the unchanged source NPZ archives. The repository does not contain the TDEP result bundle or workflow scripts.

Validation: all 21 cells have space group 194; all selected forces match the source archives with unchanged ordering and the documented unit conversion; all 21 dispersion tables contain six branches at 1080 q points; source NPZ hashes match between the workspace and repository. Regression tests cover HCP symmetry/coordination, both stacking variants, rigid translation invariance, rejection of the old Cmcm reference, and the reciprocal-space corner geometry.
