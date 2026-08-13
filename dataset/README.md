# IronCoreMD dataset

This directory contains the curated, repository-ready subset of the Fe
workspace. It intentionally includes compressed Quantum ESPRESSO trajectories,
the magnetic runs needed to reproduce current analyses, and compact dataset
summaries. Large scratch directories, wavefunctions, restart files, caches,
and generated TDEP workspaces are not included.

## Layout

- `bcc/non-mag/`: non-magnetic BCC trajectory archives, organized by lattice
  parameter and temperature.
- `bcc/magnetic-collinear/`: the 128-atom collinear magnetic trajectory and
  its ideal reference archive.
- `bcc/magnetic-noncollinear/`: the 54-atom noncollinear trajectory, local-spin
  data, and compact visualization inputs.
- `fcc/non-mag/`: non-magnetic FCC trajectory archives.
- `fcc/legacy_misclassified/`: two historical archives retained for
  provenance; their cells are BCC-like and they are excluded by the default
  FCC selectors.
- `hcp/`: HCP trajectory archives, selected QE inputs, phase-trace output,
  paramagnetic EOS data, and the staged noncollinear HCP campaign.
- `figures/`: compact CSV and visual summaries of phase coverage, dataset
  quality, and ML-scale demonstrations.

Recent curated additions include:

- `hcp/hcp_mag/lessons_learned/`: compact convergence history and operational
  lessons from the 54-atom constrained-paramagnetic HCP campaign at 5000 K;
- `figures/volume_vs_pressure_5000K_bcc_hcp_compare.{csv,png}`: the current
  BCC-HCP pressure-volume comparison at 5000 K.

## Data policy

The `.npz` files are the portable training and analysis inputs. Rebuildable
per-configuration postprocessing directories are deliberately kept outside
Git. The default ML selectors use:

```text
dataset/bcc/non-mag/*.npz
dataset/fcc/non-mag/*.npz
dataset/hcp/*.npz
```

Magnetic trajectories are opt-in because their schemas include additional
magnetization fields and they represent separate physical sampling protocols.
