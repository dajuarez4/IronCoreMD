# Phase-Aware TDEP Workflow

This directory packages the reusable TDEP postprocessing workflow for the non-magnetic iron datasets in the repository.

The workflow starts from compressed QE AIMD archives such as `dataset/bcc/*.npz` or `dataset/fcc/*.npz`, builds `tdep_*` folders, runs the harmonic TDEP fit, and regenerates the thermodynamic and phonon plots used in the repository.

## Included Scripts

- `npz_to_tdep.py`
  Generic NPZ-to-TDEP converter for the supported phases:
  - `bcc`
  - `fcc`
  It writes:
  - `infile.ucposcar`
  - `infile.ssposcar`
  - `infile.positions`
  - `infile.forces`
  - `infile.stat`
  - `infile.meta`
  - `infile.qpoints_dispersion`
  - `source_npz.txt`

- `npz_to_tdep_bcc.py`
  Backward-compatible wrapper for `npz_to_tdep.py --phase bcc`.

- `npz_to_tdep_fcc.py`
  Convenience wrapper for `npz_to_tdep.py --phase fcc`.

- `run_harmonic_tdep.py`
  Generic end-to-end driver for:
  - rebuilding selected `tdep_*` folders from NPZ archives,
  - running `extract_forceconstants`,
  - running phonon dispersion and DOS/free-energy calculations,
  - refreshing the single-temperature plots,
  - and optionally refreshing the multi-temperature comparison figures.

- `run_bcc_harmonic_tdep.py`
  Backward-compatible wrapper for `run_harmonic_tdep.py --phase bcc`.

- `run_fcc_harmonic_tdep.py`
  Convenience wrapper for `run_harmonic_tdep.py --phase fcc`.

- `summarize_free_energy.py`
  Prints a CSV-style summary of `outfile.free_energy` plus `outfile.U0`.

- `plot_free_energy_vs_volume.py`
  Writes:
  - free energy vs volume,
  - relative free energy vs volume,
  - free energy vs lattice parameter,
  - relative free energy vs lattice parameter,
  - and the thermodynamic CSV used by the pressure plot.

- `plot_volume_vs_pressure.py`
  Fits a Birch-Murnaghan EOS to the TDEP free energies and to the AIMD mean pressures, then writes:
  - a comparison pressure-volume figure,
  - an EOS-only pressure-volume figure,
  - and a CSV with the derived pressure values.

- `plot_combined_dispersion.py`
  Writes a combined phonon-dispersion and total-DOS overlay for one temperature.

- `plot_temperature_comparison.py`
  Overlays free-energy and pressure-volume curves across multiple temperatures.

- `tdep_common.py`
  Shared helper functions for folder discovery, duplicate handling, thermodynamic parsing, and default output naming.

## Data Layout

The scripts assume the repository layout:

```text
IronCoreMD/
├── codes/
│   └── tdep_workflow/
└── dataset/
    ├── bcc/
    │   ├── 2.29_5000K.npz
    │   ├── 2.52_5000-new.npz
    │   ├── 2.51_5500K.npz
    │   ├── ...
    │   └── tdep_2.29_5000K/
    └── fcc/
        ├── 3.00_5000K.npz
        ├── 3.05_5000K.npz
        ├── ...
        └── tdep_3.00_5000K/
```

The current defaults target `dataset/<phase>` relative to the repository root, and every script accepts both `--phase` and `--dataset-dir`.

## Requirements

Python-side requirements:

- `numpy`
- `scipy`
- `matplotlib`

External requirements:

- a built TDEP checkout with:
  - `extract_forceconstants`
  - `phonon_dispersion_relations`

The driver script auto-detects TDEP in either:

- `IronCoreMD/tdep/build/src`
- or a sibling checkout at `../tdep/build/src`

You can also pass `--tdep-root` explicitly.

## Quick Start

Rebuild all `5000 K` BCC TDEP folders from the NPZ archives:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD
python codes/tdep_workflow/npz_to_tdep_bcc.py --temperature-K 5000
```

Run the full harmonic TDEP workflow for one BCC temperature:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD
python codes/tdep_workflow/run_bcc_harmonic_tdep.py --temperature-label 5000
```

Run the FCC workflow with the generic phase-aware entrypoint:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD
python codes/tdep_workflow/run_harmonic_tdep.py --phase fcc --temperature-label 5000
```

or with the FCC convenience wrapper:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD
python codes/tdep_workflow/run_fcc_harmonic_tdep.py --temperature-label 5000
```

## Targeted Reruns

Rebuild and rerun only the updated `5500 K` high-volume subset:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD
python codes/tdep_workflow/run_bcc_harmonic_tdep.py \
  --temperature-label 5500 \
  2.51_5500K 2.52_5500K 2.53_5500K 2.54_5500K 2.55_5500K
```

Refresh only the plots after manually editing or rerunning existing `tdep_*` folders:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD
python codes/tdep_workflow/run_bcc_harmonic_tdep.py \
  --temperature-label 5500 \
  --no-convert \
  --no-tdep \
  --no-comparison-plots
```

## Plot and CSV Outputs

For `5000 K`, the free-energy and dispersion scripts keep the existing dataset-local naming:

- `free_energy_vs_volume.csv`
- `free_energy_vs_volume.png`
- `relative_free_energy_vs_volume.png`
- `free_energy_vs_lattice.png`
- `relative_free_energy_vs_lattice.png`
- `volume_vs_pressure_5000K_bcc.csv`
- `volume_vs_pressure_5000K_bcc.png`
- `volume_vs_pressure_5000K_bcc_eos_std.png`
- `phonon_dispersion_overlay.png`

For other phases, the pressure outputs use the phase suffix, for example:

- `volume_vs_pressure_5000K_fcc.csv`
- `volume_vs_pressure_5000K_fcc.png`
- `volume_vs_pressure_5000K_fcc_eos_std.png`

For other temperatures, the temperature is included in the file name, for example:

- `free_energy_vs_volume_5500K.csv`
- `volume_vs_pressure_5500K_bcc.csv`
- `phonon_dispersion_overlay_5500K.png`

The multi-temperature comparison script writes:

- `free_energy_vs_volume_4500K_5000K_5500K.png`
- `volume_vs_pressure_4500K_5000K_5500K_bcc.png`

## Duplicate and Unstable Point Handling

Two workflow rules are built into the discovery helpers:

1. Duplicate lattice points prefer a suffixed replacement folder over the bare folder.
   Example:
   `tdep_2.52_5000-new` is preferred over `tdep_2.52_5000K` when both exist.

2. Thermodynamic plots skip dynamically unstable TDEP points automatically.
   A folder is excluded from the thermodynamic plots when `outfile.free_energy` is non-finite or clearly pathological, such as the large positive values produced by imaginary modes.

## Recommended Commit Scope

If you want to keep repository history clean, treat this directory as the reusable workflow layer and keep large generated `tdep_*` folders or refreshed plot products in separate commits.
