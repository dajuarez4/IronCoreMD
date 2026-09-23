# IronCoreMD

First-principles iron datasets for high-pressure, high-temperature studies: BCC, FCC and HCP trajectories from Quantum ESPRESSO, finite-temperature phonons from TDEP/HELD, and archives for machine-learning workflows.

## Results

### Dataset overview

<p align="center">
  <img src="assets/iron_all_phases_dataset_space.png" alt="Dataset Atlas for the current BCC, FCC, and HCP iron trajectories" width="92%" />
</p>

Coverage of the available structures, temperatures and magnetic configurations.

<p align="center">
  <img src="assets/iron_phase_phonons_poster.png" alt="Phase-resolved finite-temperature phonon summary for BCC, FCC, and HCP iron" width="68%" />
</p>

Representative finite-temperature phonons; panels use different simulation conditions.

### BCC

<p align="center">
  <img src="assets/bcc_free_energy_vs_volume_4000K_4500K_5000K_5500K_6000K.png" alt="BCC free energy comparison for 4000 K, 4500 K, 5000 K, 5500 K, and 6000 K" width="48%" />
  <img src="assets/bcc_volume_vs_pressure_4000K_4500K_5000K_5500K_6000K.png" alt="BCC pressure-volume comparison for 4000 K, 4500 K, 5000 K, 5500 K, and 6000 K" width="48%" />
</p>

Nonmagnetic Helmholtz free energy and pressure–volume curves, 4000–6000 K.

<p align="center">
  <img src="assets/bcc_phonon_dispersion_overlay_4000K.png" alt="BCC phonon dispersion overlay 4000 K" width="19%" />
  <img src="assets/bcc_phonon_dispersion_overlay_4500K.png" alt="BCC phonon dispersion overlay 4500 K" width="19%" />
  <img src="assets/bcc_phonon_dispersion_overlay_5000K.png" alt="BCC phonon dispersion overlay 5000 K" width="19%" />
  <img src="assets/bcc_phonon_dispersion_overlay_5500K.png" alt="BCC phonon dispersion overlay 5500 K" width="19%" />
  <img src="assets/bcc_phonon_dispersion_overlay_6000K.png" alt="BCC phonon dispersion overlay 6000 K" width="19%" />
</p>

TDEP phonon and DOS overlays at 4000, 4500, 5000, 5500 and 6000 K, left to right.

<p align="center">
  <img src="assets/bcc_phonons_magnetic_vs_nonmagnetic_4000K_a2.55.png" alt="BCC collinear-magnetic and non-magnetic phonon dispersion comparison at 4000 K" width="78%" />
</p>

Collinear magnetic versus nonmagnetic BCC at 4000 K, a ≈ 2.55 Å. The reported single-state free-energy difference is +10.3 meV/Fe; this is not a converged magnetic phase boundary.

<p align="center">
  <img src="assets/bcc_noncollinear_md_dashboard.gif" alt="Synchronized 64-frame noncollinear BCC AIMD and HELD dashboard" width="92%" />
</p>

Noncollinear BCC: 64 synchronized MD frames with local moments and per-step HELD phonons. Short-trajectory diagnostic.

### Round25: constrained noncollinear BCC

![Round25 phonons and MD](docs/round25/round25_dashboard.png)

Archived September 21 analysis: 32 Fe, 4000 K target, 92 fitted configurations. Ensemble TDEP and mean HELD have no imaginary modes on the sampled path; some individual frames do. [Animation and methods](docs/round25/README.md).

![Round25 spin entropy estimate](docs/round25/round25_spin_entropy.png)

Independent-local-moment estimate: **0.95643 kB/Fe**. Constrained-spin model estimate, not measured equilibrium entropy; not added to Helmholtz free energies.

### FCC

<p align="center">
  <img src="assets/fcc_free_energy_vs_volume_4000K_4500K_5000K_5500K_6000K_6500K.png" alt="FCC free energy comparison for 4000 K, 4500 K, 5000 K, 5500 K, 6000 K, and 6500 K" width="48%" />
  <img src="assets/fcc_volume_vs_pressure_4000K_4500K_5000K_5500K_6000K_6500K.png" alt="FCC pressure-volume comparison for 4000 K, 4500 K, 5000 K, 5500 K, 6000 K, and 6500 K" width="48%" />
</p>

Nonmagnetic Helmholtz free energy and pressure–volume curves, 4000–6500 K.

<p align="center">
  <img src="assets/fcc_phonon_dispersion_overlay_4000K.png" alt="FCC phonon dispersion overlay 4000 K" width="48%" />
  <img src="assets/fcc_phonon_dispersion_overlay_4500K.png" alt="FCC phonon dispersion overlay 4500 K" width="48%" />
</p>

<p align="center">
  <img src="assets/fcc_phonon_dispersion_overlay.png" alt="FCC phonon dispersion overlay 5000 K" width="48%" />
  <img src="assets/fcc_phonon_dispersion_overlay_5500K.png" alt="FCC phonon dispersion overlay 5500 K" width="48%" />
</p>

<p align="center">
  <img src="assets/fcc_phonon_dispersion_overlay_6000K.png" alt="FCC phonon dispersion overlay 6000 K" width="48%" />
  <img src="assets/fcc_phonon_dispersion_overlay_6500K.png" alt="FCC phonon dispersion overlay 6500 K" width="48%" />
</p>

TDEP phonon and DOS overlays: 4000/4500 K, 5000/5500 K and 6000/6500 K, by row.

### HCP

<p align="center">
  <img src="assets/hcp_free_energy_vs_volume.png" alt="HCP free energy" width="48%" />
  <img src="assets/hcp_volume_vs_pressure_5000K_eos_std.png" alt="HCP pressure-volume EOS" width="48%" />
</p>

Nonmagnetic Helmholtz free energy and pressure–volume curves at 5000 K.

<p align="center">
  <img src="assets/hcp_phonon_dispersion_overlay.png" alt="HCP phonon dispersion overlay" width="78%" />
</p>

TDEP phonon and DOS overlays across the sampled HCP geometries.

## Documentation

[Dataset](dataset/README.md) · [TDEP workflow](codes/tdep_workflow/README.md) · [HPC setup and submission](hpc/README.md) · [Scripts and usage](docs/usage.md)

## References

- [Hellman et al., TDEP free energies, PRB 87, 104111 (2013)](https://doi.org/10.1103/PhysRevB.87.104111).
- [Hellman & Abrikosov, third-order force constants, PRB 88, 144301 (2013)](https://doi.org/10.1103/PhysRevB.88.144301).
- [Knoop et al., TDEP software, JOSS 9, 6150 (2024)](https://doi.org/10.21105/joss.06150).

[License](LICENSE)
