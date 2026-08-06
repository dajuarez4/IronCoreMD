# Ferromagnetic BCC Fe: spin-resolved bands and DOS

This standalone Quantum ESPRESSO workflow is tracked under `workflows/` in the
`IronCoreMD` repository. It evaluates collinear ferromagnetic BCC Fe with the
PBE PAW pseudopotential at two fixed cubic lattice parameters:

- `a_2.50`: `a = 2.50 angstrom`
- `a_2.78`: `a = 2.78 angstrom`

These are lattice parameters, not pressures. Because `tstress = .true.` in
the SCF inputs, QE will report the hydrostatic stress/pressure associated with
each fixed lattice parameter in `01_scf.out`.

Each calculation uses an ideal `4 x 4 x 4` repetition of the conventional
two-atom BCC cell: `128 Fe atoms`, `ibrav = 0`, `nspin = 2`, and a positive Fe
starting magnetization. The cubic supercell edges are `10.00 angstrom` for
`a = 2.50 angstrom` and `11.12 angstrom` for `a = 2.78 angstrom`. No spin-orbit
coupling or noncollinear magnetism is enabled.

## Calculation sequence

Each case contains:

1. `01_scf.in`: spin-polarized ground-state SCF on a `6 x 6 x 6` supercell mesh.
2. `02_bands.in`: `pw.x` band calculation along
   the folded simple-cubic supercell path `Gamma-X-M-Gamma-R-X | M-R`.
3. `03_bands_up.in`: `bands.x` extraction for `spin_component = 1`.
4. `04_bands_down.in`: `bands.x` extraction for `spin_component = 2`.
5. `05_nscf.in`: dense `8 x 8 x 8` tetrahedron NSCF calculation.
6. `06_projwfc.in`: `projwfc.x` total and orbital-projected DOS.
7. `run_jakar.sbatch`: complete Jakar calculation in the required order.

For the collinear spin-polarized calculation, `projwfc.x` writes
`E, DOSup, DOSdw, PDOSup, PDOSdw` to `*.pdos_tot`. The separate `bands.x`
inputs write spin-up and spin-down `*.gnu` band files.

## Pseudopotential

The workflow reuses:

`Fe.pbe-spn-kjpaw_psl.1.0.0.UPF`

This is the PBE scalar-relativistic PAW Fe pseudopotential already used by the
Fe simulations. Its suggested cutoffs, `ecutwfc = 71 Ry` and
`ecutrho = 496 Ry`, are used unchanged.

The pseudopotential binary is not duplicated in this workflow directory.
Before transferring to Jakar, place a verified copy at
`pseudo/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF`:

```bash
rsync -av workflows/bcc_ferromagnetic_bands_pdos_pbe/ \
  jakar:/your/scratch/path/bcc_ferromagnetic_bands_pdos_pbe/
```

The project HCP magnetic dataset contains the same named potential for its
reproducible QE inputs; verify its checksum and licensing requirements before
reusing it for production calculations.

## Submit on Jakar

From the top-level folder:

```bash
./submit_all.sh
```

Or submit one lattice parameter:

```bash
cd a_2.50
sbatch run_jakar.sbatch
```

The scripts use the established Jakar settings:

- account: `jakar_partner_physj`
- QoS: `jakar_medium_physj`
- partition: `medium`
- one exclusive node, 38 MPI tasks
- OpenMPI: `openmpi4/4.1.5`
- QE: `/shared/quantum-espresso/bin`

The scripts follow the supplied, working 128-atom Jakar template exactly:
38 MPI ranks on one exclusive node with the shared OpenMPI Quantum ESPRESSO
build. These larger supercells have enough FFT planes for all 38 ranks. All
electronic steps use the robust conjugate-gradient diagonalizer.

## Important outputs

For each case:

- `01_scf.out`: Fermi energy and converged total magnetization.
- `bands_up.dat.gnu`: spin-up electronic bands in eV.
- `bands_down.dat.gnu`: spin-down electronic bands in eV.
- `fe_a*.pdos_tot`: total DOS and summed PDOS, separated into up/down columns.
- `fe_a*.pdos_atm#1(Fe)_wfc#*(*)`: Fe orbital-resolved spin DOS.
- `06_projwfc.out`: Lowdin charges, projections, and spilling parameter.

After downloading the completed outputs, make combined band/DOS figures:

```bash
python3 plot_results.py
```

The plotting script reads the Fermi energy from `01_scf.out`, shifts band and
DOS energies to `E-E_F`, and plots the down-spin DOS with a negative sign for
visual separation.

## Convergence note

The `6^3` SCF mesh and `8^3` DOS mesh preserve the reciprocal-space density of
the earlier primitive-cell meshes after the fourfold increase in each real-space
direction. The SCF uses `nbnd = 1280`; bands and NSCF use `nbnd = 1664`.
Because a 128-atom perfect supercell folds the primitive-cell bands into a much
smaller Brillouin zone, the band plot will contain many folded bands.

Quantum ESPRESSO input references:

- https://www.quantum-espresso.org/Doc/INPUT_PW.html
- https://www.quantum-espresso.org/Doc/INPUT_BANDS.html
- https://www.quantum-espresso.org/Doc/INPUT_PROJWFC.html
