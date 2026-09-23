# round1_banstructure — Round24 and Round25 electronic bands on Jakar

This package prepares an inexpensive first electronic-band calculation for ONE
saved MD snapshot per round and includes a comparative MD dashboard. It does
not contain computed electronic bands yet, and the source wavefunctions are
not available on this Mac. Copy/import them locally on Jakar using the paths below.

## What is reused and what is calculated

1. Import a separate copy of the inactive source case's entire qe_tmp (density,
   XML, wavefunctions and supporting checkpoint files). Never use symlinks or
   run the band calculation in the active MD outdir.
2. Read positions/cell from output/atomic_structure in the checkpoint XML.
   This geometry may differ from the last printed, already-propagated MD
   position block. Do not assign it a completed MD step based on filenames.
3. Run a fixed-ion Gamma SCF using startingpot='file', startingwfc='file', and
   restart_mode='from_scratch'. This genuinely uses the saved MD wavefunctions
   as an initial guess and reconverges the density at the XML geometry. The
   checkpoint need not already have converged at that geometry.
4. Run calculation='bands' on the converged snapshot potential. The MD contains
   only Gamma wavefunctions; it does not contain wavefunctions at other k.
   New k states must be solved, so this stage deliberately uses atomic+random
   guesses rather than incorrectly matching Gamma files to new k points.
5. Run bands.x and plot the results relative to the preceding Gamma SCF Fermi
   level. There is no separate collinear up/down channel in this spinor calculation.

Same patched QE 7.5, pseudopotential, 71/496 Ry cutoffs, FD smearing 0.025334 Ry,
local atomic constraints lambda=0.1 Ry, no SOC. Gamma SCF retains the source
ELECTRONS parameters, including conv_thr=1e-4. The path eigensolver uses
CG, diago_thr_init=1e-6 and diago_full_acc=true for unoccupied bands.
Band counts are read from the checkpoint (current outputs: 922 for 48 atoms,
614 for 32 atoms); no extra unoccupied states are added by default.

Cost is limited by one snapshot and only 11 path points. This still solves
hundreds of spinor bands at each k and is not a negligible job. Neither the
Gamma potential nor this sparse path is a k-point-converged band study.
No quantitative runtime guarantee is possible without a pilot.

## Reciprocal-space interpretation

The path is Gamma–X–S–Y–Gamma–Z in fractional SUPERCELL reciprocal coordinates:
(0,0,0) -> (1/2,0,0) -> (1/2,1/2,0) -> (0,1/2,0) -> (0,0,0) -> (0,0,1/2).
There are two intervals per segment, including all endpoints = 11 k points.
Names denote geometric cell-boundary points, not exact symmetry labels of a
thermally disordered magnetic snapshot. Both structures have different supercell
Brillouin zones. Results are folded supercell spectra; they are not primitive
BCC Gamma–H–N–P dispersions and are not unfolded or thermally averaged.
The bands.x output is energy ordered (lsym=false, no_overlap=true), without
symmetry labels or claims of continuous orbital character across crossings.

## Import on Jakar after the source workflow has stopped cleanly

The supplied paths are:

- Round24: /scratch/dajuarez4/testing/non_collinear/round24_48atom_6x2x2_round23settings_nraise20_md400_jakar/cases/conv1e-4_beta1e-2_48atom
- Round25: /scratch/dajuarez4/testing/non_collinear/round25_32atom_round22final_nraise20_md400_jakar/cases/conv1e-4_beta1e-2_32atom

From the extracted main directory (Python 3 with NumPy required):

```bash
bash import_checkpoints.sh
```

The importer checks the workflow lock, presence of density and wavefunctions,
atom count/labels, cell, electron count, spinor mode, cutoffs and UPF checksum.
It checks source file size/mtime stability across copying. A lock is required;
this is designed for the supplied round24/25 workflow directories. Do not copy
from an actively running simulation or from a hard-killed, unverified checkpoint.
The importer refuses existing destination checkpoints. If only one source is
ready, run that individual python3 prepare_checkpoint.py command from the shell
script instead; successful imports do not need to be repeated.

The *.in.example files are reviewable examples from the LAST COMPLETED FORCE
geometry in the local logs. They are not executable imported checkpoints.
Importing produces the real 01_snapshot_scf.in, 02_bands.in and 03_bands_post.in
from the XML geometry, with checkpoint.json and kpath.json provenance.
Increase path resolution at import with --points-per-segment 4 (21 points).

## Submit in two allocations per round

```bash
cd round24
sbatch run_jakar.sbatch
```

The first submission runs only the snapshot SCF. Once it reports completion,
submit the SAME file again to run the band path and bands.x. Repeat for round25.
Resources: Jakar jakar_general, jakar_medium_general, medium, 1 node / 32 MPI
ranks, 48 h. The executable directory is /gpfs/scratch/dajuarez4/qe-7.5/bin.
Use the same patched pw.x as MD, not the unpatched system QE.

Each PW stage has max_seconds=160000 to allow a clean stop. An incomplete stage
is not marked complete, and existing outputs are not silently overwritten.
If a stage hits the limit, inspect and deliberately prepare its restart before
resubmitting; this initial band package does not automate interrupted-k-path
recovery. No MD files or checkpoints in the source directories are modified.

After band completion (Python 3, NumPy and Matplotlib):

```bash
cd ..
python3 plot_bands.py
```

The script produces dashboard/electronic_bands.png and round*/electronic_bands.npz.
Only computed results populate the electronic panel. Copy the dashboard and
small band files back to the Mac; you need not transfer wavefunction files.

## MD dashboard already generated

Open dashboard/index.html in a browser, or view round24_round25_dashboard.gif
and .png in that directory. Current local files contain 6 completed round24
steps and 37 completed round25 steps, not full trajectories. The animation holds
the last round24 structure after its available frames are exhausted and labels it.

Structure/spins correspond to the force-evaluated configuration; temperatures
are the post-propagation values printed for each MD step. Pressure is the
printed electronic stress P, without adding ionic kinetic stress. Full available
histories are displayed behind the moving cursor. Rounds have different starting
structures and lengths; their pressure difference cannot be attributed only to
atom count. The dashboard does not prove live scheduler status.

To refresh from a selected new single output segment:

```bash
python3 make_dashboard.py --round24-output /path/to/round24/md_400steps.run001.out --round25-output /path/to/round25/md_400steps.run001.out
```

This script reads one segment per round; it does not automatically merge restart
segments. CSVs and NPZs with source hashes accompany the dashboard.
Dependencies for dashboard generation: Python 3, NumPy, Matplotlib, Pillow.

## References and validation scope

- QE PW input: https://www.quantum-espresso.org/Doc/INPUT_PW.html
- Fixed-potential bands and geometry from the SCF save:
  https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
- bands.x options: https://www.quantum-espresso.org/Doc/INPUT_BANDS.html

Local validation covers example parsing, k-path endpoints, parameter retention,
mock checkpoint import/active-lock refusal and data/plot generation. No physical
band computation or real round24/25 checkpoint import has been executed locally.
