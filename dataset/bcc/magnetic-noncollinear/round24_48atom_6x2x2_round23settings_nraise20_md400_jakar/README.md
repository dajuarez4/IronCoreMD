# Round24: 48 atoms, 400 MD steps, nraise=20

New Jakar calculation using Round23 electronic/ionic parameters. Round23 is
retained unchanged. This is not a restart of its 32-atom checkpoint.

- BCC 6x2x2 conventional supercell; 48 Fe; a=2.39 Å.
- Cell 14.34 x 4.78 x 4.78 Å, volume per atom 6.8259595 Å³.
- Round23 first 32 Cartesian positions retained; its first 16-atom tile added
  along x. Velocities extended, COM removed, rescaled to 4000 K.
- 48 species / 48 distinct directions; first 32 target directions and labels
  retained; 16 added directions in antipodal pairs. Target moment 2 μB/Fe,
  summed target vector zero to floating-point precision.
- nstep=400, nraise=20, dt=20.67, SVR 4000 K, Verlet; no SOC.
- Same UPF, ecut 71/496 Ry, FD smearing 0.025334 Ry, atomic constraints.
- Same SCF preparation: lambda 0.005 -> 0.020 -> 0.100 Ry;
  conv_thr 4e-3 -> 1e-3 -> 1e-4 Ry. MD conv_thr=1e-4 Ry,
  mixing_beta=0.01, local-TF, mixing_ndim=20, electron_maxstep=1000.
- Requires the same patched QE 7.5 with atomic-constraint fix and ntypx=128.

On Jakar, from this directory:

```bash
sbatch run_round24_jakar.sbatch
python3 check_round24.py
```

The script retains Round23's Jakar resources: 1 node, 32 MPI ranks, 48 hours.
48 atoms does not mean 48 MPI ranks. Wall time has not been benchmarked for
this new case. The workflow reserves 2 hours for overhead and asks QE to stop
cleanly using max_seconds within a cumulative 46-hour budget. Runtime inputs
are separately numbered; original input templates remain unchanged.

After a clean time-limit stop, submit the SAME sbatch script again. It skips
completed SCF stages and continues the current stage from checkpoint. For MD,
it requests only the remaining steps up to 400 total across saved output
segments. A lock prevents concurrent launches in the same case directory.
It refuses automatic continuation after a hard interruption or SCF error.
Inspect those failures and checkpoint integrity before any manual recovery.
Never delete qe_tmp or move it separately while the calculation is incomplete.

Validation is local (input parsing, structure, target moments, velocities,
parameter comparisons and workflow control-flow tests). No remote QE run or
checkpoint restart has been performed yet. Reported SCF completion is required
at each stage; it does not independently establish phonon quality.

Increasing the x length preserves the short 4.78 Å transverse dimensions.
For the existing minimum-image HELD setup, longer interaction cutoffs are not
justified solely by moving from 32 to 48 atoms. Compare with Round23 using the
same analysis and account for changed cell size as well as trajectory length.
