# Round25: 32 atoms from Round22 final structure, 400 new MD steps

Prepared for Jakar: jakar_general, medium, 1 node, 32 MPI ranks, 48 hours.
Use the same patched QE 7.5 executable as Round23:
/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x (atomic-constraint fix and ntypx=128).

Initial positions are the final propagated coordinates printed at Round22
step 162. No final velocities or checkpoint were available locally. New
Gaussian Maxwell velocities (seed 250032) have zero COM motion and are rescaled
to 4000 K. These are not reconstructed Round22 velocities. The new run starts
from scratch and repeats staged SCF preparation; no Round22 checkpoint is needed.
The existing 162 steps are not included in the requested 400 new MD steps.

Unchanged from Round23: 32 Fe species and target spin directions, 4x2x2 cell,
a=2.39 Å, same UPF, cutoffs 71/496 Ry, FD smearing 0.025334 Ry, no SOC,
atomic constraints, final lambda=0.1 Ry, MD conv_thr=1e-4, mixing_beta=0.01,
local-TF, SVR at 4000 K, nraise=20, dt=20.67 and all other SYSTEM/ELECTRONS/IONS
entries. Initial positions/velocities, prefix and MD length are changed.

From this directory on Jakar:

```bash
sbatch run_round25_jakar.sbatch
python3 check_round25.py
```

The workflow runs staged SCF then MD, requesting a clean stop within a 46-hour
cumulative QE budget. After a clean time-limit stop, submit the same script
again; it continues the current stage and requests only remaining MD steps.
Numbered .runNNN.in/.out files preserve each segment. Hard interruptions or
SCF failures require inspection and are not automatically restarted. Keep
qe_tmp and do not submit concurrent runs of this case. A lock prevents overlap.
No remote QE or real checkpoint restart has been tested for this package.

For comparisons, treat this as a new branch from Round22 final geometry.
Analyze equilibration after the velocity initialization and thermostat change;
do not label the combined Round22/Round25 data as a single nraise=20 trajectory.
Round23 starts from a different geometry, so this is not an identical-initial-state
comparison. Compare consistent equilibrated windows and analysis settings.
