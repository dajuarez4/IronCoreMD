# Round 20: constrained-DLM ensemble diagnostic on Jakar

This folder tests a more defensible paramagnetic sampling protocol before
scaling to larger cells or calculating paramagnetic phonons. It contains eight
independent 16-atom constrained-DLM replicas. Replica 0 is the pilot; replicas
1--7 are submitted only after the pilot is accepted.

## What changed relative to Round 17

- Eight independent random magnetic textures.
- Each texture contains eight random directions and their exact antipodes.
- One species and one target direction per atom (`Fe01`--`Fe16`).
- Every replica has zero target-spin vector sum by construction.
- Independent Maxwell-Boltzmann velocities, with zero COM velocity and exactly
  4000 K, are used for every replica.
- MD is a short 100-step diagnostic.
- MD `conv_thr` is tightened from `4e-3` to `1e-3 Ry`.
- SVR `nraise` is reduced from 20 to 5 for stronger thermal coupling.

The SCF ramp remains `lambda=0.005 -> 0.020 -> 0.100` with convergence
thresholds `4e-3 -> 1e-3 -> 3e-4 Ry`. The MD constraint remains
`lambda=0.100`. This is a constrained DLM ensemble; it does not claim to be
explicit spin dynamics.

## Jakar configuration

```text
account   = jakar_general
qos       = jakar_medium_general
partition = medium
nodes     = 1 per replica
MPI tasks = 16 per replica
QE        = /gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x
```

## Submission order

Copy the complete folder to Jakar and submit only the pilot:

```bash
sbatch run_round20_pilot_jakar.sbatch
```

After it finishes:

```bash
bash check_round20_jakar.sh
```

Accept the pilot when all three SCF stages and 100 MD steps complete, the last
50-step mean temperature is reasonably close to 4000 K, the local moments are
retained, and the number of `SCF correction compared to forces is large`
warnings is zero or substantially reduced relative to Round 17.

Only after accepting the pilot, submit replicas 1--7:

```bash
sbatch run_round20_ensemble_jakar.sbatch
```

The ensemble job is an array `1-7%2`, so at most two replicas run concurrently.
Monitor with:

```bash
squeue -u dajuarez4
bash check_round20_jakar.sh
tail -f JOBID.log
tail -f logs/ARRAY_JOBID_TASKID.log
```

Do not submit the ensemble array before checking replica 0, and do not mix old
outputs or `qe_tmp` directories into these cases.
