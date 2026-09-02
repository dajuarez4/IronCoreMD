# Constrained-DLM validation and scaling plan on Lonestar6

## Current decision

The repaired Round 20 pilot is promising but is not yet an accepted production
protocol. At roughly 35 of 100 MD steps it has:

- printed target-vector norms of `2.000000 mu_B` for all 16 species;
- retained local moments of roughly `1.5--1.8 mu_B/Fe`;
- net magnetization below about one percent of the summed local moments;
- decreasing force and pressure during the initial equilibration;
- large instantaneous temperature fluctuations, with recent values returning
  toward the 4000 K target;
- some large SCF-force-correction warnings that still require a full-trajectory
  count and magnitude assessment.

Do not begin a 128-atom production trajectory from this partial result.

## Gate 1: finish and accept the 16-atom pilot

Require all of the following:

1. `100/100` dynamics steps, `JOB DONE`, and Slurm exit code `0:0`.
2. No QE fatal errors, diagonalization failures, or electronic nonconvergence.
3. Every printed constrained-vector norm remains `2.0 mu_B`.
4. The last-50-step mean temperature is reasonably close to 4000 K and has no
   sustained upward or downward drift.
5. Local moments remain substantial and the net cell magnetization remains a
   small fraction of the summed local moments.
6. SCF-force-correction ratios are acceptable. If warnings remain frequent or
   corrections are commonly above roughly 20 percent of total force, repeat a
   short pilot with a tighter MD `conv_thr` before scaling.

## Gate 2: reinforce sampling at 16 atoms

Run the remaining independent Round 20 spin/velocity seeds before increasing
cell size. It is acceptable to split replicas between Jakar and Lonestar6 only
when each cluster has independently passed the same runtime constraint test.
Record for every replica:

- QE executable SHA256 and patch provenance;
- compiler and MPI modules;
- pseudopotential SHA256;
- spin and velocity seeds;
- input-file SHA256;
- completed steps, temperature statistics, moment statistics, and force-quality
  statistics.

Compiler-level numerical differences between clusters are acceptable for an
ensemble, but outputs made with different physical constraint formulas are not.

## Gate 3: cell-size ladder

Use short pilots before long production:

| Cell | Purpose | Suggested first run |
|---|---|---|
| 16 atoms | validate thermostat, constraints, and SCF settings | 8 independent x 100 steps |
| 32 or 48 atoms | test finite-size behavior and cost scaling | 1--2 independent x 50--100 steps |
| 64 atoms | verify memory, wall time, and force quality | 1 x 20--50 steps |
| 128 atoms | benchmark first, then production | SCF ramp plus 5--10 MD steps |

Do not extrapolate a 16-core `vm-small` configuration directly to 128 atoms.
Use the full-node `normal` queue and benchmark at least one and two nodes. Start
with pure MPI and `OMP_NUM_THREADS=1`; compare wall time per SCF, memory, and
parallel efficiency before selecting production resources. Lonestar6 full CPU
nodes provide 128 cores, and MPI applications must be launched with `ibrun`.

For the 128-atom one-species-per-atom construction, verify all three facts in
the actual 128-species SCF output before MD:

1. QE parsed `number of atomic types = 128` (the separate maximum-species field
   may print `**` because QE 7.5 uses a two-character integer format there).
2. All 128 constrained-vector norms equal `2.0 mu_B` within output precision.
3. The SCF ramp completes without excessive correction forces or moment
   collapse.

## Queue selection

- Use `vm-small` only for builds and calculations using at most 16 cores and
  less than 29 GB of memory.
- Use `normal` for full-node 32/48/64/128-atom benchmarks and production.
- Run `qlimits` immediately before submission because TACC queue limits may
  change.
- Do not request specific compute-node hostnames.

The current official Lonestar6 guide is
<https://docs.tacc.utexas.edu/hpc/lonestar6/>.
