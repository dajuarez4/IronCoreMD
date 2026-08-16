# Round 10: patched-QE lambda=0.10, 50-step MD validation

## Why this round

Both Round 9 pilots completed 10 MD steps without a failed SCF, QE routine
error, or Cholesky error. `lambda = 0.10 Ry` was selected because it retained a
mean `1.573 Bohr magneton/Fe`, required the same 60 electronic iterations as
`lambda = 0.05`, and controlled the small-cell temperature more effectively.

Round 10 extends only the selected case to 50 steps. It remains a validation
trajectory, not production sampling.

## Fixed design

- BCC Fe 8-atom `2x2x1` cell, `a = 2.39 angstrom`;
- thermally displaced starting positions;
- the same 4000 K seed-`400054` velocities used in Round 9;
- noncollinear atomic target `2 Bohr magneton/Fe`;
- corrected private QE 7.5 atomic constraint, `lambda = 0.10 Ry`;
- fresh SCF ramp `0.005 -> 0.020 -> 0.100 Ry`;
- SCF thresholds `4e-3 -> 1e-3 -> 3e-4 Ry`;
- MD `conv_thr = 4e-3 Ry`, accepted for exploratory MD;
- Gamma point, Fermi-Dirac smearing, `degauss = 0.025334 Ry`;
- 50 steps at `dt = 20.67 a.u.` (approximately 0.5 fs);
- SVR thermostat, target 4000 K, `nraise = 10`.

The case performs its own three-stage preconditioning. It does not copy or
depend on a Round 9 `qe_tmp` directory.

## Copy to Lonestar6

From the Mac:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear
scp -r round10_qe75_patched_lambda0p100_md50_lonestar6 \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/
```

## Verify the patched executable

On Lonestar6:

```bash
QE_FIX="$SCRATCH/apps/qe-7.5-atomic-fixed"
test -x "$QE_FIX/bin/pw.x" && echo "pw.x is executable"
sha256sum -c "$QE_FIX/pw.x.sha256"
cat "$QE_FIX/ATOMIC_CONSTRAINT_PATCH_INFO.txt"
```

The checksum must report `OK`.

## Submit

```bash
cd /work/11381/dajuarez4/ls6/round10_qe75_patched_lambda0p100_md50_lonestar6
mkdir -p logs
find cases \( -name '*.out' -o -name qe_tmp \) -print
```

The `find` command must print nothing. Then:

```bash
sbatch run_round10_lonestar6.sbatch
```

The job uses one `vm-small` node and 8 MPI ranks. It refuses existing output or
restart state to prevent attempts from being mixed.

## Monitor

```bash
squeue -u "$USER"
bash check_round10.sh
bash plot_round10.sh --live 30
```

Replace `JOB_ID` below with the submitted job number:

```bash
tail -F logs/ls6_round10_md50.JOB_ID.out
tail -F logs/ls6_round10_md50.JOB_ID.err
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,ExitCode
```

## Acceptance criteria

Proceed toward longer/larger MD only if this run:

1. reaches 50 `Entering Dynamics` records and `JOB DONE`;
2. has no failed SCF, QE routine, Cholesky, or cdiaghg error;
3. keeps each converged MD cycle below `conv_thr = 4e-3 Ry`;
4. avoids cycles approaching `electron_maxstep = 500`;
5. maintains mean local moments preferably above `1 Bohr magneton/Fe`;
6. shows bounded temperature, force, pressure, and moment behavior.

The dashboard displays SCF accuracy, iterations per cycle, temperature, and
mean atom-resolved local moment. A single 25 fs trajectory is still too short
for production averages, but it is sufficient to expose immediate electronic
or thermostat instability.

## Lightweight result archive

After successful completion:

```bash
module reset
module load python3/3.9.7
python3 collect_round10_npz.py --require-complete
ls -lh round10_results_light.npz round10_md50_dashboard.png
```

The collector fixes the Round 9 pressure-pattern issue: it cannot mistake
`nstep=50` for `P = ...`. The safe `allow_pickle=False` archive includes SCF
histories, trajectory positions, atom-resolved forces and spin vectors,
constraint vectors, temperatures, energies, forces, and pressures. It excludes
full outputs and `qe_tmp`.

Copy the lightweight products back:

```bash
scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round10_qe75_patched_lambda0p100_md50_lonestar6/round10_results_light.npz \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round10_qe75_patched_lambda0p100_md50_lonestar6/

scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round10_qe75_patched_lambda0p100_md50_lonestar6/round10_md50_dashboard.png \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round10_qe75_patched_lambda0p100_md50_lonestar6/
```

Do not copy or commit `qe_tmp`, full `.out` files, or scheduler logs. For a
rerun, preserve the old folder and start from a fresh repository copy.

## Completed result

Round 10 reached all 50 ionic steps and `JOB DONE`. All electronic cycles
converged below `4e-3 Ry`: the mean was 6.2 electronic iterations per cycle,
the maximum was 13, and the largest final-cycle accuracy was
`3.98683e-3 Ry`. No fatal QE or Cholesky error occurred.

The atomic constraint retained the noncollinear local moments: the trajectory
mean was `1.600 Bohr magneton/Fe` and the final mean was
`1.623 Bohr magneton/Fe`. Temperature was not yet equilibrated at 4000 K: its
overall mean was about 3160 K, its final value was about 3265 K, and the final
10-record mean was about 2973 K. Consequently, the run validates the
electronic and magnetic setup but is too short for thermodynamic averages.
Round 11 extends this setup to 400 steps with the requested `nraise = 20`.
