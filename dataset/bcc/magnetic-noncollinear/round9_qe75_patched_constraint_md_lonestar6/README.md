# Round 9: patched-QE noncollinear MD comparison on Lonestar6

## Scientific objective

Round 8 established that the corrected QE 7.5 atomic constraint prints the
intended `2 Bohr magneton/Fe` target and that both `lambda = 0.05` and
`0.10 Ry` retain more than `1 Bohr magneton/Fe` in static SCF calculations.
Round 9 compares those two strengths during ionic motion.

This is a controlled pair of 10-step exploratory MD pilots:

- 8-atom BCC Fe `2x2x1` cell at `a = 2.39 angstrom`;
- identical thermally displaced positions;
- identical 4000 K velocities generated with seed `400054`;
- identical noncollinear directions and `2 Bohr magneton/Fe` targets;
- final constraint strengths `0.05` and `0.10 Ry`;
- Gamma point, Fermi-Dirac smearing, `degauss = 0.025334 Ry`;
- `dt = 20.67 a.u.` (approximately 0.5 fs), SVR thermostat at 4000 K;
- MD `conv_thr = 4e-3 Ry`, explicitly accepted only for exploratory MD.

Each case independently repeats the successful three-stage Round 8
preconditioning path (`4e-3 -> 1e-3 -> 3e-4 Ry`) before MD. This costs a few
extra minutes but makes the folder reproducible and avoids copying or mixing
fragile `qe_tmp` restart directories between rounds.

## Required executable

This round must use the private patched QE 7.5 executable built in Round 8:

```bash
QE_FIX="$SCRATCH/apps/qe-7.5-atomic-fixed"
test -x "$QE_FIX/bin/pw.x" && echo "pw.x is executable"
sha256sum -c "$QE_FIX/pw.x.sha256"
cat "$QE_FIX/ATOMIC_CONSTRAINT_PATCH_INFO.txt"
```

The checksum must report `OK`. Do not substitute the shared Lonestar6 QE 7.3
or an unpatched QE 7.5 build; those versions do not apply the intended valence
scaling to this constraint target.

## Copy to Lonestar6

From the Mac:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear
scp -r round9_qe75_patched_constraint_md_lonestar6 \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/
```

## Submit

On Lonestar6:

```bash
cd /work/11381/dajuarez4/ls6/round9_qe75_patched_constraint_md_lonestar6
mkdir -p logs
find cases \( -name '*.out' -o -name qe_tmp \) -print
```

The `find` command must print nothing for a fresh run. Then submit:

```bash
sbatch run_round9_serial_lonestar6.sbatch
```

The two cases run sequentially on one `vm-small` node with 8 MPI ranks. The
runner refuses existing output or restart state rather than mixing attempts.

## Monitor and plot

```bash
squeue -u "$USER"
bash check_round9.sh
bash plot_round9.sh --live 30
```

For scheduler output, replace `JOB_ID` with the submitted job number:

```bash
tail -F logs/ls6_round9_serial.JOB_ID.out
tail -F logs/ls6_round9_serial.JOB_ID.err
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,ExitCode
```

The Bash/gnuplot dashboard compares accuracy histories, electronic iterations
per SCF cycle, ionic temperature, and the mean of the eight atom-resolved local
moment magnitudes.

## Acceptance criteria

A useful pilot must:

1. finish all 10 ionic steps with `JOB DONE`;
2. have no `Error in routine`, Cholesky/cdiaghg error, or failed SCF cycle;
3. keep the final MD SCF accuracy below `4e-3 Ry`;
4. retain a mean local moment preferably above `1 Bohr magneton/Fe`;
5. avoid isolated SCF cycles approaching `electron_maxstep = 500`;
6. maintain reasonable temperature and force behavior for this short pilot.

Select `lambda = 0.10 Ry` for longer MD only if its moment retention improves
without a severe SCF-cost or force penalty relative to `0.05 Ry`. These 10
steps diagnose stability; they are not an equilibrated production trajectory.

## Create the lightweight archive

After both cases finish:

```bash
module reset
module load python3/3.9.7
python3 collect_round9_npz.py --require-complete
ls -lh round9_results_light.npz round9_md_dashboard.png
```

The NPZ uses no pickle/object arrays and excludes full outputs and `qe_tmp`:

```python
import numpy as np
data = np.load("round9_results_light.npz", allow_pickle=False)
print(data["case_name"])
```

Copy the lightweight results back to the repository folder:

```bash
scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round9_qe75_patched_constraint_md_lonestar6/round9_results_light.npz \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round9_qe75_patched_constraint_md_lonestar6/

scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round9_qe75_patched_constraint_md_lonestar6/round9_md_dashboard.png \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round9_qe75_patched_constraint_md_lonestar6/
```

Do not copy or commit `qe_tmp`, full `.out` files, or scheduler logs. To rerun,
preserve the old folder and start from a fresh repository copy.
