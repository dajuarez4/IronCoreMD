# Round 8: patched QE 7.5 atomic-constraint validation on Lonestar6

## Purpose and scope

This folder is a fully reproducible **static SCF validation**, not yet a
production MD run. We use a thermally displaced 8-atom BCC Fe `2x2x1` cell at
`a = 2.39 angstrom`, noncollinear magnetism, one Fe type per atom, and a Gamma
point calculation. The goal is to identify a constraint strength that both:

1. converges through the staged SCF path; and
2. preserves a physically meaningful local moment near the requested
   `2 Bohr magneton/Fe` target.

Only after this validation succeeds should the selected settings be tested in
short MD and then in a larger cell.

## What we did and learned

Round 7 proved that the shared Lonestar6 QE 7.3 constraint does not represent
the intended `2 Bohr magneton/Fe` target. The QE 7.5 release has the same
`mcons=starting_magnetization` implementation. Current upstream QE multiplies
the target by pseudopotential valence. Round 8 applies only that upstream
four-line correction to a private QE 7.5 source tree and never modifies TACC's
shared installation.

The earlier Round 6/7 scan appeared to converge for several lambda values, but
the printed constrained moment and retained local moments showed that it was
actually constraining about `0.125 Bohr magneton`, not the intended `2`. For
the Fe pseudopotential used here, `starting_magnetization = 0.125` represents
`2/16`; the affected QE implementation omitted the multiplication by the Fe
valence `zv = 16`. Therefore those runs remain useful convergence diagnostics,
but they cannot select a valid 2-Bohr-magneton constraint for MD.

Round 8 changes only this mapping in a private QE 7.5 executable. It scans
final lambda values `0.01`, `0.02`, `0.05`, and `0.10`, using gentler lambda
ramps and SCF thresholds `4e-3 -> 1e-3 -> 3e-4 Ry`. No velocities or MD steps
are included because this round must first verify the electronic constraint.

## Folder contents

- `qe75_atomic_constraint_zv.patch`: minimal four-line correction.
- `prepare_qe75_source.sh`: safely unpacks, patches, and verifies official QE
  7.5 source under `$SCRATCH/src`.
- `build_qe75_patched_lonestar6.sbatch`: builds and installs a private `pw.x`.
- `generate_round8.py`: regenerates all four cases and their staged inputs.
- `run_round8_serial_lonestar6.sbatch`: executes all cases serially on one
  `vm-small` node using 8 MPI ranks.
- `check_round8.sh`: reports convergence, printed target, measured moment, and
  retention.
- `plot_round8.sh`: creates the SCF/constraint diagnostic dashboard with
  Bash and gnuplot.
- `collect_round8_npz.py`: creates a lightweight result archive without QE
  scratch files or full text outputs.
- `cases/`: generated input decks and metadata.
- `pseudo/`: the exact Fe PAW pseudopotential used by the inputs.

## 1. Copy this folder from the repository

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear
scp -r round8_qe75_patched_constraint_lonestar6 \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/
```

## 2. Obtain and prepare QE 7.5 source

On Lonestar6, download the official source into `$SCRATCH` (or upload the same
tarball from the Mac if outbound download is unavailable):

```bash
cd "$SCRATCH"
curl -L -o q-e-qe-7.5.tar.gz \
  https://gitlab.com/QEF/q-e/-/archive/qe-7.5/q-e-qe-7.5.tar.gz
cd /work/11381/dajuarez4/ls6/round8_qe75_patched_constraint_lonestar6
bash prepare_qe75_source.sh "$SCRATCH/q-e-qe-7.5.tar.gz"
```

The preparation script refuses an existing target and verifies the source
change. It creates `$SCRATCH/src/q-e-qe-7.5-atomic-fixed`. The finished
executable is installed under your personal `$SCRATCH/apps` directory, not in
TACC's shared application tree.

## 3. Build the private executable

```bash
mkdir -p logs
sbatch build_qe75_patched_lonestar6.sbatch
```

Monitor the build with `squeue -u "$USER"`. After completion:

```bash
ls -lh "$SCRATCH/apps/qe-7.5-atomic-fixed/bin/pw.x"
sha256sum -c "$SCRATCH/apps/qe-7.5-atomic-fixed/pw.x.sha256"
cat "$SCRATCH/apps/qe-7.5-atomic-fixed/ATOMIC_CONSTRAINT_PATCH_INFO.txt"
```

Warnings such as `ifort ... ignoring unknown option
'-fallow-argument-mismatch'` and the `scan_ibrav.f90` format remarks are not
fatal. An older version of this script executed `pw.x -v` after installation;
QE interpreted it as a calculation with no input and printed `from
test_input_xml: input file not opened or empty`. That message did not establish
that compilation failed. This workflow now verifies the installed executable
with its file, patch marker, and SHA-256 checksum without launching an empty QE
calculation.

## 4. Run the four static diagnostics

```bash
find cases \( -name '*.out' -o -name qe_tmp \) -print
sbatch run_round8_serial_lonestar6.sbatch
```

The `find` command should print nothing for a fresh attempt. The four final
lambda values are `0.01`, `0.02`, `0.05`, and `0.10`; each uses a gentler
lambda ramp while keeping the corrected `2 Bohr magneton/Fe` target active.

```bash
bash check_round8.sh
bash plot_round8.sh --live 30
```

The decisive `target_print` column must be near `2.0`. If it is near `0.125`,
stop: the wrong executable was used. A scientifically useful case must also
converge and retain preferably at least `1 Bohr magneton/Fe`.

Useful Lonestar6 monitoring commands are:

```bash
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,ExitCode
tail -F logs/ls6_round8_serial.JOB_ID.out
tail -F logs/ls6_round8_serial.JOB_ID.err
bash check_round8.sh
```

The batch script intentionally refuses a case that already contains `qe_tmp`
or output files. This prevents state from separate attempts being mixed. To
rerun, preserve the old attempt in a different directory and start from a
fresh copy of this repository folder.

## 5. Collect results

```bash
module reset
module load python3/3.9.7
python3 collect_round8_npz.py --require-terminal
ls -lh round8_results_light.npz
```

Copy back only the lightweight scientific products when the run finishes:

```bash
scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round8_qe75_patched_constraint_lonestar6/round8_results_light.npz \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round8_qe75_patched_constraint_lonestar6/

scp \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/round8_qe75_patched_constraint_lonestar6/round8_dashboard.png \
  /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear/round8_qe75_patched_constraint_lonestar6/
```

Do not commit `qe_tmp`, full `.out` files, or scheduler logs. The local
`.gitignore` excludes these heavy/restart artifacts. Do not run MD or increase
the atom count until `target_print` is near `2.0`, the SCF chain converges, and
the final per-atom local moment is preferably at least `1 Bohr magneton`.

## Completed results (2026-08-15)

All 12 staged SCFs completed with no QE routine or Cholesky errors. Each output
contained one nonfatal `c_bands` warning, but all stages reached their requested
threshold. Every stage printed the corrected `2.0 Bohr magneton/Fe` target,
proving that the private patched executable—not the affected shared QE
implementation—was used.

| Final lambda (Ry) | Final-stage iterations | Accuracy (Ry) | Measured local moment (Bohr magneton/Fe) | Retention |
|---:|---:|---:|---:|---:|
| 0.01 | 9 | 2.7109e-4 | 0.626003 | 31.30% |
| 0.02 | 17 | 2.8997e-4 | 0.930252 | 46.51% |
| 0.05 | 15 | 2.9334e-4 | 1.323047 | 66.15% |
| 0.10 | 18 | 1.8290e-4 | 1.565806 | 78.29% |

The retained moment increases monotonically with lambda. Both `0.05` and
`0.10 Ry` meet the static acceptance threshold of at least
`1 Bohr magneton/Fe`; `0.10 Ry` gives the largest moment and best final
accuracy without a meaningful iteration penalty. Static SCF alone cannot show
whether the stronger penalty remains stable and affordable under ionic motion,
so Round 9 compares these two values in otherwise identical 10-step 4000 K MD
pilots.

Repository results retained here are deliberately lightweight:

- `round8_final_summary.csv`: human-readable terminal values;
- `round8_dashboard.png`: convergence and moment comparison;
- `round8_results_light.npz`: safe `allow_pickle=False` archive;
- `round8_serial_status.tsv`: Lonestar6 case completion record.
