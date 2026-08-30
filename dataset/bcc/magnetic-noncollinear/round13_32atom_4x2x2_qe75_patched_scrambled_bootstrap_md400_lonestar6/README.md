# Round 13: corrected 32-atom BCC noncollinear MD400

Round 13 replaces the failed Round-12 initialization while keeping the same
32-atom `4x2x2` thermal structure, target temperature, velocities, magnetic
directions, final constraint strength, and 400-step MD protocol.

## Corrections relative to Round 12

- Magnetic labels are deterministically shuffled with seed `130032`.
- Each BCC sublattice contains two copies of every one of the eight magnetic
  directions. Each sublattice and the full cell therefore have zero target
  moment vector sum; antipodal directions are no longer segregated by
  sublattice.
- A weak bootstrap stage at `lambda=0.001 Ry`, `degauss=0.05 Ry`, and
  `conv_thr=1e-2 Ry` is run before the established lambda ramp.
- Density mixing uses bulk-appropriate `plain` mixing, `mixing_beta=0.01`, and
  a longer 40-vector history.

The staged path is:

1. `bootstrap_l0p001_1e-2`: lambda 0.001, degauss 0.05, conv_thr 1e-2 Ry;
2. `stage1_l0p005_4e-3`: lambda 0.005, degauss 0.05, conv_thr 4e-3 Ry;
3. `stage2_l0p020_1e-3`: lambda 0.020, degauss 0.025334, conv_thr 1e-3 Ry;
4. `stage3_l0p100_3e-4`: lambda 0.100, degauss 0.025334, conv_thr 3e-4 Ry;
5. 400 MD steps at lambda 0.100 and degauss 0.025334 Ry.

All stages must converge. The submission script stops at the first failed SCF.
Do not reuse Round 12's `qe_tmp`; Round 13 starts from atomic charge and wave
functions in its own directory.

## Generate and submit

Copy the prepared package to Lonestar6:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear
scp -r round13_32atom_4x2x2_qe75_patched_scrambled_bootstrap_md400_lonestar6 \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/
```

Then submit on Lonestar6:

```bash
cd /work/11381/dajuarez4/ls6/round13_32atom_4x2x2_qe75_patched_scrambled_bootstrap_md400_lonestar6
mkdir -p logs
sbatch run_round13_lonestar6.sbatch
```

The checked-in inputs are already generated. Local regeneration is optional;
the generator refuses to overwrite calculation outputs or `qe_tmp`.

Monitor with:

```bash
bash check_round13.sh
bash plot_round13.sh --live 30
```

After a successful trajectory:

```bash
python3 collect_round13_npz.py --require-complete
```
