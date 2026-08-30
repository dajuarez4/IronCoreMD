# Round 18: 32 atoms and 32 spin species on Lonestar6 vm-small

This round scales the Round-17 setup from 16 to 32 atoms while preserving its
working QE 7.5 parameters. It requires the private Lonestar6 executable built
with the corrected atomic constraint and `ntypx=128`:

```text
$SCRATCH/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x
```

## Physical construction

- BCC Fe `4x2x2` supercell, 32 atoms, cell `9.56 x 4.78 x 4.78 A`.
- The Round-17 `2x2x2` initial state is tiled once along x.
- One artificial Fe species per atom: `Fe01` through `Fe32`.
- 32 distinct noncollinear target directions in four antipodal groups.
- The target spin-vector sum is zero by construction.
- Round-17 velocities are deterministically permuted for the second tile,
  center-of-mass motion is removed, and the full set is rescaled to 4000 K.
- All species use the same Fe PAW pseudopotential.

## QE recipe

- SCF ramp: `lambda=0.005, conv_thr=4e-3 Ry`;
  `lambda=0.020, conv_thr=1e-3 Ry`;
  `lambda=0.100, conv_thr=3e-4 Ry`.
- `mixing_mode='local-TF'`, `mixing_beta=0.010`, `mixing_ndim=20`.
- CG diagonalization, Gamma point, 71/496 Ry cutoffs.
- MD: 1000 steps, `dt=20.67 au`, SVR thermostat at 4000 K,
  `nraise=20`, `lambda=0.100`, `conv_thr=4e-3 Ry`.

Unlike Rounds 12--16, this calculation does not repeat only eight magnetic
species and does not add an unconstrained bootstrap. Its purpose is to test
the direct, successful Round-17 ramp at the next system size using 32 truly
distinct species. Because this doubles the Round-17 workload, a 48-hour
`vm-small` job may produce only a partial trajectory; the output remains
useful and can be continued in a new restart workflow if necessary.

## Copy and submit

From the Mac:

```bash
cd /Users/dajuarez4/Documents/Fe/IronCoreMD/dataset/bcc/magnetic-noncollinear
scp -r round18_32atom_4x2x2_qe75_ntyp128_lambda0p100_md1000_lonestar6_vmsmall \
  dajuarez4@ls6.tacc.utexas.edu:/work/11381/dajuarez4/ls6/
```

On Lonestar6:

```bash
cd /work/11381/dajuarez4/ls6/round18_32atom_4x2x2_qe75_ntyp128_lambda0p100_md1000_lonestar6_vmsmall
test -x "$SCRATCH/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x" && echo EXISTS
sbatch run_round18_lonestar6_vmsmall.sbatch
```

Monitor with:

```bash
squeue -u "$USER"
bash check_round18_lonestar6.sh
tail -f JOBID.log
```

Do not submit from a directory containing old `*.out` files or `qe_tmp` state.
