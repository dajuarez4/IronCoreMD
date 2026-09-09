# Round 22: 32-atom transfer pilot on Jakar

This pilot tests whether the cleanest Round-21 electronic settings transfer
from 16 to 32 atoms before committing to a long production trajectory.

- BCC Fe `4x2x2` conventional supercell: 32 atoms, `9.56 x 4.78 x 4.78 A`.
- Round-21 difficult-seed positions are tiled once along x.
- One distinct species and noncollinear target direction per atom; the target
  spin-vector sum is zero by construction.
- Velocities are decorrelated between tiles, center-of-mass motion is removed,
  and the 32-atom set is rescaled to 4000 K.
- Staged constraint ramp: `lambda=0.005`, `0.020`, then `0.100 Ry`.
- Final SCF and MD use `conv_thr=1e-4 Ry`, `mixing_beta=0.010`,
  `mixing_mode='local-TF'`, `electron_maxstep=1000`, and `nraise=5`.
- The pilot is 100 MD steps at Gamma with 32 MPI tasks.

The Jakar executable must contain both the corrected atomic constraint and
`ntypx=128`: `/gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x`.

Generate and validate locally:

```bash
python generate_round22.py
bash check_round22.sh
```

Copy this directory to Jakar, then submit with:

```bash
sbatch run_round22_jakar.sbatch
```

Advance toward production only if all 100 steps finish, no SCF cycle is
unconverged, force-correction warnings occur on no more than 5% of steps, the
temperature is controlled near 4000 K, and the SCF cost remains affordable.
The elongated 4x2x2 cell is an intermediate scaling test, not the preferred
final cell for isotropic finite-size studies.
