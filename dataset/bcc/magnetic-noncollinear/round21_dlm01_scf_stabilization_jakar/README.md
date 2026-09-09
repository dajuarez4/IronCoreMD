# Round 21: SCF-stabilization pilots on the difficult Round-20 seed

Round 20 retained stable local moments but produced 657 large-SCF-force-correction
warnings over 786 steps, and `dlm01` failed at step 86 after 500 electronic
iterations. Longer production MD is therefore deferred.

This controlled test keeps the Round-20 `dlm01` spin texture, positions,
velocities, thermostat, timestep, and constraint fixed. It compares:

1. `conv_thr=3e-4 Ry`, `mixing_beta=0.010`, 800 maximum iterations.
2. `conv_thr=1e-4 Ry`, `mixing_beta=0.010`, 1000 maximum iterations.
3. `conv_thr=3e-4 Ry`, `mixing_beta=0.005`, 1000 maximum iterations.

Generate locally with `python generate_round21.py`. On Jakar, submit with
`sbatch run_round21_array_jakar.sbatch` and monitor with `bash check_round21.sh`.

Select a production setting only if all 100 steps complete, force-correction
warnings fall sharply, and the mean electronic iterations per MD step remains
affordable. Do not submit longer trajectories solely because a case finishes.
