# Input-parameter comparison

| Parameter | try1 | try2 |
|---|---:|---:|
| Starting magnetization (all 128 types) | 0.15 | 0.05 |
| mixing_beta | 0.100 | 0.010 |
| mixing_mode | TF | local-TF |
| mixing_ndim | 12 | 20 |
| diago_thr_init | QE default | 1.0e-6 |
| diago_cg_maxiter | QE default | 200 |
| calculation / restart | md / from_scratch | md / from_scratch |
| nat / ntyp | 128 / 128 | 128 / 128 |
| ecutwfc / ecutrho (Ry) | 100 / 496 | 100 / 496 |
| occupations / smearing / degauss (Ry) | smearing / m-v / 0.04 | smearing / m-v / 0.04 |
| noncolin / lspinorb / nosym | true / false / true | true / false / true |
| constraint / lambda (Ry) | atomic direction / 0.2 | atomic direction / 0.2 |
| conv_thr (Ry) / electron_maxstep | 1.0e-4 / 500 | 1.0e-4 / 500 |
| diagonalization | cg | cg |
| MD nstep / dt (a.u.; fs) | 400 / 20.670; 0.500 fs | 400 / 20.670; 0.500 fs |
| thermostat / T / nraise | SVR / 5000 K / 20 | SVR / 5000 K / 20 |
| potential / wavefunction extrapolation | second_order / second_order | second_order / second_order |
| cell vectors (Angstrom) | (8.56,0,0); (4.28,7.41318,0); (0,0,13.68) | same |
| k-point grid / shift | 1x1x1 / 0 0 0 | 1x1x1 / 0 0 0 |
| pseudopotential | Fe PBE-spn KJPAW PSL 1.0.0 | same |

The 128 polar and azimuthal spin-direction angles are identical in both inputs; only their common starting magnitude differs.
