# Round 23: Round 22 con nraise=20

Copia de las entradas iniciales de round22. Único cambio físico: `nraise=5` → `nraise=20`.
El prefijo QE cambia a `Fe_bcc32_r23_conv1e4_beta1e2` para identificar esta corrida.

Se conservan las posiciones, velocidades, 32 direcciones de espín, pseudopotencial,
celda BCC 4×2×2 de 32 átomos, 4000 K, termostato SVR, dt=20.67,
100 pasos MD y las tres etapas SCF de preparación. Inicio desde las entradas
originales, no desde el reinicio del paso 63 de round22.

En Jakar, desde este directorio:

```bash
sbatch run_round23_jakar.sbatch
bash check_round23.sh
```

Recursos: 1 nodo, 32 procesos MPI, 48 horas, cuenta jakar_general.
Ejecutable: /gpfs/scratch/dajuarez4/qe-7.5/bin/pw.x (mismo QE modificado que round22).

Comparar los resultados con round22 para evaluar el efecto sobre los fonones;
el cambio de nraise por sí solo no demuestra una mejora.
