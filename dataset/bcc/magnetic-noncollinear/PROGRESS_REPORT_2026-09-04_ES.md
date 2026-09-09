# Informe de progreso: Fe BCC paramagnético no colineal

**Fecha:** 4 de septiembre de 2026  
**Estado:** metodología aún en validación; todavía no hay una trayectoria de
producción apta para un ajuste TDEP cuantitativo.

## Objetivo

Construir dinámica molecular *ab initio* de Fe BCC a 4000 K con desorden
magnético local (DLM), momentos no colineales restringidos y fuerzas lo
suficientemente convergidas para producir configuraciones útiles para TDEP y
el cálculo posterior de fonones.

## Lo que hemos establecido

### 1. La implementación de la restricción magnética importa

Cada átomo necesita su propia especie artificial de QE para imponer una
dirección magnética independiente. Por ello, los cálculos de 16, 32 y 48
átomos requieren respectivamente `ntyp=16`, `32` y `48`, aunque todas las
especies utilicen el mismo pseudopotencial de Fe.

La versión estándar de QE 7.5 no es suficiente para este protocolo. El
ejecutable debe incluir:

- `ntypx=128` para aceptar más de diez especies;
- la corrección de la restricción atómica
  `mcons = Zval * starting_magnetization * direction`.

Se identificó un ejecutable antiguo que imponía solamente `0.125 mu_B` en vez
de los `2.0 mu_B` solicitados. Los resultados creados con ese binario no son
compatibles con los cálculos corregidos. El ejecutable reparado de Jakar fue
validado con 16 de 16 normas de restricción iguales a `2.0 mu_B`.

### 2. La construcción DLM está balanceada

Las texturas magnéticas utilizan direcciones aleatorias junto con sus
antípodas. La suma vectorial de los momentos objetivo es cero por construcción,
sin obligar a que todos los átomos compartan una misma dirección. Las
velocidades se generan o transforman de forma reproducible, se elimina el
movimiento del centro de masa y se reescalan exactamente a 4000 K.

Este protocolo representa MD con momentos DLM restringidos; no constituye una
simulación explícita de dinámica de espines.

### 3. Una restricción fuerte conserva bien las direcciones, pero dificulta el SCF

Con `lambda=0.100 Ry`, Round 20 conservó muy bien las direcciones objetivo:

- error angular medio: `0.66--0.85 grados`;
- percentil 95: `1.38--1.78 grados`;
- máximo observado: `2.62 grados`.

Sin embargo, esa precisión magnética introduce rigidez electrónica. Round 20
produjo 657 advertencias de corrección SCF grande sobre 786 pasos acumulados,
y la réplica difícil `dlm01` falló en el paso 86 después de alcanzar 500
iteraciones electrónicas. Por esta razón no se justificó extender las
trayectorias ni iniciar todavía una producción grande.

### 4. Las fuerzas requieren un umbral SCF más estricto

Round 21 aisló el problema usando la misma configuración difícil y comparó
tres combinaciones:

| Caso | `conv_thr` (Ry) | `mixing_beta` | Resultado observado |
|---|---:|---:|---|
| `conv3e-4_beta1e-2` | `3e-4` | `0.010` | Falló en el paso 67; 12 advertencias |
| `conv3e-4_beta5e-3` | `3e-4` | `0.005` | Falló en el paso 71; 11 advertencias y mayor costo SCF |
| `conv1e-4_beta1e-2` | `1e-4` | `0.010` | Seguía activo en el paso 46; cero advertencias |

La evidencia disponible rechaza `conv_thr=3e-4 Ry` para esta semilla difícil.
Reducir `mixing_beta` de `0.010` a `0.005` tampoco resolvió el problema y
aumentó el número de iteraciones. Hasta ahora, la mejor combinación es:

```text
conv_thr          = 1.0e-4 Ry
mixing_mode       = local-TF
mixing_beta       = 0.010
mixing_ndim       = 20
electron_maxstep  = 1000
```

El caso activo había usado 3224 iteraciones SCF en 46 pasos, unas 70 por paso.
Ese costo es alto pero el comportamiento de las fuerzas es más limpio. La
temperatura media observada era 4847 K y la instantánea 4463 K; por lo tanto,
la trayectoria aún no podía declararse térmicamente equilibrada. La energía,
la presión y las fuerzas también continuaban evolucionando.

### 5. Una trayectoria parcial no basta para fonones cuantitativos

El intento Round 19 de 48 átomos produjo solamente 59 configuraciones
(`0.059 ps`) y además utilizó la escala de restricción incorrecta. Presentó
advertencias de corrección de fuerza en 57 de 59 pasos. Sus dispersiones HELD
son útiles únicamente como diagnóstico visual; no constituyen fonones
validados ni una base defendible para TDEP.

La lección general es que terminar algunos pasos de MD no es suficiente. Para
TDEP se necesitan fuerzas convergidas, una restricción magnética verificada,
una trayectoria equilibrada y un número adecuado de configuraciones
decorrelacionadas.

## Nuevo experimento preparado: Round 22

Se creó un piloto de transferencia a 32 átomos para comprobar si la receta
electrónica de Round 21 sigue siendo estable al duplicar el sistema:

- supercelda BCC convencional `4x2x2`, 32 átomos;
- celda `9.56 x 4.78 x 4.78 A` y el mismo volumen por átomo;
- 32 especies y 32 direcciones magnéticas distintas;
- suma vectorial magnética objetivo numéricamente cero;
- velocidades con centro de masa cero y temperatura exacta de 4000 K;
- rampa `lambda=0.005 -> 0.020 -> 0.100 Ry`;
- SCF y MD finales con `conv_thr=1e-4 Ry`, `mixing_beta=0.010`,
  `electron_maxstep=1000` y `nraise=5`;
- piloto de 100 pasos MD en Gamma y 32 tareas MPI.

Las validaciones locales confirmaron 32 posiciones y especies únicas, balance
DLM, velocidad del centro de masa nula, temperatura inicial de 4000 K y la
presencia correcta de todas las opciones críticas en los archivos de entrada.

Esta supercelda alargada es apropiada como prueba intermedia de escalamiento,
pero no es la geometría final ideal para estudiar efectos de tamaño isotrópicos.
Si funciona, una celda más isotrópica de 54 átomos (`3x3x3` convencional) será
una mejor candidata física para la siguiente etapa.

## Criterio de aceptación

Round 21 debe terminar antes de declarar definitiva su receta. Round 22 debe
considerarse exitoso solamente si cumple simultáneamente:

1. 100 de 100 pasos MD y ningún ciclo SCF no convergido;
2. advertencias de corrección SCF grande en no más del 5% de los pasos;
3. error angular medio no mayor de 3 grados y percentil 95 no mayor de 7 grados;
4. ausencia de colapso sistemático de los momentos locales;
5. temperatura controlada alrededor de 4000 K después del transitorio;
6. costo medio de iteraciones SCF compatible con una trayectoria de producción.

Si la convergencia sigue fallando, la siguiente variable controlada debe ser
`lambda=0.050 Ry`, manteniendo `conv_thr=1e-4 Ry` y `mixing_beta=0.010`. No se
debe relajar nuevamente el umbral de fuerzas para obtener una trayectoria más
larga a costa de su calidad.

## Conclusión actual

El principal cuello de botella no es generar una celda grande, sino obtener
fuerzas confiables bajo restricciones magnéticas fuertes. Ya se corrigió el
problema del binario, se validó una construcción DLM balanceada y se identificó
que `conv_thr=1e-4 Ry` con `mixing_beta=0.010` es la opción más prometedora.
Round 22 constituye la prueba mínima necesaria para saber si esa mejora se
transfiere de 16 a 32 átomos. Aún no debe iniciarse una producción larga ni un
ajuste TDEP cuantitativo.

## Rutas reproducibles en el repositorio

- Round 20: `dataset/bcc/magnetic-noncollinear/round20_16atom_dlm8_conv1e3_nraise5_md100_jakar/`
- Round 21: `dataset/bcc/magnetic-noncollinear/round21_dlm01_scf_stabilization_jakar/`
- Round 22: `dataset/bcc/magnetic-noncollinear/round22_32atom_4x2x2_round21settings_md100_jakar/`
- Auditoría de 48 átomos: `dataset/bcc/magnetic-noncollinear/round19_48atom_6x2x2_qe75_ntyp128_lambda0p100_md1000_jakar/`
- Documentación del QE reparado: `hpc/jakar/README.md`
