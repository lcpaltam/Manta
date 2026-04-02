# Ajuste intermedio aplicado

Se generó una versión intermedia del modelo a partir de `mod_courtemanche_WTsync.f90` con estos cambios:

| Parámetro / tasa | Antes (`WTsync`) | Ahora (`tuned_v1`) | Efecto buscado |
|---|---:|---:|---|
| `p_gNa` en `get_parameter_CRN` | `p_gNa` (desde `.inc`, típicamente 7.8) | `15.0` | Aumentar `INa_peak` y `dV/dtmax` |
| `actshift` | `0.0 mV` | `-20.0 mV` | Facilitar activación/reclutamiento de canales |
| `a2` | `1.0 * WT` | `0.25 * WT` | Reducir salida rápida del estado abierto |

Notas:
- Se mantiene el control por concentración (`p_verna_conc_nM=0` para control sin fármaco).
- No se tocaron las demás tasas WT ya sincronizadas en el bloque Markov.
