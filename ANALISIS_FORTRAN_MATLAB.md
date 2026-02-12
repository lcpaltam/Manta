# Comparación Fortran (`mod_courtemanche.f90`) vs MATLAB (`Courtemanche_model_v2.m`)

## Diferencias relevantes en código y parámetros

1. **Concentración de fármaco en Markov de INa**
   - MATLAB toma `settings.drug_compound` y `settings.drug_concentration` dinámicamente por corrida.
   - Fortran usa el parámetro fijo `p_verna_conc_nM` (actualmente en `0.0`) para construir `conc` en el bloque Markov.

2. **Escalado de conductancia de sodio (`gNa`)**
   - MATLAB HH usa `g_Na = 7.8`.
   - Fortran define `p_gNa = 7.8*0.58`.

3. **Integración numérica del bloque Markov (INa)**
   - MATLAB delega la rigidez al solver ODE.
   - Fortran ahora usa substeps + actualización **semi-implícita diagonal** estado-a-estado:
     \[ X_{n+1} = \frac{X_n + \Delta t\,P_n}{1 + \Delta t\,L_n} \]
     con `dX/dt = P - L*X`.

## Mejor ajuste recomendado

1. Mantener el criterio adaptativo por rigidez (`max_rate`) para `D_conc > 0`.
2. Mantener protección de exponenciales (`safe_exp`) para evitar overflow.
3. Si se pierden picos de INa, subir gradualmente el piso de `n_sub` (por ejemplo 6→8→10) antes de subir mucho el máximo.
4. Para homologar Fortran vs MATLAB, alinear `gNa` y el manejo de concentración de fármaco en ambos lados.

## Cambios implementados en esta versión

- Se mantuvo `safe_exp` en las tasas Markov sensibles.
- Se reemplazó la actualización explícita por una **semi-implícita** para los 27 estados Markov (nativo, Dp, D).
- Se conservó renormalización y fallback de seguridad para NaN.
