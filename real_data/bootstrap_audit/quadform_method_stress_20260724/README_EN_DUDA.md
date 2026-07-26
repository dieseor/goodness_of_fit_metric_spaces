# EN DUDA — 2026-07-24

Esta batería es una auditoría numérica y matemática de los evaluadores de

\[
Q=\sum_j \lambda_j\chi^2_{h_j}(\delta_j),\qquad \lambda_j>0.
\]

No modifica el dispatcher de producción, los experimentos de la Sección 6 ni
ningún resultado del artículo. Sus objetivos son:

1. comprobar si \(\log a_0\) y
   \(\rho=1-\lambda_{\min}/\lambda_{\max}\) explican el coste y los
   `ifault` de Farebrother;
2. comparar Farebrother, Davies e Imhof con tolerancias declaradas;
3. usar Monte Carlo sólo como contraste independiente para casos seleccionados
   con discrepancias macroscópicas; no como referencia de precisión fina.

El script `scripts/audit_quadform_stratified_timing.R` y todos los resultados
de `stratified_48_cases/` y
`stratified_144_cases_preventive_rules/` también quedan explícitamente
**EN DUDA**. En particular, la clasificación de este último directorio mide
eventos operativos (`ifault` de Farebrother y desacuerdos cuando Imhof satisface
su propio `abserr`); no identifica por sí sola la CDF matemática verdadera.

También queda explícitamente **EN DUDA** la batería
`preventive_432_analytic_rules/`. Para esa batería concreta, el único evento
llamado “problema” en las tablas de reglas preventivas es `ifault != 0` de
Farebrother; no se usa Imhof para esa etiqueta. La rejilla tiene 432 celdas de
estrés deterministas y no representa todavía una distribución de parámetros
ajustados en los experimentos empíricos.

No se adoptará ninguna regla de selección sin revisión matemática y aprobación
explícita del autor.
