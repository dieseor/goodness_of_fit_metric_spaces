# EN DUDA — 2026-07-24

Los scripts y resultados de `bootstrap_audit/` son diagnósticos exploratorios.
No modifican los resultados del paper ni establecen todavía una regla de
selección entre HBE, Farebrother, Imhof y Davies. Cualquier uso en producción
queda pendiente de revisión matemática y aprobación explícita del autor.

Quedan expresamente **EN DUDA** los scripts de auditoría:

- `scripts/audit_ks_cvm_bootstrap_real_data.R`;
- `scripts/audit_fitted_parameter_stability.R`;
- `scripts/audit_joint_ks_cvm_percentile_gaps.R`;
- `scripts/audit_logistic_gaussian_quadform_real_profiles.R`;
- `scripts/refine_logistic_gaussian_quadform_disagreements.R`;
- `scripts/analyze_hbe_quadform_error_features.R`.

Y quedan **EN DUDA** todos los CSV, RDS y figuras bajo este directorio. Son
evidencia para revisar, no resultados que deban incorporarse a las tablas del
artículo ni parámetros de una regla automática.

TODO HBE: recuperar y verificar el teorema de cota de Buckley--Eagleson,
comprobar exactamente sus hipótesis frente a la suma no central usada aquí y
decidir sólo entonces si existe una condición automatizable para emplear HBE.
