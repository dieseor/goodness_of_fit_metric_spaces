# TODO EN DUDA — Dispatcher para formas cuadráticas no centrales

Estado a 2026-07-26: el candidato se ha implementado como infraestructura
compartida, marcado `EN DUDA`, por petición explícita del autor. Sigue pendiente
la revisión matemática final y su validación formal antes de considerar la
decisión cerrada. No se han modificado el artículo ni resultados ya guardados.

Para

\[
Q=\sum_j\lambda_j\chi^2_{h_j}(\delta_j),\qquad \lambda_j>0,
\]

la cadena candidata es:

1. Aplicar el pre-screen de Farebrother
   \(\kappa(\lambda)>10^4\) o
   \(\log a_0<\log(\texttt{double.xmin})\), una vez validado en parámetros
   empíricos. Si se activa, no ejecutar Farebrother.
2. Rama con pre-screen: Imhof a precisión habitual; después aceptar la misma
   salida si ya satisface la precisión de rescate \(10^{-4}\), o reintentar
   Imhof con ese objetivo; después Davies; finalmente Monte Carlo.
3. Rama sin pre-screen: Farebrother; si `ifault != 0`, Imhof habitual;
   después el mismo rescate Imhof, Davies y Monte Carlo. No reintentar
   Farebrother tras haberlo descartado preventivamente.
4. Aceptar Imhof con recorte sobre \([0,1]\) sólo cuando `abserr <= eps` y
   la salida cruda pertenece a \([-abserr, 1 + abserr]\). No recortar una
   salida con error grande.
5. Davies sólo se acepta con `ifault = 0`. Para el rescate propuesto,
   estudiar `acc = 10^{-4}` y un `lim` alto; validar el compromiso coste /
   precisión antes de fijar valores.
6. Si se llega a Monte Carlo, imprimir un mensaje al inicio y al final de cada
   rutina, incluyendo identificador de caso, tolerancia, número de muestras,
   error/intervalo alcanzado y tiempo. La rutina debe informar también de una
   interrupción.
7. No usar HBE como fallback de seguridad.

Evidencia exploratoria, no concluyente: en cuatro casos de estrés que no
alcanzaron \(10^{-4}\) con el primer Imhof \((10^{-8},\ \texttt{limit}=20000)\),
una segunda llamada con objetivo \(10^{-4}\) rescató dos. Elevar `limit` a
200000 no modificó los dos restantes. Por tanto, el reintento Imhof puede
compensar, pero no hay evidencia de que aumentar `limit` sea una solución
general.
