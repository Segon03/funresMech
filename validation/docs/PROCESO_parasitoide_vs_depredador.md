# Proceso parasitoide frente a depredador en D2 y D3 (exploratorio)

Fecha: 2026-10-19 · Script: `validation/R/17_parasitoide_vs_depredador.R` · funresMech 1.2.0-rc2 · preset exploratory, 2 semillas, grilla de z 0.5–5 (paso 0.5, con refinamiento). Tablas: `validation/tablas/proceso_ajustes.csv`, `proceso_comparacion.csv`.

- **Parasitoide:** el hospedador atacado sigue disponible. Si se reencuentra, la hembra gasta tiempo de manipulación en él (sin discriminación).
- **Depredador:** el atacado se retira. Equivale a una discriminación perfecta y sin costo de los hospedadores ya parasitados.

Ambos modelos tienen 5 parámetros, así que el ΔAIC es 2 × ΔNLL.

## Resultados

| Datos | ΔAIC (depredador − parasitoide), semillas 1 y 2 | z (parasitoide) | IC de z (parasitoide) | z (depredador) | IC de z (depredador) |
|---|---|---|---|---|---|
| D2 (73 ensayos) | +6.3, +6.4 | 2.5, 3.0 | [1.09, abierto], [1.29, abierto] | 1.0, 1.5 | [0.86, 2.28], [0.86, 2.57] |
| D2 sin el ensayo atípico (72) | +5.5, +5.6 | 2.0, 2.5 | [1.60, 3.95], [1.51, 3.73] | 1.0, 1.25 | [0.83, 1.67], [0.87, 1.76] |
| D3 (64) | +5.8, +4.6 | 2.5, 2.5 | [2.05, 4.22], [1.81, abierto] | 1.2, 1.25 | [0.81, 1.61], [0.78, 1.79] |

- h es similar en los dos procesos (0.58–0.62 h en el parasitoide y 0.59–0.76 h en el depredador).
- k es muy chico en D2 con el parasitoide (0.03–0.06), y el depredador lo lleva a 0.2–0.3.

## Lectura

1. **El proceso parasitoide ajusta mejor en los tres conjuntos y en las dos semillas** (ΔAIC entre +4.6 y +6.4, 6 de 6 en el mismo sentido). La diferencia de NLL (2.3–3.2) supera varias veces el ruido Monte Carlo del óptimo (0.1–0.6), así que el sentido es robusto. La magnitud es aproximada (preset exploratorio).
2. **No depende del ensayo atípico de D2**: sin él, el ΔAIC sigue en +5.5.
3. **La conclusión sobre z depende del proceso supuesto.** Con el parasitoide, z ≈ 2–3 y el IC excluye z = 1 (respuesta de tipo III). Con el depredador, z ≈ 1–1.25 y el IC contiene a 1 (tipo II) en los tres conjuntos. Es la misma conclusión del brazo B de la Fase 5: el proceso pesa más que la familia de distribución.
4. **El depredador no queda descartado del todo.** Un ΔAIC de 5–6 equivale a un peso relativo de ~5–8% (exp(−Δ/2)). Además, la discriminación real puede ser parcial o tener costo, un caso intermedio que el paquete no modela.

## Implicación biológica (texto sugerido para la discusión)

Los datos de *Apanteles myeloenta* son más compatibles con un parasitoide que vuelve a encontrar y manipula hospedadores ya parasitados que con uno que los evita por completo. Eso **no prueba la ausencia de discriminación**: el modelo depredador solo representa el extremo de una discriminación perfecta y sin costo. La evidencia de una respuesta de tipo III (z > 1) está condicionada a ese supuesto de proceso, y debería contrastarse con observaciones directas de superparasitismo, como la disección de hospedadores o el registro del comportamiento de rechazo.

## Si se informa en el texto principal

Repetir con el preset standard (unas 4–5 h para los 6 ajustes en la PC de 23 núcleos). Para el material suplementario alcanza con esta corrida.
