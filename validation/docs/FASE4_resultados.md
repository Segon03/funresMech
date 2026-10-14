# Fase 4: comparación de optimizadores y estimación de z (2026-10-13)

Corrida: `validation/R/11_comparacion_optimizadores.R` en la PC de 23-24
núcleos; 2.73 h de reloj (56 h-núcleo; ~190 s por punto, unas 2 veces más
rápido por núcleo que la PC del piloto). Presupuesto de cada ajuste: preset
exploratory (2000 simulaciones por densidad, 50 x 40 evaluaciones).

## Verificación de las salidas

* Parte A: 540 filas (9 escenarios x 5 datos x 4 estrategias x 3 corridas),
  sin faltantes ni duplicados; la brecha y el éxito recalculados coinciden
  (diferencias < 1e-12, redondeo del CSV).
* Parte B: 54 filas (9 x 3 x 2); la cobertura recalculada coincide.

## A. Estrategias de optimización (mismo presupuesto)

Brecha = NLL precisa (20000 simulaciones, 5 reevaluaciones) menos la mejor
encontrada para el mismo conjunto de datos. Éxito = brecha < 0.5.

| Estrategia | Éxito (< 0.5) | < 1 | < 2 | Brecha media | Percentil 95 | Máxima | En cota | s por ajuste |
|---|---|---|---|---|---|---|---|---|
| DEoptim (1.1.0) | 43% | 75% | 95% | 4.90 | 2.06 | 423 | 0% | 198 |
| CMA-ES | 78% | 86% | 92% | 2.51 | 14.4 | 68 | 1.5% | 185 |
| CRS2 | 67% | 70% | 72% | 10.5 | 59.1 | 280 | 16% | 164 |
| DEoptim + CRN + pulido | **80%** | **96%** | **99%** | **0.36** | **0.94** | **2.2** | 0% | 209 |

* DEoptim con números aleatorios comunes y pulido Nelder-Mead es la más
  robusta: nunca falla de forma grave y es la más reproducible entre corridas
  (rango mediano de la NLL entre corridas 0.36 frente a 0.82 de DEoptim).
  Cuesta 6% más tiempo.
* CMA-ES encuentra a menudo el mejor óptimo, pero tiene fallas graves
  ocasionales (percentil 95 = 14 unidades de NLL).
* CRS2 se descarta: colas largas y 16% de estimaciones en una cota.
* Incluso la mejor estrategia queda a más de 0.5 unidades del mejor óptimo en
  el 20% de las corridas con el presupuesto exploratorio: la verificación
  (corridas independientes, más esfuerzo) es necesaria.

## B. Estimación de z: perfil frente a conjunta

| Estrategia | Cobertura IC 95% | IC abierto | Ancho mediano | Excluye z = 1 si z = 1 | Excluye z = 1 si z > 1 | h-núcleo |
|---|---|---|---|---|---|---|
| Perfil (grilla 0.5-5, paso 0.5) | 26/27 (96%) | 11% | 1.15 | 0% | 100% | 14.9 |
| Conjunta + perfil adaptativo | 18/27 (67%) | 4% | 0.73 | 11% | 94% | 12.8 |

* La estimación conjunta subcubre: con 5 parámetros libres y el mismo
  presupuesto, los puntos vecinos del perfil quedan peor optimizados que el
  punto conjunto y el intervalo resulta artificialmente estrecho. Coincide con
  Okuyama (2026): la optimización es más estable con menos parámetros libres.
* El chequeo de fiabilidad detectó 4 de los 9 intervalos conjuntos que no
  cubren (falla) y el único intervalo de perfil que no cubre (falla).
* Advertencia: la grilla del perfil contiene los valores verdaderos de z, lo
  que favorece al perfil en el error de z; la comparación justa es la
  cobertura.

## Decisiones propuestas

1. Método por defecto: DEoptim + números aleatorios comunes + pulido.
2. z por perfil sobre grilla como método por defecto; la estimación conjunta
   queda como opción experimental, documentando su subcobertura.
3. La Fase 5 evalúa también el preset standard, porque con el exploratorio
   el 20% de los ajustes no llega al mejor óptimo.
