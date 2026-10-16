# Fase 5: validación ampliada, verificación e interpretación

Corrida: `validation/R/12_validacion_ampliada.R`, PC de 23-24 núcleos, terminada el 2026-10-16 (679 de 679 trabajos, sin errores).
Análisis: `validation/R/13_analisis_fase5.R`, que genera las tablas `fase5_analisis_*.csv` y las figuras `figuras/fase5_*.png`.

Salvo que se indique otra cosa, todos los números de A a E son del preset **exploratory** (2000 simulaciones por densidad, 50 generaciones × 40 individuos) con grilla de z de 0.5 a 5 en pasos de 0.5.

---

## 1. Verificación (Paso 1)

| Chequeo | Resultado |
|---|---|
| 679 resultados guardados; conteos A 360, B 160, C 180, D 72 + 12, E 4680, F 3 | Pasa |
| C (exploratory) y D (1.2.0) son los mismos ajustes que A | Pasa (idénticos) |
| Todo intervalo contiene a z_hat | Pasa |
| En F, el mínimo del perfil es z_hat; D2 sin atípico = 72 ensayos, fila 13 quitada | Pasa |
| Falla de optimización (NLL precisa del estimado > NLL precisa en la verdad + 1) | 1.9% de los ajustes de A (7 de 360). Es plausible. |
| La NLL precisa en la verdad es estable | **Con reparos.** En 2 ajustes de A (y en una evaluación de la 1.1.0), la evaluación en la verdad tocó el piso de probabilidad: un resultado observado no apareció en 20 000 simulaciones y sumó +708. Eso explica los déficits muy negativos (−421 y −283) y el rango de 140 en C. No afecta las métricas basadas en cobertura, solo los promedios de déficit; en el análisis se usan medianas. |
| D3 con preset final frente a Okuyama (2026) y frente a la validación de la 1.1.0 | Pasa. z = 2.25 [1.45, 5.10]; con la 1.1.0, [1.4, 4.8]. |
| D2 sin la fila 13 frente a la auditoría del 2026-10-05 | Pasa. z = 2.25 [1.45, 4.04]; antes, 2.0 [1.40, 4.28]. |
| D2 completo | Pasa. z = 3.75 [1.18, abierto]: perfil plano hasta z = 20, como en la validación de la 1.1.0. |

**Conclusión:** los resultados son internamente coherentes y reproducen lo ya conocido para D2 y D3. Aparecieron dos problemas del método, que no son errores de la corrida y se documentan como hallazgos H16 y H17.

---

## 2. Interpretación por brazo (Paso 2)

### A. Simulación amplia (360 ajustes, modelo correcto)

| Medida | Valor | IC 95% |
|---|---|---|
| Cobertura del IC de z (nominal 95%) | **86.7%** | 82.8–89.8 |
| · z verdadero = 1 (en la grilla) | 91.7% | 85.3–95.4 |
| · z = 1.75 (a mitad de camino entre 1.5 y 2.0) | **76.7%** | 68.3–83.3 |
| · z = 2.75 (a mitad de camino entre 2.5 y 3.0) | 91.7% | 85.3–95.4 |
| Error de tipo I (IC excluye 1 cuando z = 1) | 6.7% | 3.4–12.6 |
| Potencia para detectar z > 1 (z = 1.75 / 2.75) | 85.8% / 100% | — |
| IC abierto | 2.8% (10 ajustes; todos con límite superior abierto, ninguno con el inferior abierto) | 1.5–5.0 |
| Falla de optimización (déficit > 1) | 1.9% | 0.9–4.0 |

- La mediana de z_hat coincide con el punto de grilla más cercano al valor verdadero: 1.0, 2.0 y 3.0.
- h se recupera sin sesgo (mediana de ĥ/h entre 0.95 y 1.08), salvo con pocos encuentros (E = 5: +22% a +32%).
- La extensión de la grilla hacia abajo se activó en 28 ajustes y ninguno quedó con el límite inferior abierto. En la Fase 4 eran 2 de 3 los intervalos abiertos por abajo.
- El diseño pequeño (30 ensayos) cubrió mejor (94%) que el grande (86%). Encaja con la causa descrita en H17: los intervalos más anchos sufren menos la grilla gruesa.

**La subcobertura es el resultado más importante de la Fase 5 y tiene dos causas, ambas corregibles en el código:**

**H16. Saltos por piso de probabilidad en el perfil.**
- En 65 de 3979 puntos de perfil (33 de 360 ajustes, 9%), alguna de las 3 reevaluaciones no simuló un resultado observado. Eso suma +708 por resultado y deja picos de cientos de unidades en el perfil.
- 13 de esos picos caen justo en z = 1, lo que infla el ΔAIC frente a z = 1: 3 de los 9 falsos positivos a z = 1 se deben a esto.
- El chequeo "Observed outcomes produced by the simulation" no los detecta, porque mira solo el óptimo.

**H17. Grilla gruesa e interpolación lineal.**
- El IC se interpola linealmente entre puntos separados 0.5. Cuando el perfil es curvo y el intervalo es angosto (ancho mediano ~0.9), la interpolación lineal acorta el intervalo.
- Es peor cuando el z verdadero cae a mitad de camino entre dos puntos de la grilla (z = 1.75).
- Hay casos degenerados, como [1.496, 1.502] o [2.496, 2.504]: un vecino de la grilla queda muy por encima del umbral y el límite cae pegado a z_hat.
- El refinamiento actual solo agrega puntos en tramos más anchos que 1.5 veces el paso, así que nunca refina la grilla base.

Recalculé el IC sobre **los mismos perfiles** (tabla `fase5_analisis_A_cobertura_alternativas.csv`):

| | Actual | Sin los puntos con piso | Sin piso y con spline monótona |
|---|---|---|---|
| Global | 86.7% | 87.8% | **90.6%** |
| z = 1.75 | 76.7% | 78.3% | 85.8% |

Esto no basta para llegar al 95%, porque la spline no puede recuperar la forma del perfil entre puntos separados 0.5. Hacen falta **más puntos alrededor del mínimo y de los cruces del umbral**. Eso requiere volver a ajustar.

**El presupuesto no es la causa** (brazo C): con el preset standard, la cobertura a z = 1.75 sigue en 58% (7 de 12). La subcobertura viene de cómo se construye el intervalo, no de la optimización.

### B. Desvíos del modelo (160 ajustes; k = 0.3; z = 1 y 2.25; 8 conjuntos por celda)

| Desvío | Ajuste por defecto | Ajuste con el modelo correcto | ¿El AIC elige el correcto? |
|---|---|---|---|
| Búsqueda Lognormal | z sin sesgo relevante (+0.12); cobertura 88–100% | Similar | 62% (ΔAIC > 2 en 38%) |
| Búsqueda Weibull | +0.25 a +0.31; cobertura 88% | Similar (75%) | 50% (ΔAIC > 2 en 6%) |
| Manipulación Gamma (s/h = 1) | +0.38; cobertura 88% | Similar | 75% (ΔAIC > 2 en 12%) |
| **Proceso depredador** | **z sesgado: +2.94 si z = 1, −0.44 si z = 2.25; h subestimado (ĥ/h = 0.66)** | Corregido: +0.25 y −0.06; ĥ/h = 1.01 | 88% (ΔAIC > 2 en 50%) |
| Heterogeneidad individual en a (σ = 0.5) | +0.50 a +0.56; cobertura 75–100% | No existe en el paquete | — |

Lectura:
- **La familia de distribuciones importa poco para z con estos diseños.** Los datos casi no la distinguen: el AIC acierta al nivel del azar en búsqueda.
- **El proceso sí importa.** Ajustar un proceso parasitoide a datos de depredador sesga mucho z y h. El modo depredador lo corrige y el AIC lo detecta la mitad de las veces.
- La heterogeneidad entre individuos sesga z hacia arriba alrededor de 0.5. Se declara como limitación.
- Con 8 conjuntos por celda, cada proporción tiene una incertidumbre de ±15 puntos. Las conclusiones son cualitativas.

### C. Presupuesto (36 conjuntos; núcleo, diseño grande)

| Variante | Éxito ±0.5 NLL | Éxito ±1 NLL | Brecha media | Cobertura | Tiempo mediano (1 núcleo) |
|---|---|---|---|---|---|
| exploratory (1 corrida) | 78–83% | 92–94% | 0.30–0.42 | 72–78% | ~35 min |
| mejor de 3 exploratory | 83% | **100%** | 0.25 | 75% | ~2.2 h |
| standard | **100%** | **100%** | 0.07 | 81% | ~2.7 h |

- Para la **optimización**, standard es claramente superior y cuesta poco más que el mejor de 3. **El mejor de 3 no aporta frente a standard.**
- La **cobertura no mejora con el presupuesto**, lo que confirma H17.

### D. 1.1.0 frente a 1.2.0

Con los mismos 36 conjuntos simulados (núcleo, diseño grande) y la misma grilla:

| | 1.1.0 | 1.2.0 |
|---|---|---|
| Cobertura | 80.6% | 77.8% |
| IC abiertos | **16.7%** | **0%** |
| Error absoluto de z (media) | 0.39 | 0.32 |
| Brecha a la mejor solución (NLL) | 0.33 | 0.12 |
| Detección de z > 1 | 61% | 61% |
| Tiempo (1 núcleo) | ~56 min (3000 simulaciones) | ~35 min de ajuste (2000 simulaciones); 58 min con `verify_fit()` |

Con D2 y D3 reales (3 semillas por versión, grilla con paso 0.25):
- Las estimaciones y las NLL precisas son equivalentes entre versiones: D3 z = 2.33 de media en las dos.
- El límite inferior de D3 es más estable en la 1.1.0 (rango 0.02 contra 0.06), pero en los dos casos es menor que el paso de la grilla.
- La 1.2.0 tardó entre 35% y 50% menos con D2 y D3, en parte porque usa 2000 simulaciones en vez de 3000.

**Lectura honesta:**
- La 1.2.0 optimiza mejor y elimina los intervalos abiertos (cotas adimensionales con ampliación automática, y extensión de la grilla hacia los dos lados).
- **No mejora la cobertura**, porque las dos versiones comparten la construcción del intervalo (H17).
- La mejora para el manuscrito existe, pero es moderada. La mejora importante tiene que venir de corregir H16 y H17.

### E. Chequeos de fiabilidad

| Chequeo | Alertas | AUC frente a no cubre | AUC frente a falla de optimización |
|---|---|---|---|
| Independent runs agree (NLL) | 55% | 0.65 | 0.67 |
| Monte Carlo noise of the NLL | 44% | 0.56 | 0.70 |
| Optimiser settled before the end | 54% | 0.51 | 0.51 |
| Optimiser value vs re-evaluated NLL | 9% | 0.51 | 0.69 |
| More effort does not improve | 3% | 0.50 | 0.70 |
| Profile of z regular | 5% | 0.56 | 0.70 |
| Semáforo global: alguna falla | 8% | sensibilidad 17%, especificidad 93% | sensibilidad 43%, especificidad 93% |
| Semáforo global: alguna alerta | 89% | sensibilidad 92%, especificidad 12% | sensibilidad 86%, especificidad 11% |

Lectura:
- Los chequeos **detectan moderadamente las fallas de optimización** (AUC ~0.7), aunque solo hubo 7, así que la estimación es imprecisa.
- **No pueden detectar la subcobertura**, y es lógico: no es una falla de optimización, es un problema de cómo se construye el intervalo. La solución es corregir la causa, no agregar chequeos.
- Con los umbrales provisorios, **tres chequeos dan alerta en la mitad de los ajustes**, casi sin discriminar. Eso es fatiga de alarmas: hay que subir esos umbrales.
- El chequeo del piso de probabilidad no mira los puntos fuera del óptimo (H16).
- Umbral sugerido, que deja un 10% de falsas alarmas entre los ajustes buenos (`p90_buenos`):
  - Independent runs agree (NLL): advertencia a 1.0 en vez de 0.5;
  - Monte Carlo noise: 0.4 en vez de 0.25;
  - mejora tardía: 3.4 en vez de 1;
  - brecha optimizador–reevaluación: 2 (queda igual).

  Antes de fijarlos hay que recalibrar con los datos de la corrida de confirmación (H18).

**H18.** Los umbrales provisorios generan demasiadas alertas, y la calibración sobre 7 fallas es pobre. Se recalibran después de corregir H16 y H17.

### F. D2 y D3 con preset final (material suplementario)

| Datos | Ensayos | z | IC 95% | k | ΔAIC frente a z = 1 | ΔAIC frente al beta-binomial | Semáforo |
|---|---|---|---|---|---|---|---|
| D2 | 73 | 3.75 | [1.18, abierto hasta 20] | 0.036 | 3.27 | **+8.72** (prefiere el mecanístico) | warning (IC abierto, k chico, a y h variables entre corridas) |
| D2 sin la fila 13 | 72 | 2.25 | [1.45, 4.04] | 0.51 | 5.91 | −2.24 | ok |
| D3 | 64 | 2.25 | [1.45, 5.10] | 0.88 | 4.71 | −0.99 | ok |

- Los perfiles son suaves: el error estándar de la NLL es ≤ 0.17 y no hubo pisos.
- Los resultados son coherentes con Okuyama (2026) y con la validación de la 1.1.0.
- La preferencia por el mecanístico frente al beta-binomial en D2 **depende del ensayo atípico**: sin él, los dos modelos son equivalentes o gana el beta-binomial. Es el mismo diagnóstico de la auditoría del 2026-10-05.
- Costo con el preset final: 37–51 h por conjunto con 2 núcleos, incluida `verify_fit()`. Es un dato para la sección de costo de cómputo.

---

## 3. ¿Responde a Okuyama?

| C | Evidencia de la Fase 5 | Estado |
|---|---|---|
| C2 (distribuciones elegibles) | B: muestra cuándo importa la elección (proceso sí, familia poco) y que el paquete permite el modelo correcto | Respondido |
| C3 y C6 (guía, corridas independientes) | E: la verificación detecta fallas de optimización, pero no la subcobertura; los umbrales deben recalibrarse | **Parcial** hasta corregir H16–H18 |
| C5 (cotas y unidades) | D: 0% de IC abiertos contra 17% en la 1.1.0; ninguna estimación pegada a una cota | Respondido |
| C8 (validación amplia) | A y B: amplio en parámetros y desvíos, acotado en diseños (2) | **Parcial**: la cobertura de 87% hay que corregirla o declararla; faltan diseños (D8) |
| C9 (valores por defecto) | C: exploratory acierta 80% (±0.5) y 94% (±1); standard 100% | Respondido: recomendar standard para informar |
| C7 (z conjunta) | Fase 4 | Respondido |

**Hoy no podemos decir que el intervalo de z tiene cobertura nominal.** Okuyama lo va a mirar primero. Con 87% de cobertura, el texto tendría que decir que el intervalo es anticonservador. Con H16 y H17 corregidos y confirmados, la validación queda sólida.

---

## 4. Propuesta para las decisiones

| ID | Propuesta |
|---|---|
| **D11 (nueva): corregir H16 y H17 antes de congelar** | (a) **Reevaluación adaptativa**: si alguna reevaluación toca el piso, se repite con 4 veces n_sim; si persiste, se marca el punto y no se usa para el IC. El chequeo del piso pasa a mirar todo el perfil. (b) **Refinamiento local**: agregar siempre puntos a mitad de paso (0.25 y luego 0.125) alrededor de z_hat y en los dos tramos donde se cruza el umbral, hasta que el tramo sea ≤ ~1/4 del ancho del intervalo. (c) Interpolación por spline monótona en vez de lineal. Costo estimado: +4 a 8 puntos por ajuste (~+40–60% de tiempo). |
| **D8: corrida de confirmación (Fase 5b)** | Necesaria para validar D11. Repetir el núcleo de A (18 escenarios × 8 conjuntos) con el código corregido, más un factor nuevo de diseño (T en minutos o días, y 13 densidades × 5 réplicas) para cubrir la brecha de C8. Estimación: ~1 día en la PC de 23 núcleos. |
| D9 (mejor de N) | **No implementarlo.** standard es mejor y cuesta parecido. |
| D10 (preset por defecto) | Mantener standard por defecto (ya lo es) y exploratory solo para una primera mirada; está justificado por C. |
| D7 (ABC) | Justificarlo en el texto. Este análisis no cambia esa recomendación. |
| H18 (umbrales) | Recalibrar con la Fase 5b, usando el percentil 90 de los ajustes buenos como umbral de advertencia. |

---

## 5. Correcciones sobre los datos guardados (script 14, sin reajustar)

Se re-evaluó, con 4 veces `n_sim` (hasta 3 rondas) y en los parámetros ya estimados, cada punto del perfil que había tocado el piso de probabilidad (65 puntos, 33 ajustes de 360). Después se recalculó el IC de z.

| z verdadero | Actual | Sin puntos con piso | Re-evaluado, lineal | Re-evaluado, spline |
|---|---|---|---|---|
| 1 | 91.7% | 92.5% | 94.2% | 90.8% |
| 1.75 | 76.7% | 78.3% | 80.0% | 85.0% |
| 2.75 | 91.7% | 92.5% | 91.7% | 95.0% |
| **Global** | **86.7%** | 87.8% | 88.6% | **90.3%** |

- La re-evaluación resuelve el piso en 31 de 33 ajustes (persiste en 2). La mejora de cobertura es la misma que al descartar los puntos, así que la re-evaluación es una alternativa válida y no pierde información.
- El spline suma ~2 puntos, sobre todo a z = 1.75. A z = 1 no ayuda (94 → 91%), dentro del error de 120 ajustes por escenario (IC de Wilson de ±5 puntos).
- **No llega a 95%.** z = 1.75 sigue en 85%. Con lo guardado no se puede mejorar más, porque faltan puntos de perfil cerca de z_hat y de los cruces del umbral. Eso lo corrige D11(b).

**Chequeos (H18).** Con la cobertura corregida y "bueno" = IC cubre, déficit ≤ 1 y |error de z| ≤ 1 (307 buenos, 53 malos), el AUC de los 12 chequeos va de 0.43 a 0.63: ninguno discrimina bien. Es esperable para la cobertura, porque en cada ajuste la probabilidad de que el IC no cubra es ~5% por azar y no deja huella en la optimización. Además, A usó el preset `exploratory`, y chequeos como "Optimiser settled" o "Independent runs agree" dependen del preset. Recalibrar con estos datos mezclaría el preset con la calidad del ajuste. **La recalibración se hace con la Fase 5b usando el preset `standard`,** que es el que se recomienda para informar.
