# Fase 5b: confirmación de funresMech 1.2.0-rc2

Fecha: 2026-10-18 · Script: `validation/R/15_confirmacion_fase5b.R` · Análisis: `validation/R/16_analisis_fase5b.R`
Corrida: 288 trabajos, terminó el 2026-10-18 a las 06:04, sin errores. Código: `v1.2.0-rc2` (D11 con interpolación PCHIP).

| Brazo | Diseño | Preset | Conjuntos |
|---|---|---|---|
| G1 | Núcleo de la Fase 5 (8 y 5 densidades) | exploratory | 144 (los mismos datos que la Fase 5) |
| G2 | Núcleo, conjuntos 1–4 | standard | 72 (los mismos datos que G1) |
| G3 | 13 densidades × 5 réplicas | exploratory | 72 |

Todos con `verify_fit()`. z verdadero ∈ {1, 1.75, 2.75}, k ∈ {0.2, 1, 5}.

---

## 1. Verificación

| Chequeo | Resultado |
|---|---|
| Trabajos completos | 288/288 |
| Cada IC contiene su z_hat | Sí (el error de `monoH.FC` quedó corregido) |
| G1 usa los mismos datos que la Fase 5 | 144/144 idénticos |
| G2 usa los mismos datos que G1 | 72/72 idénticos |
| Fallas de optimización (déficit > 1) | 0 de 288 |
| Ajustes con piso de probabilidad persistente | 1 de 288 (40 re-evaluaciones lo resolvieron en los demás) |

**Pasa.**

---

## 2. Resultados

### 2.1 Cobertura del IC de z

| Brazo | z = 1 | z = 1.75 | z = 2.75 | **Global** [IC 95%] |
|---|---|---|---|---|
| G1 (núcleo, exploratory) | 89.6% | 85.4% | 95.8% | **90.3%** [84.3, 94.1] |
| G2 (núcleo, standard) | 91.7% | 79.2% | 95.8% | **88.9%** [79.6, 94.3] |
| G3 (13 densidades) | 95.8% | 91.7% | 95.8% | **94.4%** [86.6, 97.8] |
| **Todos** | | | | **91.0%** [87.1, 93.8] (n = 288) |

- La cobertura global queda **por debajo del 95%** (prueba binomial, p = 0.004). El intervalo es levemente anticonservador.
- Los fallos son **simétricos**: el z verdadero queda por debajo del IC tantas veces como por encima. No hay sesgo de dirección.
- Por factor (G1 + G3, exploratory): k = 0.2 → 87.5%; k = 1 → 97.2%; k = 5 → 90.3%. Por diseño: 8 densidades 87.5%, 5 densidades 93.1%, 13 densidades 94.4% (diferencias no significativas con 72 ajustes por nivel; Fisher p = 0.13).

### 2.2 Comparaciones pareadas (mismos datos)

| Comparación | Cobertura | Discordantes | Ancho mediano del IC | Tiempo por ajuste |
|---|---|---|---|---|
| Fase 5 (1.2.0 sin D11) → G1 (rc2) | 90.3% → 90.3% | 5 y 5 (McNemar p = 1) | 1.04 → 0.98 | 0.80 → 0.93 h (+16%) |
| G1 exploratory → G2 standard | 88.9% → 88.9% | 0 y 0 | 0.99 → 1.00 | 0.93 → 4.18 h (×4.5) |

- **D11 no cambia la cobertura** en el núcleo. Corrige casos patológicos (pisos de probabilidad, perfiles con pocos puntos cerca del límite) y lo hace a bajo costo.
- **El preset standard no cambia la cobertura**: los mismos 72 conjuntos se cubren o fallan igual. Su NLL es mejor por más de 0.5 en solo el 8% de los ajustes.

**Conclusión principal:** la subcobertura (~91%) **no es un problema numérico** (optimización, ruido Monte Carlo, grilla o interpolación). Esas cuatro causas quedaron descartadas, una por una, entre la Fase 5 y la 5b. Es una propiedad estadística del intervalo de razón de verosimilitud perfilada con muestras finitas: el umbral χ²₁/2 = 1.92 es asintótico, y con 30–65 ensayos y parámetros débilmente identificados (k chico, colinealidad z–k) el perfil es más estrecho de lo que debería. Mejora con más densidades (G3: 94.4%).

### 2.3 Error tipo I, potencia y sesgo

| | z = 1 (error tipo I: IC excluye 1) | z = 1.75 (potencia) | z = 2.75 (potencia) |
|---|---|---|---|
| G1 | 8.3% | 87.5% | 100% |
| G2 | 8.3% | 83.3% | 100% |
| G3 | 4.2% | 100% | 100% |

- El error tipo I es coherente con la subcobertura (nominal 5%, observado 4–8%). El criterio ΔAIC > 2 frente a z = 1 da cifras parecidas (4–10%).
- El sesgo de z es nulo (mediana del error entre −0.03 y +0.25). h se recupera sin sesgo (h_hat/h mediano entre 0.99 y 1.12).
- **k se sobreestima** (k_hat/k mediano entre 1.1 y 1.6): es un parámetro de forma estimado con pocas réplicas por densidad. Se informa como limitación conocida y no afecta a z ni a h.

### 2.4 Chequeos de confiabilidad (H18): resuelto sin cambiar umbrales

Con el **preset standard** (G2, 72 ajustes) ningún chequeo marca `fail`, y el porcentaje de advertencias en ajustes buenos es bajo:

| Chequeo | Advertencias en ajustes buenos (standard) | En exploratory |
|---|---|---|
| Independent runs agree (NLL) | 6.3% | 50% |
| Optimiser settled before the end of the run | 0% | 58% |
| Monte Carlo noise of the NLL | 3.1% | 30% |
| Independent runs agree (a, h) | 3.1% | 7% |
| Interval of z closed | 1.6% | 0.5% |
| Los demás | 0% | 0–4% |

- **Los umbrales actuales están bien calibrados para el preset que se recomienda informar.** Las advertencias masivas de la Fase 5 venían del preset exploratory, y son correctas: ese presupuesto no alcanza para informar resultados, y los chequeos lo detectan.
- Ningún chequeo predice si el IC cubre el z verdadero (AUC 0.33–0.68). Es lo esperable: el ~5–9% de no cobertura es azar muestral y no deja huella en la optimización. Los chequeos miden la **confiabilidad numérica** del ajuste, no la cobertura, y así se los describe en el manuscrito.
- **Decisión H18: mantener los umbrales.** Se documenta su tasa de falsas alarmas con standard (≤ 6%).

---

## 3. Estado de la respuesta a Okuyama

| C | Evidencia | Estado |
|---|---|---|
| C3 y C6 (guía de fiabilidad, corridas independientes) | Con standard, ≤ 6% de falsas alarmas y 0 fallas; con exploratory los chequeos avisan que hay que subir el presupuesto | **Respondido** |
| C8 (validación amplia) | 5 brazos en la Fase 5 + 5b, 3 diseños (5, 8 y 13 densidades), 9 combinaciones z × k, desvíos del modelo, presupuesto, 1.1.0 vs 1.2.0. Cobertura 91% (87–94%), causas numéricas descartadas | **Respondido**, declarando que el IC es levemente anticonservador |
| C9 (valores por defecto) | standard no mejora la cobertura frente a exploratory, pero sí la confiabilidad numérica (chequeos limpios) | Respondido: recomendar standard o final para informar |

---

## 4. Propuestas

| ID | Propuesta |
|---|---|
| **D12. Cómo informar la cobertura** | Declararla en el manuscrito: "the profile-likelihood interval for z had 91% coverage (95% CI 87–94%) across 288 simulated data sets, lower with highly variable search times (k = 0.2) and higher with 13 densities (94%)". Recomendar diseños con ≥ 8 densidades y, si se necesita cobertura exacta, un IC calibrado por bootstrap paramétrico (no implementado, por costo). |
| D13. Bootstrap paramétrico del IC (opcional) | No implementarlo en la 1.2.0: cada IC costaría unas 200 veces un ajuste. Mencionarlo como trabajo futuro. |
| D14. Congelar `v1.2.0-rc2` como código final | No queda ningún cambio de método pendiente. Lo que sigue es documentación (NEWS, README, viñetas, `R CMD check`), manuscrito y respuesta a Okuyama. |

Figura: `validation/figuras/fase5b_cobertura.png`.
