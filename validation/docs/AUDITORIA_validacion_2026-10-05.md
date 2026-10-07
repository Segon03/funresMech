# Auditoría de la validación de funresMech contra Okuyama (2026)

Proyecto: `funresMech_validacion` · Fecha: 05-oct-2026 · Motor: v1.0 (C++/Rcpp, escala log)

Revisión completa de los análisis, resultados y scripts antes de escribir el
manuscrito y de llevar los cambios al paquete. Complementa a
`docs/MOTOR_RCPP_cambios_y_validacion.md`, que registra los cambios del motor
y la reparametrización.

---

## 1. Veredicto

- **La validación está completa y es correcta para escribir el manuscrito.**
  El motor reproduce exactamente al de Okuyama, los dos ensayos están cerrados
  y todas las corridas son internamente consistentes.
- **Todavía no está validado el paquete funresMech como tal**: lo validado es
  el código de este proyecto (`fun/motor_rcpp.R` + `fun/motor_okuyama.cpp`).
  La versión 1.0.4 del paquete **no** puede reproducir a Okuyama en D2/D3 por
  sus cotas de búsqueda fijas (sección 4). Hace falta portar los cambios y
  verificar la equivalencia antes de afirmar que "funresMech reproduce a
  Okuyama".
- Se encontró una limitación de los datos con consecuencias para los usuarios
  del paquete (un ensayo atípico que hace a z no identificable por arriba en
  D2). Se agregó un diagnóstico previo para detectarlo (sección 5, L1).

---

## 2. Inventario: qué está listo

| Componente | Scripts | Tablas | Figuras | Estado |
|---|---|---|---|---|
| Pruebas del motor | `00b_motor_rcpp_pruebas.R` | `motor_rcpp_pruebas.csv` | — | Listo (8/8 PASA en Windows, R 4.6.1, 28-sep) |
| Ensayo 1: recuperación de parámetros | `02c_ensayo1_rcpp.R`, `03b_ensayo1_rcpp_resumen.R` | `ensayo1_rcpp_resumen.csv`, `_ajustes.csv`, `_perfiles.csv` | `fig_ensayo1_rcpp_perfiles.png`, `_parametros.png` | Listo (1200/1200 puntos) |
| Ensayo 2: reproducción de Okuyama (D2, D3) | `04h`, `04i` | `ensayo2_rcpp_vs_okuyama.csv`, `_resultado.csv`, `_perfiles.csv` | `fig_ensayo2_rcpp_perfiles.png` | Listo (224/224) |
| Ensayo 2: extensión de D2 hasta z = 20 | `04j` (+ `04i`) | `ensayo2_rcpp_extension_*.csv` | `fig_ensayo2_rcpp_extension.png` | Listo (7/7) |
| Ensayo 2: sensibilidad de D2 a la réplica 2/10 | `04k` | `ensayo2_rcpp_D2_sensibilidad*.csv` | `fig_ensayo2_rcpp_D2_sensibilidad.png` | Listo (23/23) |
| Diagnóstico de atípicos (nuevo) | `06_diagnostico_atipicos.R`, `fun/diagnostico_datos.R` | `diagnostico_cribado_D2_D3.csv`, `_ensayo1.csv`, `diagnostico_residuos_pit.csv` | — | Listo (05-oct) |
| Síntesis y reporte | `05_sintesis_y_figuras.R`, `report/validacion_funresMech.Rmd` | — | — | **Desactualizado**: lee las tablas del motor R viejo |

Archivados (no van al manuscrito): Ensayo 1 viejo (`02`, `02b`, `03`,
`results/ensayo1*`, `tablas/ensayo1_*` sin "rcpp") y Ensayo 2 con motor R
(`04`–`04g`, `results/ensayo2`, `results/ensayo2_corr_OKU` salvo `bb_*`).

---

## 3. Verificaciones realizadas

### 3.1 El motor reproduce el modelo de Okuyama

`tablas/motor_rcpp_pruebas.csv` (Windows, R 4.6.1, Rcpp 1.1.2):

- **Idéntico a Okuyama** con la misma semilla: Script S5 (4/4 casos, prueba 1a)
  y Script S4 completo (NLL D3 = 132.04994217 en ambos, prueba 8).
- **Idéntico al motor R de funresMech v1.0.4** (prueba 2: NLL D2 7944.851530 en
  ambos), así que el cambio de motor no cambia el modelo del paquete.
- El corte por saturación no altera la distribución (prueba 3, χ² p > 0.7).
- Hallazgo a documentar: con Rcpp actual, el Script S5 **literal** de Okuyama
  reutiliza números aleatorios cuando s > 0 (`RNGScope` anidado) y da NLL ≈ 940
  en vez de ≈ 132 para D3 (prueba 1c). El motor del proyecto no tiene ese
  problema. Conviene consultarlo con el autor antes de mencionarlo.

### 3.2 Consistencia de las corridas

| Corrida | Puntos | Configuración | NA | Archivos rotos | Parámetro en cota |
|---|---|---|---|---|---|
| Ensayo 2, escala log | 112 | n_sim 10000, 200 gen. | 0 | 0 | 0 |
| Ensayo 2, escala okuyama (ablación) | 112 | igual | 0 | 0 | 67 (`a`, esperado: es la falla de la búsqueda lineal) |
| Extensión D2 | 7 | igual | 0 | 0 | 0 |
| Sensibilidad D2 | 23 | igual | 0 | 0 | 0 |
| Ensayo 1 | 1200 | n_sim 5000, 150 gen. | 0 | 0 | 3 (`k`, cota superior 5) |

### 3.3 Reproducción de la Tabla S1 de Okuyama (escala log, NLL de DEoptim)

| | z_hat (Ok.) | IC95 (Ok.) | AIC (Ok.) | ΔAIC BB − sim (Ok.) | Veredicto revisado |
|---|---|---|---|---|---|
| D2 | 2.8 (3.1) | [1.2, ≥ 20] ([1.1, 4.6]) | 336.572 (336.589) | 10.35 (10.33) | Reproduce, salvo IC superior no acotado |
| D3 | 2.5 (2.5) | [1.5, 4.7] ([1.4, 4.8]) | 272.982 (272.890) | 0.37 (0.62) | Reproduce |

- h coincide a 3 decimales en ambos datasets.
- El AIC beta-binomial de D2 coincide a 3 decimales (346.917 vs 346.916). En D3
  el nuestro es 0.15 menor (273.353 vs 273.507): el ajuste con varios
  arranques (`fit_betabinom`) encuentra un óptimo algo mejor. La conclusión
  (ΔAIC ≈ 0, modelos equivalentes en D3; simulación muy preferida en D2) no
  cambia.
- El z_hat de Okuyama en D2 (3.1) queda a 0.26 unidades de NLL del mínimo,
  dentro del ruido Monte Carlo (tolerancia 0.62).
- El límite superior 4.6 de Okuyama en D2 se explica por su búsqueda lineal de
  `a` en [0, 5]: la ablación en escala okuyama lo reproduce (4.2–4.3), con `a`
  en la cota desde z = 2.5 y la NLL deteriorándose desde z ≈ 4.3. En escala
  log el perfil sigue plano hasta z = 20 (ΔNLL máx. 0.80).

### 3.4 Explicación de D2 (confirmada)

| Variante de D2 (escala log) | n | z_hat | IC95 de z | k en z_hat | ΔNLL en z = 6 / 10 / 20 |
|---|---|---|---|---|---|
| Completo | 73 | 2.8 | [1.12, ≥ 20] | 0.056 | 0.48 / 0.62 / 0.80 |
| Sin la réplica 2/10 (N = 10) | 72 | 2.0 | **[1.40, 4.28]** | **0.60** | 2.78 / 3.60 / 4.40 |
| Control: sin una réplica 10/10 (N = 10) | 72 | — | no acotado | 0.056 | 0.63 / 0.83 / 0.89 |

Quitar ese único ensayo baja la NLL 7.2 unidades, lleva k al rango de D3 y
cierra el intervalo en valores parecidos a los de D3 ([1.4, 4.8]). Quitar otro
ensayo de la misma densidad no cambia nada. El Ensayo 1 (escenario C, k chico)
muestra que el IC superior no acotado es lo esperable con k chico: 80% de las
réplicas (40% en B, 0% en A).

Origen de los datos: D2 y D3 son los datasets 57 y 58 de FoRAGE, ambos de
Kishani Farahani & Goldansaz (2013, Eur. J. Entomol.), *Apanteles myeloenta*
sobre *Ectomyelois ceratoniae*, T = 24 h. **Pendiente: revisar en el trabajo
original si el ensayo 2/10 en N = 10 es un dato real o un error de
digitalización.**

### 3.5 Ensayo 1 (recuperación de parámetros, 20 réplicas por escenario)

| | A: tipo II, z = 1 | B: "D3", z = 2.5 | C: "D2", z = 2.75 |
|---|---|---|---|
| Mediana de z_hat | 1.00 | 2.50 | 2.75 |
| Cobertura del IC95 | 20/20 | 18/20 | 19/20 |
| IC superior no acotado | 0% | 40% | 80% |
| IC inferior > 1 | — | 19/20 | 12/20 |
| h / verdadero (mediana) | 1.02 | 1.00 | 0.99 |
| k / verdadero (mediana) | 1.20 | 4.05 | 0.97 |

### 3.6 Reproducibilidad entre plataformas

Los 60 datasets del Ensayo 1 se regeneraron en Linux (R 4.3.3) con las mismas
semillas y resultaron **idénticos** a los de Windows (R 4.6.1). Las semillas
fijadas en los scripts alcanzan para reproducir los datos simulados.

---

## 4. Adaptación al paquete funresMech (v1.0.4, revisión del código fuente)

Revisado en `github.com/Segon03/funresMech` (rama main, DESCRIPTION 1.0.4).

| # | Problema en v1.0.4 | Dónde | Consecuencia | Solución (validada en este proyecto) |
|---|---|---|---|---|
| P1 | Cotas fijas: a ∈ [0.001, 2], h ∈ [0.001, **0.5**], k ∈ [**0.5**, 5], s (sdlog) ∈ [0.001, **0.5**], z ∈ [0.5, **3**] | `R/fit_full.R`, `R/server.R` (perfil) | **Excluyen los valores reales**: h ≈ 0.61 en D2/D3, k = 0.056 en D2, s de D2 ≈ 0.87 en sdlog; a z alto se necesita a < 10⁻⁵. El paquete no puede reproducir a Okuyama en D2/D3. | Búsqueda en escala log con cotas amplias (`COTAS_LOG`, `theta_a_nat()`), con aviso si un parámetro queda en la cota |
| P2 | Umbral del IC: `min_nll + qchisq(0.95, 1)` (+3.84) | `R/server.R` | IC de ~99.5%, no de 95% | `+ qchisq(0.95, 1) / 2` (+1.92), verificado con el IC beta-binomial de D2 |
| P3 | `plan(multisession)` sin efecto: ni DEoptim ni el perfil usan `future` | `R/server.R` | La opción de núcleos no acelera nada; el perfil es secuencial | `future_lapply` sobre la grilla de z (como 04h) |
| P4 | z_hat y AIC del ajuste con z libre (`fit_full`), ΔAIC contra el perfil | `R/server.R` | Con ruido Monte Carlo el ajuste libre puede quedar peor que el perfil (ΔAIC incoherente) | z_hat = mínimo del perfil (Okuyama, Sec. 2.3) |
| P5 | NLL = `bestval` de DEoptim (mínimo de evaluaciones ruidosas) | `R/server.R` | Sesgo optimista; perfiles "serruchados" | Re-evaluar en el óptimo con semillas nuevas e informar ambas (`ajustar_z_fijo`) |
| P6 | Motor en R puro (~43× más lento) | `R/simulate_trial.R`, `R/simulate_distribution.R` | n_sim y grillas reducidas por costo | `src/motor.cpp` compilado en el paquete; también elimina el problema de la cache de `sourceCpp` (L7) |
| P7 | `simulate_trial(x = 0)` devuelve 0 (hay un test) | `tests/testthat/test-modelos.R` | El motor C++ hace `stop()` con H < 1 | Devolver 0 en la capa R cuando x = 0 |
| P8 | s es el desvío en escala log (sdlog); Okuyama usa el desvío natural | `R/simulate_trial.R` | Los valores de s no son comparables con la Tabla S1 | Adoptar el de Okuyama (`modo = "okuyama"`) o documentar sdlog = √log(1 + (s/h)²) |
| P9 | Fin del ensayo `t >= T` (Okuyama: `t > T`) | `R/simulate_trial.R` | Diferencia de medida cero | Mantener, documentado (el motor tiene los dos modos) |
| P10 | Pedidos de CRAN: sin `set.seed` dentro de funciones; tests de funciones internas | — | — | El motor usa el RNG de R vía `RNGScope`, el usuario controla la semilla; portar las pruebas de `00b` a testthat |

**Verificación de equivalencia (falta):** después de portar, ajustar con el
paquete algunos puntos de z de D2 y D3 (p. ej. z = 2.5 y 4.0) con la misma
semilla y comprobar que la NLL coincide con `nll_motor()` de este proyecto.
Con eso el manuscrito puede afirmar que funresMech v1.1.0 reproduce a Okuyama.

---

## 5. Limitaciones encontradas y soluciones

### L1. Ensayos atípicos e influyentes (enmascaramiento)

- **Evidencia.** Un solo ensayo de D2 (2 de 10 parasitados en N = 10, cuando
  los vecinos parasitan 85–95%) obliga a k ≈ 0.056 y deja a z sin límite
  superior (3.4). Con el modelo ajustado a todos los datos, ese ensayo **no
  parece raro**: P(Y ≤ 2 | N = 10) = 0.014, residuo cuantílico −2.3. Sin él, P =
  5×10⁻⁶ (residuo −4.5). El ajuste "absorbe" al atípico achicando k: los
  residuos posteriores al ajuste no alcanzan para detectarlo.
- **Solución (implementada: `fun/diagnostico_datos.R`, `R/06_diagnostico_atipicos.R`).**
  1. **Antes del ajuste**, `cribado_atipicos()`: Beta-Binomial con la media de
     Okuyama en forma cerrada; para cada ensayo, probabilidad predictiva de su
     resultado con el modelo ajustado *sin* ese ensayo (leave-one-out), con
     umbral de Bonferroni α/n, y cambio relativo en la sobredispersión. Tarda
     ~10 s por dataset.
     - D2: marca **solo** el ensayo 2/10 (p = 5.3×10⁻⁷, umbral 6.8×10⁻⁴; sin
       él, la sobredispersión cae 97%). D3: ninguno.
     - Falsos positivos en los 60 datasets del Ensayo 1: 0/20 en A, 0/20 en B
       y 4/20 en C. En C (k chico, colas pesadas por construcción) hay
       resultados extremos legítimos; el ensayo de D2 es más extremo que todos
       ellos.
  2. **Después del ajuste**, `residuos_pit()`: residuos cuantílicos
     aleatorizados con la predictiva simulada por el motor.
  3. **Análisis de sensibilidad**: re-ajustar sin los ensayos marcados (como
     `04k`) e informar los dos resultados.
- **Recomendación para el paquete.** Nunca descartar datos automáticamente.
  Correr el cribado al cargar los datos, mostrar los ensayos marcados y ofrecer
  el ajuste de sensibilidad. Un ensayo marcado puede ser un error de carga o
  una hembra defectuosa, pero también un resultado real de un proceso con
  mucha variabilidad (escenario C).

### L2. z no identificable por arriba cuando k es chico (cresta z–k)

- **Evidencia.** A z alto el modelo compensa con a ∝ H_ref^(−z) y k más chico
  (z·k ≈ constante); con k chico la probabilidad de encuentro casi no depende
  de a·N^z. Hay IC superior no acotado en 80% de las réplicas de C, 40% de B y
  en D2 (perfil plano hasta z = 20).
- **Solución para el paquete.**
  - Informar el IC como censurado (`[1.1, ≥ z_max]`) en vez de cortarlo en el
    extremo de la grilla.
  - Extender la grilla automáticamente (p. ej. hasta z = 10–20) cuando el
    perfil no cruzó el umbral en el extremo.
  - Mostrar k a lo largo del perfil y avisar si k_hat < ~0.1, que es la señal
    de esta cresta.
  - Seguir reportando la conclusión robusta: el límite inferior (z > 1 = tipo
    III) sí está identificado en la mayoría de los casos.

### L3. k mal identificado cuando es grande

- **Evidencia.** En B (k = 0.82) k se sobreestima ~4 veces, a menudo cerca de
  la cota superior (5); 3 de 1200 puntos quedaron en la cota.
- **Solución.** Subir la cota superior de k en escala log (p. ej. a 50: no
  cuesta nada en escala log), informar k como "≥ cota" cuando la toca y
  advertir en la documentación que k ≳ 1 equivale a búsqueda casi regular.
  No afecta a z ni a h.

### L4. Ruido Monte Carlo y sesgo del mínimo de DEoptim

- **Evidencia.** La NLL re-evaluada supera a la de DEoptim en 0.3–0.8 unidades;
  el ruido de una evaluación con n_sim = 10000 es ~0.2–0.3 unidades. En
  perfiles planos, z_hat depende de ese ruido.
- **Solución (implementada).** Re-evaluación en el óptimo; z_hat comparado por
  ΔNLL frente a una tolerancia de ruido estimada del propio perfil (veredicto
  revisado de `04i`). **Opción para el paquete:** números aleatorios comunes
  (la misma semilla en cada evaluación de la NLL, pasada por el usuario como
  argumento), que suaviza el perfil sin cambiar el modelo.

### L5. Inestabilidad numérica por frecuencias simuladas nulas

- **Evidencia.** Si un conteo observado nunca aparece en las simulaciones, su
  probabilidad se reemplaza por `.Machine$double.xmin` y la NLL salta
  (re-evaluaciones de 705 y 564 en la ablación de D2, error estándar 140).
- **Solución.** Mantener `minp` igual al de Okuyama para la comparación, pero
  avisar cuando ocurre en el óptimo y sugerir más simulaciones. Alternativa
  documentada: un piso de 1/(10·n_sim).

### L6. Costo computacional

- **Evidencia.** Un punto del perfil con la configuración de Okuyama (n_sim
  10000, NP 40, 200 generaciones) tarda 1–2.3 h por núcleo; un perfil
  completo de 56 puntos, ~7 h de pared con 23 núcleos (3 tandas).
- **Solución para el paquete.**
  - Motor C++ y paralelización sobre la grilla de z.
  - Dos configuraciones: "exploración" (n_sim 5000, 150 generaciones, paso
    0.25, validada en el Ensayo 1) y "final" (la de Okuyama).
  - Grilla en dos etapas: gruesa y después refinada (paso 0.1) cerca del
    mínimo y de los cruces del umbral.

### L7. Carrera en la cache de compilación de Rcpp (resuelto)

- **Evidencia.** `sourceCpp()` compartida entre workers dejó vacío su índice y
  cortó dos corridas largas (30-sep y 03-oct).
- **Solución.** Cache por proceso (`fun/motor_rcpp.R`, 04-oct). En el paquete
  el problema desaparece porque el C++ se compila una sola vez, al instalar.

### L8. Búsqueda lineal de los parámetros (hallazgo metodológico central)

- **Evidencia.** Con la búsqueda lineal de Okuyama ([0, 5]⁴), DEoptim casi
  nunca explora `a` < 0.005; desde z ≈ 4.3 el ajuste se deteriora y aparece un
  límite superior artificial. La escala log nunca dio una NLL peor (dentro del
  ruido) y fue mejor en 28 de 112 puntos (hasta 18 unidades).
- **Solución.** La escala log es la recomendada para el paquete. La escala
  lineal se conserva solo como ablación.

### L9. h y s poco informados cuando la tasa de encuentro es baja

- **Evidencia.** En A, h varía entre 0.2 y 1.45 veces el valor verdadero y s se
  subestima (mediana 0.53).
- **Solución.** Es una propiedad del diseño experimental, no del método.
  Mencionarlo en la discusión y en la ayuda del paquete: el tiempo de
  manipulación se estima bien solo si las densidades altas saturan al
  parasitoide.

---

## 6. Pendientes

**Antes de enviar el manuscrito**

1. Portar la sección 8 de `MOTOR_RCPP_cambios_y_validacion.md` y P1–P10 al
   paquete (v1.1.0), con el cribado de atípicos (L1) y los avisos de L2–L5.
2. Verificar la equivalencia del paquete con el proyecto (sección 4).
3. Revisar en Kishani Farahani & Goldansaz (2013) el ensayo 2/10 en N = 10 de
   D2.

**Para el reporte reproducible**

4. Actualizar `R/05_sintesis_y_figuras.R` y `report/validacion_funresMech.Rmd`
   a las tablas `ensayo1_rcpp_*`, `ensayo2_rcpp_*` y `diagnostico_*`.
5. Versiones finales de las figuras, con estilo uniforme y en inglés.

**Opcional**

6. Consultar con Okuyama el problema de `RNGScope` del Script S5 (3.1).
