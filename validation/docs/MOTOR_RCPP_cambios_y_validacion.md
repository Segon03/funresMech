# Motor Rcpp + reparametrización log: cambios, hallazgos y plan para funresMech

Proyecto: `funresMech_validacion` · Fecha: 28-sep-2026 · Motor: v1.0

Este documento registra (1) qué se cambió en el proyecto de validación y por qué,
(2) los hallazgos técnicos que lo motivaron, con la evidencia, (3) la lista de
cambios a llevar al paquete **funresMech** una vez cerrada la validación, y (4)
texto base para el manuscrito. El paquete funresMech **no** se modifica hasta que
la validación esté completa.

---

## 1. Resumen

La validación contra Okuyama (2026) con el motor R puro de funresMech v1.0.4
(scripts 04d/04f/04g) no reproducía la Tabla S1: el perfil de verosimilitud de z
quedaba "serruchado" e inflado a z altos (p. ej. D2, z = 4.5: NLL 169.0 cuando lo
alcanzable es ~165), y el límite superior del IC de D2 quedaba censurado. Se
identificaron tres causas, ninguna de ellas del modelo en sí:

1. **Velocidad**: el motor R puro es ~65–80× más lento que el C++ de Okuyama. Eso
   obligó a usar n_sim = 3000 (Okuyama: 10000) y cortes de tiempo por evaluación.
2. **Trampa computacional**: con h ≈ 0 y k ≈ 0 (o a·H^z enorme) el bucle
   `while (t < T)` da decenas de millones de vueltas por ensayo sin cambiar el
   resultado. En R eso "cuelga" un worker; en 04d/04f se esquivó poniendo un piso
   artificial h ≥ 0.01, que recorta el espacio de parámetros.
3. **Búsqueda en escala lineal**: DEoptim muta sumando diferencias en la escala de
   las cotas. Con cotas lineales [0, 5], valores como a ≈ 10⁻⁵ (necesarios a z alto)
   o k ≈ 0.01 ocupan una fracción ínfima del rango y casi no se exploran.

Solución adoptada:

| Cambio | Qué resuelve | ¿Altera el modelo? |
|---|---|---|
| Motor en C++ vía Rcpp (`fun/motor_okuyama.cpp`) | velocidad (causa 1) | No: reproduce bit a bit al motor de Okuyama y al de funresMech |
| Corte por saturación | trampa (causa 2), sin pisos artificiales | No: exacto (prueba 3) |
| Búsqueda de DEoptim en escala log | exploración (causa 3) | No: misma verosimilitud, otra parametrización del espacio de búsqueda |
| n_sim = 10000, grilla de z con paso 0.1, s como desvío natural | fidelidad a Okuyama | No |
| Re-evaluación de la NLL en el óptimo | sesgo optimista del mínimo de DEoptim | No (se informa aparte) |

"Seguirá siendo un paquete de R": sí. Rcpp compila una función en C++ que R llama
como a cualquier otra; el usuario no ve ninguna diferencia salvo la velocidad.

---

## 2. Archivos del proyecto de validación

Nuevos (no se borró ni se modificó ningún archivo anterior; 04d–04g quedan como
registro del motor R puro):

| Archivo | Contenido |
|---|---|
| `fun/motor_okuyama.cpp` | (A) copia literal del Script S5 de Okuyama; (A') S5 con la corrección del RNGScope (§4); (B) el motor del proyecto `motor_okuyama_cpp()` con `modo` okuyama/funresMech y `corte_saturacion`; utilidades de diagnóstico |
| `fun/motor_rcpp.R` | `cargar_motor_rcpp()` (compila una vez, con caché en `fun/.cache_rcpp/`), `nll_gamma_lnorm_S4()` (copia literal del Script S4), `nll_motor()`, `negloglik_fixed_z_rcpp()` (reemplazo directo de la función del paquete), reparametrización log (`theta_a_nat()`, `nat_a_theta()`, `COTAS_LOG`) y `ajustar_z_fijo()` |
| `R/00b_motor_rcpp_pruebas.R` | 10 pruebas del motor (§5). Escribe `tablas/motor_rcpp_pruebas.csv` |
| `R/04h_ensayo2_rcpp_log.R` | Perfiles de z de D2 y D3 (0.5–6.0, paso 0.1) con el motor Rcpp; escala "log" (principal) y "okuyama" (ablación). Checkpoint por punto en `results/ensayo2_rcpp/<escala>/` |
| `R/04i_ensayo2_rcpp_reensamblar.R` | z_hat, IC95 (grilla e interpolado), parámetros, AIC, ΔAIC vs beta-binomial, veredicto vs Tabla S1, figura |
| `docs/MOTOR_RCPP_cambios_y_validacion.md` | este documento |

Salidas de 04i: `tablas/ensayo2_rcpp_perfiles.csv`, `tablas/ensayo2_rcpp_resultado.csv`,
`tablas/ensayo2_rcpp_vs_okuyama.csv`, `figuras/fig_ensayo2_rcpp_perfiles.png`.

### Cómo correr

Requisitos (Windows): **Rtools** de la misma versión que R
(<https://cran.r-project.org/bin/windows/Rtools/>). Verificar con
`Sys.which("make")` (tiene que devolver una ruta) y `install.packages("Rcpp")`.

1. `R/00b_motor_rcpp_pruebas.R` — compila el motor (~1 min la primera vez) y corre las
   pruebas (~2–5 min). Todas deben dar **PASA**.
2. (Opcional, recomendado) `R/04h_ensayo2_rcpp_log.R` con `MODO_PRUEBA <- TRUE`:
   corrida de humo de ~2 min. Luego `04i` con `MODO_PRUEBA <- TRUE`. Solo verifica
   que todo funcione; los números no significan nada.
3. `R/04h_ensayo2_rcpp_log.R` con `MODO_PRUEBA <- FALSE`. Imprime al inicio el
   tiempo estimado. Referencia: ~0.4 s por evaluación de la NLL → ~0.9 h-núcleo por
   punto del perfil; 224 puntos (2 escalas × 2 datasets × 56 z) ≈ 200 h-núcleo ≈
   **8–9 h con 23 workers** (la mitad si se deja solo `ESCALAS <- "log"`). Las tareas
   "log" se despachan primero.
4. `R/04i_ensayo2_rcpp_reensamblar.R` — segundos; se puede correr con 04h a medio
   terminar para mirar resultados parciales.

---

## 3. Diagnóstico del motor R puro (por qué fallaba la estimación)

- **Mismo modelo, misma verosimilitud.** funresMech v1.0.4 implementa el mismo
  proceso que Okuyama: tiempo de búsqueda Gamma(forma k, tasa k·a·H^z), huésped
  elegido al azar, manipulación lognormal de media h, resultado = huéspedes
  distintos parasitados; NLL = −Σ log p̂(H1 | H0) con p̂ de n_sim ensayos y piso
  `.Machine$double.xmin`. Verificado bit a bit (prueba 2).
- **Pero no llegaba al óptimo** a z alto: a z = 4.5 el `a` óptimo es ~10⁻⁵–10⁻⁴,
  que en la escala lineal [10⁻⁶, 5] de 04d/04f casi no se explora. Re-escalando
  solo `a` a mano desde los parámetros de z = 3.1 se obtenían NLL 1.6–4 unidades
  mejores que las de DEoptim en z = 3.7–4.9 (chequeo del 25-sep).
- **La trampa computacional** (h ≈ 0, k ≈ 0): 10 ensayos a H = 50 con h = k = 10⁻⁶
  requieren 34.5 millones de vueltas del bucle; con el corte por saturación, 2238
  (prueba 4). En el simulador de Okuyama la trampa existe igual, pero el C++ la
  hace tolerable dentro de su rango de z; al llevarlo a z = 5.9 también se "colgó"
  (>14 min sin terminar una corrida, prueba hecha el 26-sep).

---

## 4. Hallazgo: el Script S5 literal reutiliza números aleatorios con Rcpp actual

`rlnorm0_cpp()` abre su propio `RNGScope` y se llama **dentro** del bucle de
`parasitism_gamma_lnorm_cpp()`, que ya tiene uno abierto. Al menos desde Rcpp 1.0.0 (2018),
`RNGScope` no tiene contador: cada apertura llama a `GetRNGstate()`, que relee
`.Random.seed` — que el bucle externo todavía no actualizó. Resultado: cada tiempo
de manipulación reutiliza números aleatorios ya usados para la búsqueda y la
elección del huésped. Se verificó en el código fuente de Rcpp (`src/api.cpp`,
versiones 1.0.0 a 1.1.0: `enterRNGScope()` llama siempre a `GetRNGstate()`) y
empíricamente (prueba 1c):

- Demostración mínima: un `RNGScope` anidado devuelve el mismo número que la
  llamada anterior.
- Distribución (H = 2, parámetros de D3 de Okuyama), P(0/1/2 parasitados):
  R puro independiente 0.201/0.491/0.308; motor del proyecto 0.198/0.494/0.308
  (χ² p = 0.46); **S5 literal 0.200/0.635/0.165** (χ² p ≈ 0).
- NLL de D3 en los parámetros publicados (z = 2.5; la implícita en el AIC de
  Okuyama es 131.44): S5 literal ≈ 940; S5 corregido ≈ 132.

Interpretación prudente: los resultados publicados de Okuyama son consistentes con
un simulador **correcto** (su NLL coincide con la de nuestro motor en sus
parámetros, prueba 6), así que sus estimaciones no parecen afectadas; lo afectado
es el código suplementario tal como se distribuye, al compilarlo con Rcpp actual
(solo cuando s > 0). Con s = 0 el S5 literal y nuestro motor son idénticos
(prueba 1b). La corrección es de una línea (no abrir un `RNGScope` dentro del
bucle). **Sugerencia**: comunicarlo al autor antes de mencionarlo en el paper.

---

## 5. Pruebas del motor (`R/00b_motor_rcpp_pruebas.R`)

Resultado en Linux, R 4.3.3, Rcpp 1.0.12 (28-sep-2026): **10/10 PASA**. Repetir en
la PC de cómputo (Windows) antes de 04h.

| # | Prueba | Resultado de referencia |
|---|---|---|
| 1a | Motor = S5 (con RNGScope corregido), misma semilla | 4/4 casos idénticos |
| 1b | Motor = S5 literal cuando s = 0 | idéntico |
| 1c | Diagnóstico RNGScope (§4) | confirmado |
| 2 | Modo funresMech = motor R de funresMech v1.0.4, misma semilla | `simulate_trial` 4/4 idénticos; `negloglik_fixed_z` idéntica |
| 3 | Corte por saturación no cambia la distribución | χ² p = 0.72, 0.89, 0.73; NLL D2 164.10 vs 164.16 (p = 0.68) |
| 4 | Trampa h≈0/k≈0 | 34.5 M vueltas → 2238; NLL completa en 0.8 s |
| 5 | Velocidad | 0.41 s por NLL (D2, n_sim = 10000) vs ~27 s en R puro (~65×) |
| 6 | Anclaje con Tabla S1 | D2: 164.3 vs 163.3 (a publicado con 1 cifra); D3: 132.0–132.2 vs 131.4 |
| 7 | Reparametrización ida y vuelta | error relativo 1e-15 |
| 8 | NLL del Script S4 = `nll_motor()` | idéntica a 8 decimales |

### 5b. Prueba de la reparametrización con la configuración completa (D2, z = 4.5)

Un punto del perfil donde el motor R puro fallaba, con la configuración definitiva
de 04h (n_sim = 10000, NP = 40, itermax = 200, 5 re-evaluaciones), mismo motor
Rcpp, cambiando solo la escala de búsqueda (corrido en la nube, 1 núcleo cada uno):

| Búsqueda | NLL DEoptim | NLL re-evaluada | a | h | k | s | Tiempo |
|---|---|---|---|---|---|---|---|
| Motor R v1.0.4, lineal (04f, n_sim = 3000) | 169.01 | — | — | — | — | — | — |
| Rcpp, escala lineal de Okuyama | 166.25 | 166.31 ± 0.15 | 6.7e-5 (pegado al piso de la escala lineal) | 0.665 | 0.023 | 0.831 | 63 min |
| **Rcpp, escala log** | **163.62** | **164.49 ± 0.10** | 3.9e-5 | 0.626 | 0.026 | 0.669 | 57 min |

La escala log mejora la NLL en ~1.8 unidades respecto de la búsqueda lineal con el
mismo motor, y en ~5 respecto del motor R anterior. Referencia: el mínimo de
Okuyama para D2 es 163.29 y su IC llega a z = 4.6, o sea que en su perfil z = 4.5
está dentro del IC (NLL < 163.29 + 1.92 = 165.2); con la escala lineal quedaría
afuera (comparación aproximada: el perfil completo de 04h es el que decide). Además h, k y s en z = 4.5 con la escala log quedan muy cerca de los
publicados en z_hat (0.606 / 0.037 / 0.682), como se espera de un perfil suave.

---

## 6. Diferencias funresMech v1.0.4 vs Okuyama (2026) a tener en cuenta

| Aspecto | funresMech v1.0.4 | Okuyama (2026) | Decisión en la validación |
|---|---|---|---|
| Motor | R puro | C++ (Rcpp) | C++ |
| n_sim | configurable (app) | 10000 | 10000 |
| `s` | desvío en escala **log** (sdlog) | desvío en escala **natural**, buscado como s = u² | como Okuyama (comparable con Tabla S1). Equivalencia: sdlog = √log(1 + (s/h)²) |
| Cotas de búsqueda | fijas en el código: a∈[0.001,2], h∈[0.001,0.5], z∈[0.5,3], k∈[0.5,5], s∈[0.001,0.5] | [0,5] para a, h, k, √s | escala log amplia (§7) |
| Estimador de z | ajuste con z libre (5 parámetros) | mínimo del perfil sobre grilla | mínimo del perfil |
| Umbral del IC95 | NLL_min + qchisq(.95,1) = +3.84 (**error**) | +qchisq(.95,1)/2 = +1.92 | +1.92 (verificado a 3 decimales con el IC beta-binomial de D2) |
| T | la app lo pide; los scripts viejos usaban 1 | 24 h | 24 |
| Fin del ensayo | `t >= T` | `t > T` | como Okuyama (diferencia de medida cero) |

---

## 7. Reparametrización log (detalle)

DEoptim busca en θ = (log λ_ref, log h, log k, u) con

- λ_ref = a·H_ref^z (tasa de encuentro a la densidad de referencia; H_ref = media
  geométrica de las densidades distintas: D2 ≈ 15.4). Con z fijo es una
  reparametrización 1 a 1 de `a`, pero hace que el rango útil de búsqueda **no
  dependa de z** (a z alto `a` debe achicarse ~H_ref^z veces y λ_ref lo absorbe).
- cotas: λ_ref ∈ [10⁻³, 10³], h ∈ [10⁻³, 5], k ∈ [10⁻⁵, 5], u ∈ [0, 5] (s = u²
  ∈ [0, 25], igual que Okuyama).
- u queda lineal porque s = 0 (manipulación fija) debe ser alcanzable.

La escala "okuyama" (lineal, [0, 5]⁴, s = p₄²) se conserva en 04h como ablación.

---

## 8. Cambios propuestos para funresMech (próxima versión)

A aplicar **solo después** de cerrar la validación (04h/04i con veredicto
REPRODUCE). Orden sugerido:

1. **Motor en C++**
   - Agregar `src/motor.cpp` con `motor_okuyama_cpp()` (sin las funciones de
     diagnóstico ni la copia literal de S5).
   - DESCRIPTION: `LinkingTo: Rcpp`, `Imports: Rcpp`; NAMESPACE:
     `useDynLib(funresMech, .registration = TRUE)` e `importFrom(Rcpp, sourceCpp)`;
     correr `Rcpp::compileAttributes()`.
   - Reescribir `simulate_distribution()` y `negloglik_fixed_z()` sobre el motor
     (modelo: `negloglik_fixed_z_rcpp()` en `fun/motor_rcpp.R`), manteniendo sus
     firmas para no romper `server.R`.
   - Mantener `simulate_trial()` en R solo como referencia interna para los tests.
2. **Corte por saturación** activado por defecto.
3. **Búsqueda en escala log** en `fit_full()` y en el perfil de `server.R`
   (λ_ref, log h, log k, u). Cotas amplias y, opcionalmente, configurables.
4. **Umbral del IC**: `min_nll + qchisq(0.95, 1)/2` en `server.R` (hoy sin `/2`).
5. **z_hat** = mínimo del perfil (o informar ambos: perfil y z libre, advirtiendo si
   el libre es peor que el perfil).
6. **Parametrización de s**: decidir entre (a) adoptar la de Okuyama (desvío
   natural) o (b) mantener sdlog y documentar la conversión en la ayuda y en los
   reportes. Recomendación: (a), para que los usuarios comparen directo con Okuyama.
7. **Tests (testthat)** — también responden al pedido de la revisora de CRAN de
   testear funciones no exportadas:
   - motor C++ idéntico a `simulate_trial()` en R con la misma semilla (prueba 2);
   - corte por saturación: misma distribución (prueba 3, con n chico para CRAN);
   - ida y vuelta de la reparametrización (prueba 7);
   - umbral del IC (caso beta-binomial D2 → [1.310, 2.552]).
8. **CRAN**: las funciones no deben fijar semilla internamente (pedido de la
   revisora); el motor usa el RNG de R vía `RNGScope`, así que `set.seed()` del
   usuario lo controla. Revisar que el código C++ no abra `RNGScope` anidados
   (§4). `R CMD check --as-cran` en Windows (win-builder) por la compilación.
9. **NEWS.md**: registrar todo lo anterior como v1.1.0.

---

## 9. Texto base para el manuscrito (borrador, en inglés)

**Methods — simulation engine and optimisation.** The simulation likelihood of
the mechanistic model was computed with a C++ implementation of the stochastic
foraging process (Okuyama 2026), called from R through Rcpp (Eddelbuettel &
François 2011). With identical random seeds, the engine reproduces exactly the
output of the reference R implementation in funresMech and of Okuyama's (2026)
supplementary simulator. Each trial is stopped once all hosts have been
parasitised; because the response variable is the number of distinct parasitised
hosts, this rule does not change the simulated distribution (verified with χ²
tests on 2×10⁵ simulated trials), but it removes pathological run times when
handling time and the Gamma shape parameter approach zero. For each value of z on
a grid from 0.5 to 6.0 (step 0.1), the remaining parameters were estimated with
differential evolution (DEoptim; NP = 40, 200 generations, 10,000 simulated
trials per density level), searching on log(a·H_ref^z), log h, log k and √s, where
H_ref is the geometric mean of the host densities. The log-scale search keeps the
search range independent of z and improves exploration of small values of a and
k. The point estimate of z was the minimum of the profile, and its 95% confidence
interval comprised the z values with NLL ≤ NLL_min + χ²₁,₀.₉₅/2.

**Results** — completar con `tablas/ensayo2_rcpp_vs_okuyama.csv` y la figura
`figuras/fig_ensayo2_rcpp_perfiles.png` (z_hat, IC, AIC, ΔAIC vs beta-binomial de
D2 y D3 frente a Okuyama 2026, y la comparación escala log vs escala lineal).

**Nota posible (discusión o material suplementario)**: el problema del `RNGScope`
anidado del Script S5 (§4), solo después de consultarlo con el autor.

---

## 10. Estado

- [x] Motor C++ y capa R escritos y probados (00b: 10/10 PASA en Linux).
- [x] 04h y 04i probados de punta a punta en modo prueba.
- [x] Escala log vs lineal en D2, z = 4.5, configuración completa (§5b).
- [x] 00b en la PC de cómputo (Windows).
- [x] 04h completo (224 perfiles) → 04i → veredicto (30-sep / 01-oct-2026):
      D3 REPRODUCE en las dos escalas. D2: z_hat y AIC reproducen; el límite
      superior del IC no está identificado (perfil plano). 04i agrega un
      veredicto revisado (ver encabezado de 04i).
- [x] 04j: extensión de D2 en escala log hasta z = 20 (02-oct-2026): perfil
      plano (ΔNLL máx. 0.80), ningún parámetro en cota → IC de z en D2 =
      [1.1, ≥20]. El 4.6 de Okuyama se explica por su búsqueda lineal de a.
- [ ] Ensayo 1 definitivo (R/02c → R/03b): 3 escenarios × 20 réplicas,
      ~50 h. Piloto en la nube (02-oct): B y C se comportan como se esperaba;
      A se corrigió de a = 0.65 a a = 0.05 (con 0.65 saturaba las densidades
      bajas y no rechazaba z = 2).
      04-oct: la corrida se cortó a las ~35 h por la cache de Rcpp compartida
      entre workers ("archivo de entrada vacío"; el mismo error del 30-sep con
      04h). Corregido en fun/motor_rcpp.R (cache por proceso, en tempdir()) y
      02c ahora tolera fallas de tareas sueltas y escribe los .rds de forma
      atómica. Se retoma desde donde quedó.
- [x] Ensayo 1 completo (1200/1200, 04-oct): medianas de z_hat = valor verdadero en A, B y C;
      cobertura 20/20, 18/20, 19/20.
- [x] 04k (04-oct): sin la réplica 2/10, k sube de 0.056 a 0.60 y el IC de z se cierra en
      [1.40, 4.28]; el control (sin una réplica 10/10) sigue plano. Hipótesis confirmada.
- [x] Auditoría completa (05-oct): ver docs/AUDITORIA_validacion_2026-10-05.md. Nuevo
      diagnóstico de atípicos: fun/diagnostico_datos.R + R/06_diagnostico_atipicos.R.
- [ ] Actualizar R/05_sintesis_y_figuras.R y report/validacion_funresMech.Rmd
      a las tablas nuevas (ensayo1_rcpp_*, ensayo2_rcpp_*).
- [ ] Llevar §8 a funresMech (solo con la validación cerrada) y verificar que
      el paquete reproduce lo validado (misma semilla, algunos z).
