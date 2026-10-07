# funresMech_validacion

Proyecto de RStudio para la sección **"Validación del paquete funresMech con
el modelo de Okuyama como motor"** del manuscrito *Núñez-Campero (2026) —
funresMech: An R Package for Mechanistic Functional Response Analysis using
the Okuyama Model*.

Reemplaza y unifica lo que originalmente eran las Secciones 2.4 (validación
con datos simulados) y 2.5 (comparación con datos reales de FoRAGE) en dos
ensayos livianos: **Ensayo 1** (recuperación de parámetros conocidos, 60
ajustes) y **Ensayo 2** (reproducción de las conclusiones de Okuyama 2026 en
2 datasets reales, D2 y D3, 4 ajustes). Ver el documento
`Validacion_funresMech_Okuyama.docx` (carpeta `Analisis_FoRAGE`, un nivel
arriba de este proyecto) para el detalle completo del diseño.

## ACTUALIZACIÓN 28-sep-2026 — Motor Rcpp (corrida de validación definitiva)

El Ensayo 2 se rehace con un motor de simulación en C++ (Rcpp), idéntico al de
Okuyama (2026), con corte exacto por saturación y búsqueda de DEoptim en escala
log. Detalle, hallazgos y plan para el paquete en
**`docs/MOTOR_RCPP_cambios_y_validacion.md`**. Los scripts 04d–04g (motor R puro)
se conservan como registro, no hace falta volver a correrlos.

Requisito extra en Windows: **Rtools** de la misma versión que R
(`Sys.which("make")` debe devolver una ruta) y `install.packages("Rcpp")`.

Orden:

1. `R/00b_motor_rcpp_pruebas.R` — compila el motor y corre 10 pruebas (todas PASA).
2. `R/04h_ensayo2_rcpp_log.R` — perfiles de z de D2 y D3 (paso 0.1, n_sim = 10000).
   Primero conviene una corrida de humo con `MODO_PRUEBA <- TRUE` (~2 min).
3. `R/04i_ensayo2_rcpp_reensamblar.R` — tabla vs Okuyama, veredicto y figura.

## ACTUALIZACIÓN 02-oct-2026 — Ensayo 1 definitivo (motor Rcpp)

El Ensayo 1 se rehace con el motor Rcpp y la búsqueda en escala log, con datos
simulados por el **propio modelo de Okuyama** sobre el diseño experimental de D3
(13 densidades 2–50, 64 ensayos, T = 24 h). Tres escenarios (A: tipo II, z = 1,
a = 0.05; B: tipo III "como D3", z = 2.5; C: tipo III con k chico "como D2", z = 2.75),
20 réplicas cada uno, perfil de z de 0.25 a 5.0 (paso 0.25), n_sim = 5000,
150 generaciones. Detalle en el encabezado de `R/02c_ensayo1_rcpp.R`.

1. `R/02c_ensayo1_rcpp.R` — simula los 60 datasets y ajusta los perfiles
   (1200 tareas, ~50 h con 23 workers). Primero `MODO_PRUEBA <- TRUE` (~1 min).
2. `R/03b_ensayo1_rcpp_resumen.R` — sesgo, RMSE, cobertura del IC95 de z y
   figuras (`tablas/ensayo1_rcpp_*.csv`, `figuras/fig_ensayo1_rcpp_*.png`).

**Ensayo 1 viejo ARCHIVADO** (no va al manuscrito): `R/02_ensayo1_simulacion.R`,
`R/02b_ensayo1_confirmacion.R`, `R/03_ensayo1_resumen.R`, `results/ensayo1/`,
`results/ensayo1_confirmacion/` y `tablas/ensayo1_*.csv` (sin "rcpp") se
conservan solo como registro: motor R puro, n_sim = 300, búsqueda lineal y
datos Beta-Binomial (densidades 5–100, T = 1).

## ACTUALIZACIÓN 05-oct-2026 — Auditoría y diagnóstico de atípicos

Estado completo, verificaciones y lista de cambios para el paquete en
**`docs/AUDITORIA_validacion_2026-10-05.md`**. Nuevo:
`R/06_diagnostico_atipicos.R` (usa `fun/diagnostico_datos.R`): cribado
leave-one-out de ensayos atípicos antes del ajuste y residuos cuantílicos
después (~10–15 min).

## Instalación

```r
install.packages(c("remotes", "DEoptim", "extraDistr", "readr", "dplyr",
                    "ggplot2", "rmarkdown", "knitr", "future", "future.apply"))
remotes::install_github("Segon03/funresMech", ref = "v1.0.4")

# Generar renv.lock la primera vez (el que viene en el repo es un
# placeholder minimo, ver mas abajo):
# install.packages("renv")
# renv::init()
# renv::snapshot()
```

## Orden de ejecución

1. **`R/00_benchmark.R`** — CORRER PRIMERO. Mide cuánto tarda un ajuste real
   antes de comprometerse a los barridos completos. Incluye un chequeo
   casi instantáneo con los mismos valores "chicos" que usa la propia
   suite de tests del paquete (`itermax=20, NP=20, n_sim=50`).
2. **`R/01_extract_forage_d2_d3.R`** — extrae D2 y D3 del CSV de FoRAGE ya
   copiado en `data_raw/` y los deja en `data_clean/D2.rds`, `data_clean/D3.rds`
   con las columnas `dens`/`par` que espera `funresMech:::fit_full()`.
3. **`R/02_ensayo1_simulacion.R`** — barrido de 60 ajustes simulados
   (ρ = 0, 0.05, 0.1 × 20 réplicas), **en paralelo** (`future.apply`,
   usa todos los cores menos uno). Configuración real: `n_sim=300`,
   `NP=20` (reducidos desde 1000/40 — ver "Cómputo e historial de
   tiempos" más abajo). El IC95% de z por réplica (perfil local) sólo
   se calcula para 3 réplicas de control por nivel de ρ, no las 60 —
   ver el comentario al principio del script. Guarda cada ajuste en
   `results/ensayo1/` apenas termina (checkpointing): si se corta,
   correr el script de nuevo retoma donde quedó, tarea por tarea.
4. **`R/03_ensayo1_resumen.R`** — calcula sesgo y cobertura de z por nivel
   de ρ a partir de `results/ensayo1/*.rds`. Escribe `tablas/ensayo1_*.csv`.
5. **`R/04_ensayo2_datos_reales.R`** — ajusta D2 y D3 con el modelo
   mecanístico (`funresMech:::`) y con el modelo beta-binomial propio
   (`fun/betabinom_model.R`), **en paralelo**. El ajuste PRINCIPAL (punto
   estimado de z) usa la configuración real del manuscrito sin reducir
   (`n_sim=3000`, `NP=40`) — es el número citable. El perfil de z para
   el IC95% usa una grilla más gruesa (paso 0.2, 13 puntos en vez de 26)
   y `n_sim=1500` — ver "Cómputo e historial de tiempos". Cada ajuste
   principal y cada punto del perfil se guarda por separado en
   `results/ensayo2/` (checkpointing). Escribe
   `tablas/ensayo2_D2_D3_resultado.csv`.
6. **`R/05_sintesis_y_figuras.R`** — genera las figuras finales en `figuras/`.
7. **`report/validacion_funresMech.Rmd`** — arma el reporte HTML final a
   partir de las tablas y figuras generadas.

## Notas sobre la estructura (respecto de la lista original del plan)

La carpeta `fun/` termina con **tres** archivos en vez de uno solo:

- `betabinom_model.R` — la única pieza de *modelo estadístico* nueva (el
  paquete no tiene beta-binomial, confirmado leyendo el código fuente).
- `simulate_ensayo1.R` — generador de datos sintéticos del Ensayo 1; se
  separó para no duplicar código entre `00_benchmark.R` y
  `02_ensayo1_simulacion.R`.
- `profile_utils.R` — réplica exacta de cómo `R/server.R` del paquete
  calcula el IC95% de z y el AIC_full vs. AIC_restricted (ver más abajo).

`renv.lock` se incluye como un placeholder JSON mínimo (sin paquetes
listados) porque un lockfile real depende de las versiones exactas
instaladas en esta máquina — no tiene sentido inventarlo. Correr
`renv::init()` + `renv::snapshot()` la primera vez para generarlo de verdad,
tal como se indica en "Instalación".

## Hechos verificados contra el código fuente real (github.com/Segon03/funresMech, v1.0.4)

Todo el código de este proyecto está escrito contra el código real del
paquete, no contra supuestos. Puntos clave verificados leyendo
`R/simulate_trial.R`, `R/simulate_distribution.R`, `R/negloglik_fixed_z.R`,
`R/fit_full.R`, `R/server.R`, `R/ui.R`, `NAMESPACE` y
`tests/testthat/test-modelos.R`:

- **Ninguna función mecanística está exportada.** Sólo `run_app()` lo está.
  Por eso todo este proyecto llama a las funciones internas con tres dos
  puntos (`funresMech:::fit_full(...)`, etc.) — exactamente como lo hace
  el propio archivo de tests del paquete.
- **Nombres de columnas:** `dens` (densidad) y `par` (nº parasitados), no
  `x`/`y`.
- **`fit_full(data_spp, T_exp, itermax, NP, reltol, n_sim_profile)`** no
  tiene valores por defecto: los 6 argumentos son obligatorios.
- **Límites internos del optimizador** (hard-coded, no configurables desde
  la app): a∈[0.001,2], h∈[0.001,0.5], z∈[0.5,3.0] (sólo en el ajuste
  completo de 5 parámetros), k∈[0.5,5.0], s∈[0.001,0.5]. z=1, a=0.3,
  h=0.02, k=1 del Ensayo 1 caen cómodamente adentro; s=0 no es alcanzable
  exactamente como estimación (el piso es 0.001).
- **Parametrización de la Lognormal ya resuelta en el código:**
  `meanlog = log(h) - 0.5*s^2`, `sdlog = s`; si `s=0` no se usa `rlnorm`
  en absoluto, se usa `h` directo. No hace falta reinterpretarla.
- **No existe un modelo beta-binomial en el paquete** (no está en
  `DESCRIPTION`/`Imports` ni en ningún archivo de `R/`). El mecanismo de
  comparación de modelos que sí trae la app es interno al modelo
  mecanístico: AIC con z libre (5 parámetros) vs. AIC con z=1 fijo
  (4 parámetros, con el NLL en z=1 interpolado del perfil). Se replicó
  esa lógica en `fun/profile_utils.R::aic_full_vs_restricted()` como
  chequeo adicional "gratis" (sale del mismo perfil que ya hace falta
  calcular para el IC de z).
- **Atención — posible discrepancia a confirmar con el autor antes de
  publicar los IC en el paper:** el umbral de verosimilitud que usa
  `R/server.R` para el IC95% de z es `min_nll + qchisq(0.95, 1)`, **sin**
  dividir por 2. El criterio estándar de razón de verosimilitudes es
  `NLL(z) <= NLL(z_hat) + qchisq(0.95,1)/2`. Omitir el "/2" da un
  intervalo más ancho que el estándar. `fun/profile_utils.R` replica el
  código del paquete tal cual está por defecto (`use_standard_factor =
  FALSE`), pero permite calcular también la versión estándar
  (`use_standard_factor = TRUE`) para comparar antes de decidir cuál
  usar en el paper.
- **Unidades de tiempo:** la app etiqueta `T_exp` como horas, pero el
  modelo es agnóstico a la unidad — lo que importa es que a, h y T estén
  en la misma unidad. Tanto el manuscrito como Okuyama (2026) trabajan
  con T=1 "unidad de tiempo". Por eso todos los scripts usan `T_exp = 1`
  (incluso para los datos reales de 24 h), nunca `T_exp = 24` — ver el
  comentario en `04_ensayo2_datos_reales.R`.

## Cómputo e historial de tiempos (importante)

Una corrida real de `00_benchmark.R` con la configuración ORIGINAL
(`n_sim=1000`, `NP=40`, `itermax=50`) midió **121.9 minutos para un solo
ajuste**. Eso hacía inviable el diseño original tal cual estaba escrito:
Ensayo 1 completo (60 ajustes) ≈ 122 h, Ensayo 2 (perfil completo de 26
puntos x 2 datasets) ≈ 317 h, ambos secuenciales.

En respuesta a eso, esta versión del proyecto:

- Reduce `n_sim`/`NP` para el barrido de 60 ajustes del Ensayo 1
  (`n_sim=300`, `NP=20`).
- Calcula el IC95% de z (perfil local) sólo para 3 réplicas de control
  por nivel de ρ, no las 60 — al implementar el cómputo en paralelo se
  detectó que ese perfil (originalmente calculado para las 60 réplicas)
  multiplicaba por ~10 el costo total del Ensayo 1, algo que la primera
  versión de `00_benchmark.R` no reflejaba (sólo cronometraba el ajuste
  principal). El sesgo y la dispersión de z-hat — el resultado central
  del ensayo — sí se calculan sobre las 60 réplicas completas.
- Para el Ensayo 2, mantiene `n_sim=3000`/`NP=40` (configuración real del
  manuscrito, sin reducir) para el AJUSTE PRINCIPAL de cada dataset — es
  el número que se cita en el paper — pero usa una grilla más gruesa
  (paso 0.2, 13 puntos en vez de 26) y `n_sim=1500` sólo para el PERFIL
  de z (afecta el ancho/forma del IC, no el punto estimado).
- Paraleliza ambos ensayos con `future`/`future.apply`
  (`plan(multisession, workers = parallel::detectCores() - 1)`),
  manteniendo el checkpointing por tarea (cada ajuste o punto de perfil
  se guarda en su propio archivo apenas termina, así que un corte no
  hace perder el trabajo ya hecho).

`00_benchmark.R` mide directamente un ajuste con la configuración YA
REDUCIDA del Ensayo 1 y extrapola (por fórmula de costo, no medición
directa) cuánto tardaría cada ensayo completo, en paralelo, con los
cores detectados en esta máquina — correrlo primero y mirar esa
estimación antes de lanzar `02_ensayo1_simulacion.R` o
`04_ensayo2_datos_reales.R`. Si igual da muy alto, `n_sim_profile` (en
`04_ensayo2_datos_reales.R`) y `n_ci_subset` (en
`02_ensayo1_simulacion.R`) son las dos palancas más directas para bajarlo
más, sin tocar el punto estimado principal de ninguno de los dos
ensayos.

## Estructura de carpetas

```
funresMech_validacion/
|-- funresMech_validacion.Rproj
|-- renv.lock                    # placeholder, ver arriba
|-- README.md
|-- data_raw/                    # CSV de FoRAGE (ya copiado)
|-- data_clean/                  # D2.rds, D3.rds (generados por 01)
|-- R/                           # scripts 00-05, correr en orden
|-- fun/                         # funciones reutilizables (ver arriba)
|-- results/
|   |-- ensayo1/                 # checkpoints por (rho, replica) del Ensayo 1
|   |-- ensayo2/                 # checkpoints por ajuste/punto de perfil del Ensayo 2
|   `-- ensayo2_mech_*.rds, ensayo2_bb_*.rds   # resumen reensamblado (usa 05)
|-- tablas/                      # CSV finales
|-- figuras/                     # PNG para el paper
`-- report/                      # Rmd + HTML final
```
