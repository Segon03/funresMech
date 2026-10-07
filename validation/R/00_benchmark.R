# ============================================================
# R/00_benchmark.R  -  CORRER PRIMERO (Seccion 5.1 del plan)
#
# Mide cuanto tarda UN ajuste con la configuracion REAL que van a usar
# 02_ensayo1_simulacion.R y 04_ensayo2_datos_reales.R, y con esa medicion
# extrapola (formula de costo, no medicion directa) cuanto va a tardar
# cada ensayo completo, en paralelo, con los cores disponibles.
#
# HISTORIAL IMPORTANTE: una corrida anterior de este mismo script, con
# la configuracion "grande" original (n_sim=1000, NP=40, itermax=50),
# midio 121.9 MINUTOS para un solo ajuste. Eso hacia inviable el diseno
# original (Ensayo 1 completo ~122 h, Ensayo 2 ~317 h). Los ajustes de
# abajo (n_sim y NP reducidos para el barrido de Ensayo 1, grilla mas
# gruesa + n_sim reducido solo para el PERFIL de Ensayo 2) son la
# respuesta a eso. Este script ahora mide directamente la configuracion
# YA REDUCIDA, en vez de extrapolar desde el numero viejo.
# ============================================================

if (!requireNamespace("funresMech", quietly = TRUE)) {
  stop(
    "funresMech no esta instalado. Correr:\n",
    "  install.packages('remotes')\n",
    "  remotes::install_github('Segon03/funresMech', ref = 'v1.0.4')"
  )
}
if (!requireNamespace("extraDistr", quietly = TRUE)) install.packages("extraDistr")
if (!requireNamespace("DEoptim", quietly = TRUE))    install.packages("DEoptim")
if (!requireNamespace("future", quietly = TRUE))       install.packages("future")
if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")

library(funresMech)   # carga el namespace para que ::: funcione
source(file.path("fun", "betabinom_model.R"))
source(file.path("fun", "simulate_ensayo1.R"))

n_cores   <- parallel::detectCores()
n_workers <- max(1, n_cores - 1)
cat(sprintf("Cores detectados: %d -> se usaran %d workers en paralelo (02 y 04).\n\n",
            n_cores, n_workers))

# ---- Paso 0a: chequeo casi instantaneo -------------------------------
# Mismos valores "chicos" que usa tests/testthat/test-modelos.R del
# propio paquete para pasar los checks de CRAN rapido (verificado en
# el codigo fuente real).
cat("Paso 0a: chequeo rapido con datos de juguete...\n")
data_toy <- data.frame(
  dens = rep(c(5, 10), each = 5),
  par  = c(1, 2, 1, 3, 2, 4, 5, 3, 6, 4)
)

t0 <- Sys.time()
ajuste_toy <- funresMech:::fit_full(
  data_spp = data_toy, T_exp = 1,
  itermax = 20, NP = 20, reltol = 1e-2, n_sim_profile = 50
)
t_toy <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf(
  "  OK en %.1f s. par = a=%.3f h=%.3f z=%.3f k=%.3f s=%.3f\n",
  t_toy, ajuste_toy$par["a"], ajuste_toy$par["h"],
  ajuste_toy$par["z"], ajuste_toy$par["k"], ajuste_toy$par["s"]
))

# ---- Paso 0b: UN ajuste con la configuracion REAL (ya reducida) de
#      Ensayo 1: n_sim=300, NP=20, itermax=50 ---------------------------
cat("\nPaso 0b: un ajuste con n_sim=300, NP=20 (configuracion real de Ensayo 1, ya reducida)...\n")
set.seed(1)
df_test <- simulate_dataset_ensayo1(rho = 0)

t0 <- Sys.time()
fit_test <- funresMech:::fit_full(
  data_spp = df_test, T_exp = 1,
  itermax = 50, NP = 20, reltol = 1e-2, n_sim_profile = 300
)
t_fit <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf(
  "  OK en %.1f min. par = a=%.3f h=%.3f z=%.3f k=%.3f s=%.3f\n",
  t_fit / 60, fit_test$par["a"], fit_test$par["h"],
  fit_test$par["z"], fit_test$par["k"], fit_test$par["s"]
))
cat(sprintf(
  "  (para referencia: con n_sim=1000/NP=40 esto habia medido 121.9 min; si el\n  escalado es ~lineal en NP x n_sim, se esperaba aca ~%.1f min)\n",
  121.9 * (300/1000) * (20/40)
))

# ---- Extrapolacion por formula de costo --------------------------------
# Costo de UN ajuste DEoptim ~ NP x (itermax+1) x n_sim (lineal en cada
# factor: cada generacion evalua NP individuos, cada evaluacion corre
# n_sim simulaciones por nivel de densidad). t_fit de Paso 0b es la
# referencia (NP=20, n_sim=300, itermax=50).
costo_relativo <- function(NP, n_sim, itermax = 50) {
  (NP / 20) * (n_sim / 300) * ((itermax + 1) / 51)
}

# --- Ensayo 1 (02_ensayo1_simulacion.R) ---
# 60 ajustes principales (NP=20, n_sim=300) +
# un subconjunto de replicas (n_ci_subset=3 por rho = 9 de 60) con un
# perfil LOCAL de z de 5 puntos (grilla paso 0.2) para chequear
# cobertura del IC - ver comentario en 02_ensayo1_simulacion.R sobre por
# que NO se hace esto para las 60 replicas (multiplicaria el costo x10).
n_ensayo1_fits    <- 60
n_ci_subset       <- 9   # 3 replicas x 3 niveles de rho
n_puntos_perfil_1 <- 5
t_e1_principal <- n_ensayo1_fits * costo_relativo(20, 300) * t_fit
t_e1_perfil    <- n_ci_subset * n_puntos_perfil_1 * costo_relativo(20, 300) * t_fit
t_e1_seq       <- t_e1_principal + t_e1_perfil

# --- Ensayo 2 (04_ensayo2_datos_reales.R) ---
# 2 ajustes principales "grandes" (NP=40, n_sim=3000, config real del
# manuscrito - SIN reducir, es el numero citable) + grilla de perfil mas
# gruesa (paso 0.2 = 13 puntos) x 2 datasets, con n_sim=1500 (reducido
# solo para el perfil).
n_datasets_e2  <- 2
n_puntos_e2    <- length(seq(0.5, 3.0, by = 0.2))
t_e2_principal <- n_datasets_e2 * costo_relativo(40, 3000) * t_fit
t_e2_perfil    <- n_datasets_e2 * n_puntos_e2 * costo_relativo(40, 1500) * t_fit
t_e2_seq       <- t_e2_principal + t_e2_perfil

cat(sprintf(paste0(
  "\n--- Estimacion aproximada (formula de costo, no medicion directa) ---\n",
  "Ensayo 1 (%d ajustes + %d puntos de perfil en %d replicas de control):\n",
  "  secuencial: ~%.1f h   |   paralelo (%d workers): ~%.1f h\n",
  "Ensayo 2 (2 ajustes grandes n_sim=3000 + %d puntos de perfil n_sim=1500):\n",
  "  secuencial: ~%.1f h   |   paralelo (%d workers): ~%.1f h\n"
),
  n_ensayo1_fits, n_puntos_perfil_1, n_ci_subset,
  t_e1_seq / 3600, n_workers, t_e1_seq / 3600 / n_workers,
  n_datasets_e2 * n_puntos_e2,
  t_e2_seq / 3600, n_workers, t_e2_seq / 3600 / n_workers
))

cat(
  "\nEstos numeros son una extrapolacion por formula (costo ~ NP x (itermax+1)\n",
  "x n_sim), no una medicion directa de cada escenario - tratarlos como orden\n",
  "de magnitud, no como promesa exacta. Si el paralelo no acompana (por RAM,\n",
  "por ejemplo con muchos workers cada uno corriendo DEoptim), el tiempo real\n",
  "puede ser mayor. Si el Ensayo 2 sigue pareciendo muy largo, bajar\n",
  "'n_sim_profile' en 04_ensayo2_datos_reales.R (esta como variable al\n",
  "principio del script, comentada) o 'n_ci_subset' en 02_ensayo1_simulacion.R.\n",
  sep = ""
)
