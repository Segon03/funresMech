# ============================================================================
# 11_comparacion_optimizadores.R
#
# Fase 4 del plan de respuesta a Okuyama:
#   A. Comparación de estrategias de optimización con el MISMO presupuesto de
#      evaluaciones (preset exploratory: 2000 simulaciones por densidad,
#      50 generaciones x 40 individuos), en un punto del perfil (z verdadero):
#        DEoptim | CMA-ES | CRS2 | DEoptim + números aleatorios comunes + pulido
#      3 corridas independientes por estrategia y conjunto de datos.
#      Cada óptimo se reevalúa con 20000 simulaciones (5 veces) para obtener una
#      NLL "precisa"; una corrida tiene éxito si su NLL precisa queda a menos
#      de 0.5 de la mejor encontrada para ese conjunto de datos.
#   B. Estimación de z: perfil sobre grilla (0.5-5, paso 0.5) frente a
#      estimación conjunta + perfil adaptativo. Cobertura del IC, error de z y
#      tiempo.
#
# Escenarios (con verdad conocida, diseño genérico, NO el de D2/D3):
#   z en {1, 2, 3} x k en {0.1, 1, 5}; h = 0.6 h, s = 0.5 h, T = 24 h,
#   a tal que se esperan 20 encuentros por ensayo a la densidad de referencia.
#   Diseño: densidades 2, 4, 8, 12, 16, 24, 32, 48 con 8 réplicas (64 ensayos).
#
# Uso (PC de 23-24 núcleos, R en Windows, desde funresMech_dev):
#   devtools::install()          # versión 1.2.0 de desarrollo
#   source("validation/R/11_comparacion_optimizadores.R")
# Para una prueba corta antes de la corrida larga:  RAPIDO <- TRUE; source(...)
#
# Es reanudable: cada trabajo se guarda en validation/results/fase4/ y, si se
# interrumpe, al volver a correr el script solo hace los que faltan.
# Duración estimada con 22 núcleos: A ~3 h, B ~4 h (total ~7 h).
# Resultados: validation/tablas/fase4_A_optimizadores.csv,
#             validation/tablas/fase4_B_z_conjunto_vs_perfil.csv
#             validation/tablas/fase4_A_resumen_por_escenario.csv
#             validation/tablas/fase4_B_resumen.csv
# ============================================================================

library(funresMech)
library(future)
stopifnot(packageVersion("funresMech") >= "1.2.0")
for (p in c("cmaes", "nloptr")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

if (!exists("RAPIDO")) RAPIDO <- FALSE
N_WORKERS <- if (RAPIDO) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RES <- file.path("validation", "results", if (RAPIDO) "fase4_rapido" else "fase4")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE)
dir.create("validation/tablas", showWarnings = FALSE)

T_EXP <- 24
DENS <- rep(c(2, 4, 8, 12, 16, 24, 32, 48), each = 8)
H_REF <- exp(mean(log(unique(DENS))))
ESC <- expand.grid(z = c(1, 2, 3), k = c(0.1, 1, 5))
if (RAPIDO) ESC <- ESC[c(2, 5), ]
ESC$escenario <- sprintf("z%g_k%g", ESC$z, ESC$k)
ESC$a <- 20 / (T_EXP * H_REF^ESC$z)
ESC$h <- 0.6; ESC$s <- 0.5
N_DATOS_A <- if (RAPIDO) 1 else 5
N_DATOS_B <- if (RAPIDO) 1 else 3
N_REP <- if (RAPIDO) 1 else 3

modelo <- funres_model()
ctrl_base <- if (RAPIDO) {
  funres_control("exploratory", n_sim = 300, itermax = 6, NP = 12, n_reeval = 2, bb = FALSE)
} else {
  funres_control("exploratory", bb = FALSE)
}
estrategias <- list(
  DEoptim     = list(method = "DEoptim"),
  CMAES       = list(method = "CMAES"),
  CRS2        = list(method = "CRS2"),
  DEoptim_crn = list(method = "DEoptim", crn = TRUE, polish = 200))

# ---- conjuntos de datos simulados (semillas fijas) ---------------------------
datos_de <- function(i_esc, j) {
  e <- ESC[i_esc, ]
  set.seed(1000 * i_esc + j)
  s <- simulate_response(modelo, dens = unique(DENS), T = T_EXP,
                         par = c(a = e$a, h = e$h, z = e$z, k = e$k, s = e$s), n_sim = 8)
  data.frame(dens = s$dens, par = s$response)
}

# ---- trabajos ----------------------------------------------------------------
trabajos_A <- do.call(rbind, lapply(seq_len(nrow(ESC)), function(i)
  expand.grid(parte = "A", i_esc = i, dato = seq_len(N_DATOS_A),
              estrategia = names(estrategias), rep = seq_len(N_REP),
              stringsAsFactors = FALSE)))
trabajos_B <- do.call(rbind, lapply(seq_len(nrow(ESC)), function(i)
  expand.grid(parte = "B", i_esc = i, dato = seq_len(N_DATOS_B),
              estrategia = c("perfil", "conjunto"), rep = 1,
              stringsAsFactors = FALSE)))
trabajos <- rbind(trabajos_A, trabajos_B)
trabajos$archivo <- file.path(DIR_RES, sprintf("%s_%s_d%d_%s_r%d.rds", trabajos$parte,
                                               ESC$escenario[trabajos$i_esc], trabajos$dato,
                                               trabajos$estrategia, trabajos$rep))
pendientes <- trabajos[!file.exists(trabajos$archivo), ]
cat(sprintf("Trabajos: %d en total, %d pendientes, %d núcleos.\n",
            nrow(trabajos), nrow(pendientes), N_WORKERS))

correr <- function(tr) {
  library(funresMech)
  future::plan(future::sequential)          # dentro de cada trabajo: secuencial
  e <- ESC[tr$i_esc, ]
  d <- datos_de(tr$i_esc, tr$dato)
  verdad <- c(a = e$a, h = e$h, z = e$z, k = e$k, s = e$s)
  t0 <- proc.time()[["elapsed"]]
  if (tr$parte == "A") {
    est <- estrategias[[tr$estrategia]]
    ctrl <- do.call(funres_control, c(list("exploratory", n_sim = ctrl_base$n_sim,
                                           itermax = ctrl_base$itermax, NP = ctrl_base$NP,
                                           n_reeval = ctrl_base$n_reeval, bb = FALSE), est))
    set.seed(tr$rep * 7919 + tr$i_esc * 101 + tr$dato)
    pt <- funresMech:::.fit_point(e$z, d, T_EXP, modelo, ctrl,
                                  funresMech:::.make_bounds(list(), T_EXP),
                                  funresMech:::.h_ref(d))
    seg <- proc.time()[["elapsed"]] - t0
    nll_prec <- mean(vapply(1:5, function(i)
      funresMech:::.nll_sim(pt$par, modelo, d, T_EXP, if (RAPIDO) 2000 else 20000), 0))
    out <- data.frame(parte = "A", escenario = e$escenario, z_true = e$z, k_true = e$k,
                      dato = tr$dato, estrategia = tr$estrategia, rep = tr$rep,
                      nll_opt = pt$nll_opt, nll_reeval = pt$nll_reeval, nll_precisa = nll_prec,
                      a_rel = pt$par[["a"]] / e$a, h_rel = pt$par[["h"]] / e$h,
                      k_rel = pt$par[["k"]] / e$k, s_rel = pt$par[["s"]] / e$s,
                      en_cota = any(pt$at_bound), nfeval = pt$nfeval, segundos = seg)
  } else {
    set.seed(tr$i_esc * 31 + tr$dato)
    f <- suppressWarnings(fit_funres(
      d, T_EXP, modelo, z_grid = if (RAPIDO) c(1, 2, 3) else seq(0.5, 5, by = 0.5),
      control = ctrl_base, extend_z = TRUE, refine = FALSE,
      z_method = if (tr$estrategia == "conjunto") "joint" else "profile", z_step = 0.25))
    seg <- proc.time()[["elapsed"]] - t0
    ci <- f$ci
    out <- data.frame(parte = "B", escenario = e$escenario, z_true = e$z, k_true = e$k,
                      dato = tr$dato, estrategia = tr$estrategia, rep = 1,
                      z_hat = f$par[["z"]], z_low = ci$z_low, z_high = ci$z_high,
                      cubre = (is.na(ci$z_low) || ci$z_low <= e$z) && (is.na(ci$z_high) || e$z <= ci$z_high),
                      abierto = ci$low_censored || ci$high_censored,
                      excluye_1 = !is.na(ci$z_low) && ci$z_low > 1,
                      n_puntos = nrow(f$profile), fiabilidad = attr(reliability(f), "overall"),
                      segundos = seg)
  }
  saveRDS(out, tr$archivo)
  TRUE
}

if (nrow(pendientes)) {
  plan(multisession, workers = N_WORKERS)
  t_ini <- Sys.time()
  fs <- lapply(seq_len(nrow(pendientes)), function(i)
    future(correr(pendientes[i, ]), seed = TRUE,
           globals = list(correr = correr, pendientes = pendientes, i = i, ESC = ESC,
                          datos_de = datos_de, estrategias = estrategias, modelo = modelo,
                          ctrl_base = ctrl_base, T_EXP = T_EXP, DENS = DENS, RAPIDO = RAPIDO)))
  for (i in seq_along(fs)) {
    ok <- tryCatch(value(fs[[i]]), error = function(e) {
      message("Trabajo ", i, " falló: ", conditionMessage(e)); FALSE })
    if (i %% 10 == 0 || i == length(fs))
      cat(sprintf("[%s] %d/%d trabajos listos\n", format(Sys.time(), "%H:%M"), i, length(fs)))
  }
  plan(sequential)
  cat("Tiempo total:", format(Sys.time() - t_ini), "\n")
}

# ---- tablas ------------------------------------------------------------------
leer <- function(parte) {
  f <- trabajos$archivo[trabajos$parte == parte & file.exists(trabajos$archivo)]
  do.call(rbind, lapply(f, readRDS))
}
A <- leer("A"); B <- leer("B")
sufijo <- if (RAPIDO) "_rapido" else ""
if (!is.null(A)) {
  mejor <- tapply(A$nll_precisa, paste(A$escenario, A$dato), min)
  A$brecha <- A$nll_precisa - mejor[paste(A$escenario, A$dato)]
  A$exito <- A$brecha < 0.5
  write.csv(A, sprintf("validation/tablas/fase4_A_optimizadores%s.csv", sufijo), row.names = FALSE)
  resA <- aggregate(cbind(exito, brecha, en_cota, segundos) ~ estrategia + escenario, A, mean)
  resA_tot <- aggregate(cbind(exito, brecha, en_cota, segundos) ~ estrategia, A, mean)
  cat("\nA. Estrategias de optimización (promedio sobre escenarios):\n"); print(resA_tot, digits = 3)
}
if (!is.null(B)) {
  write.csv(B, sprintf("validation/tablas/fase4_B_z_conjunto_vs_perfil%s.csv", sufijo), row.names = FALSE)
  B$error_z <- abs(B$z_hat - B$z_true)
  resB <- aggregate(cbind(cubre, abierto, error_z, n_puntos, segundos) ~ estrategia, B, mean)
  cat("\nB. Estimación de z:\n"); print(resB, digits = 3)
}
if (!is.null(A)) write.csv(resA, sprintf("validation/tablas/fase4_A_resumen_por_escenario%s.csv", sufijo), row.names = FALSE)
if (!is.null(B)) write.csv(resB, sprintf("validation/tablas/fase4_B_resumen%s.csv", sufijo), row.names = FALSE)
cat("\nListo. Enviar a Claude los archivos validation/tablas/fase4_*.csv\n")
