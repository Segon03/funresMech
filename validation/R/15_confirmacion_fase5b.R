# ============================================================================
# 15_confirmacion_fase5b.R
#
# Fase 5b: confirmación de funresMech 1.2.0 (candidata rc1) con el método
# corregido en D11 (reevaluación adaptativa del piso de probabilidad,
# refinamiento local de la grilla de z e interpolación por spline monótona).
# Parámetros conocidos; la cobertura del IC de z debe llegar a >= 90%
# (objetivo 95%).
#
#   G1. Núcleo de A de la Fase 5 (18 escenarios: z {1, 1.75, 2.75} x k {0.2, 1, 5}
#       x diseño {grande 64 ensayos, pequeño 30}), los MISMOS 8 conjuntos de
#       datos que A, preset exploratory + verify_fit(). Permite comparar con la
#       Fase 5 sobre datos idénticos (efecto de D11).
#   G2. El mismo núcleo con el preset STANDARD (el que se recomienda para
#       informar) + verify_fit(), 4 conjuntos por escenario. Es la base para
#       recalibrar los umbrales de reliability() (H18).
#   G3. Diseño nuevo: 13 densidades x 5 réplicas (65 ensayos), z x k (9
#       escenarios), 8 conjuntos, exploratory + verify_fit(). Cubre la brecha
#       de diseños de C8.
#
# Uso (PC de 23-24 núcleos, R en Windows, desde funresMech_dev):
#   # reiniciar R (Ctrl+Shift+F10)
#   devtools::install()
#   rm(list = ls())
#   RAPIDO <- TRUE;  source("validation/R/15_confirmacion_fase5b.R")   # prueba (~15 min)
#   rm(list = ls())
#   source("validation/R/15_confirmacion_fase5b.R")                    # corrida larga
#
# Reanudable: cada trabajo se guarda en validation/results/fase5b/. Si se
# interrumpe, volver a correr el script hace solo lo que falta.
# Resultados: validation/tablas/fase5b_*.csv (enviar a Claude).
# ============================================================================

library(funresMech)
library(future)
stopifnot(packageVersion("funresMech") >= "1.2.0")

if (!exists("RAPIDO")) RAPIDO <- FALSE
N_WORKERS <- if (RAPIDO) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RES <- file.path("validation", "results", if (RAPIDO) "fase5b_rapido" else "fase5b")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE)
dir.create("validation/tablas", showWarnings = FALSE)
SUFIJO <- if (RAPIDO) "_rapido" else ""

# ---- diseño (igual que la Fase 5 + el diseño nuevo) ---------------------------
T_EXP <- 24
DISENOS <- list(grande  = rep(c(2, 4, 8, 12, 16, 24, 32, 48), each = 8),   # 64 ensayos
                pequeno = rep(c(2, 6, 12, 24, 48), each = 6),              # 30 ensayos
                trece   = rep(c(2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48), each = 5))  # 65
GRILLA <- if (RAPIDO) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.5)          # = Fase 5
N_G1 <- if (RAPIDO) 1L else 8L
N_G2 <- if (RAPIDO) 1L else 4L
N_G3 <- if (RAPIDO) 1L else 8L
N_PREC <- if (RAPIDO) 2000L else 20000L

ctrl_de <- function(preset, ...) {
  if (RAPIDO) {
    p <- switch(preset, exploratory = list(n_sim = 200, itermax = 5),
                standard = list(n_sim = 300, itermax = 8))
    funres_control(preset, n_sim = p$n_sim, itermax = p$itermax, NP = 12, n_reeval = 2, ...)
  } else funres_control(preset, ...)
}

base <- list(z = 1, k = 1, sh = 0.5, hT = 0.025, E = 20, diseno = "grande",
             gen = "correcto", sigma_a = 0)
esc <- function(...) { e <- modifyList(base, list(...)); as.data.frame(e, stringsAsFactors = FALSE) }
nombre_esc <- function(e) sprintf("%s_z%g_k%g_sh%g_hT%g_E%g_%s", e$gen, e$z, e$k, e$sh, e$hT,
                                  e$E, e$diseno)
Z_A <- c(1, 1.75, 2.75); K_A <- c(0.2, 1, 5)
ESC_NUCLEO <- do.call(rbind, apply(expand.grid(z = Z_A, k = K_A, diseno = c("grande", "pequeno"),
                                               stringsAsFactors = FALSE), 1, function(r)
  esc(z = as.numeric(r[["z"]]), k = as.numeric(r[["k"]]), diseno = r[["diseno"]], bloque = "nucleo")))
ESC_TRECE <- do.call(rbind, apply(expand.grid(z = Z_A, k = K_A), 1, function(r)
  esc(z = r[["z"]], k = r[["k"]], diseno = "trece", bloque = "trece")))
if (RAPIDO) { ESC_NUCLEO <- ESC_NUCLEO[c(2, 14), ]; ESC_TRECE <- ESC_TRECE[5, ] }
ESC <- rbind(ESC_NUCLEO, ESC_TRECE)
ESC$id <- nombre_esc(ESC)
stopifnot(!anyDuplicated(ESC$id))

verdad <- function(e) {
  dens <- DISENOS[[e$diseno]]
  H_ref <- exp(mean(log(unique(dens))))
  c(a = e$E / (T_EXP * H_ref^e$z), h = e$hT * T_EXP, z = e$z, k = e$k, s = e$sh * e$hT * T_EXP)
}
# Mismos datos que la Fase 5 (misma semilla) para los diseños grande y pequeño.
datos_sim <- function(e, j) {
  dens <- DISENOS[[e$diseno]]
  p <- verdad(e)
  set.seed(sum(utf8ToInt(e$id)) * 1000 + j)
  out <- do.call(rbind, lapply(unique(dens), function(x)
    simulate_response(funres_model(), x, T_EXP, p, n_sim = sum(dens == x))))
  data.frame(dens = out$dens, par = out$response)
}

nll_precisa <- function(par, model, d) {
  mean(vapply(1:5, function(i) funresMech:::.nll_sim(par, model, d, T_EXP, N_PREC), 0))
}
aligerar <- function(f) { f$verification_runs <- NULL; f$joint <- NULL; f }
fila_fit <- function(f, z0) {
  ci <- f$ci; p <- f$par; pr <- f$profile
  data.frame(z_hat = p[["z"]], a_hat = p[["a"]], h_hat = p[["h"]], k_hat = p[["k"]],
             s_hat = p[["s"]], z_low = ci$z_low, z_high = ci$z_high,
             abierto_bajo = ci$low_censored, abierto_alto = ci$high_censored,
             cubre = (is.na(ci$z_low) || ci$z_low <= z0) && (is.na(ci$z_high) || z0 <= ci$z_high),
             excluye_1 = !is.na(ci$z_low) && ci$z_low > 1,
             nll = f$nll, nll_se = f$nll_se, dAIC_z1 = f$aic$delta_AIC,
             dAIC_BB = f$aic$delta_AIC_BB, n_puntos = nrow(pr),
             puntos_piso = sum(pr$n_floor > 0, na.rm = TRUE),
             reintentos_piso = sum(pr$n_retry, na.rm = TRUE),
             extendido = f$extended, refinado = f$refined,
             cotas_ampliadas = paste(f$widened, collapse = ","),
             fiabilidad = attr(reliability(f), "overall"))
}
rel_larga <- function(f) {
  r <- as.data.frame(reliability(f))
  data.frame(check = r$check, status = r$status, value = suppressWarnings(as.numeric(r$value)))
}

# ---- lista de trabajos -------------------------------------------------------
# costo en unidades de "punto exploratory de 8 densidades" (~200 s); incluye
# +30% de puntos por la extensión y el refinamiento, solo para ordenar y estimar.
tr <- list()
add <- function(brazo, id, j, preset, costo)
  tr[[length(tr) + 1]] <<- data.frame(brazo = brazo, id = id, dato = j, variante = preset,
                                      costo = costo, stringsAsFactors = FALSE)
pts <- (length(GRILLA) + 3) * 1.3
for (i in seq_len(nrow(ESC))) {
  e <- ESC[i, ]; f_dis <- if (e$diseno == "pequeno") 5 / 8 else 1
  if (e$bloque == "nucleo") {
    for (j in seq_len(N_G1)) add("G1", e$id, j, "exploratory", f_dis * (pts + 7.5))
    for (j in seq_len(N_G2)) add("G2", e$id, j, "standard", f_dis * 5 * (pts + 7.5))
  } else {
    for (j in seq_len(N_G3)) add("G3", e$id, j, "exploratory", 1.1 * (pts + 7.5))
  }
}
trabajos <- do.call(rbind, tr)
trabajos$archivo <- file.path(DIR_RES, gsub("[^A-Za-z0-9_.-]", "_", sprintf(
  "%s__%s__d%s__%s.rds", trabajos$brazo, trabajos$id, trabajos$dato, trabajos$variante)))
trabajos <- trabajos[order(-trabajos$costo), ]
pend <- trabajos[!file.exists(trabajos$archivo), ]

seg_pt <- 200
cat(sprintf("Trabajos: %d en total, %d pendientes; %d núcleos.\n", nrow(trabajos), nrow(pend), N_WORKERS))
if (nrow(pend))
  print(aggregate(cbind(trabajos = 1, horas_nucleo = costo * seg_pt / 3600) ~ brazo, pend, sum), digits = 3)
if (nrow(pend) && RAPIDO) {
  cat("Prueba rápida (RAPIDO): ~15-30 min. Las horas de la tabla corresponden a la corrida completa.\n")
} else if (nrow(pend)) {
  cat(sprintf("Duración estimada: %.0f h (el trabajo más largo: %.1f h). Orientativa: puede variar +-30%%.\n",
              max(sum(pend$costo) * seg_pt / 3600 / N_WORKERS, max(pend$costo) * seg_pt / 3600),
              max(pend$costo) * seg_pt / 3600))
}

correr <- function(t) {
  library(funresMech)
  future::plan(future::sequential)
  t0 <- proc.time()[["elapsed"]]
  semilla <- sum(utf8ToInt(paste(t$brazo, t$id, t$dato, t$variante))) * 7 + 13
  e <- ESC[ESC$id == t$id, ]
  d <- datos_sim(e, t$dato); p0 <- verdad(e); m <- funres_model()
  set.seed(semilla)
  f <- suppressWarnings(fit_funres(d, T_EXP, m, z_grid = GRILLA,
                                   control = ctrl_de(t$variante, bb = TRUE)))
  f <- verify_fit(f)
  res <- fila_fit(f, e$z)
  res$nll_prec_hat <- nll_precisa(f$par, m, d)
  res$nll_prec_true <- nll_precisa(p0, m, d)
  out <- list(res = res, rel = rel_larga(f), fit = aligerar(f))
  out$segundos <- proc.time()[["elapsed"]] - t0
  out$trabajo <- t[, c("brazo", "id", "dato", "variante")]
  saveRDS(out, t$archivo)
  TRUE
}

if (nrow(pend)) {
  plan(multisession, workers = N_WORKERS)
  t_ini <- Sys.time()
  glob <- c("correr", "pend", "ESC", "DISENOS", "T_EXP", "GRILLA", "N_PREC", "RAPIDO", "ctrl_de", "verdad",
            "datos_sim", "nll_precisa", "aligerar", "fila_fit", "rel_larga")
  fs <- vector("list", nrow(pend))
  for (i in seq_len(nrow(pend))) {
    fs[[i]] <- future(correr(pend[i, ]), seed = TRUE, globals = c(mget(glob), list(i = i)))
    cat(sprintf("[%s] iniciado %d/%d: %s %s %s %s\n", format(Sys.time(), "%d %H:%M"), i,
                nrow(pend), pend$brazo[i], pend$id[i], pend$dato[i], pend$variante[i]))
  }
  errores <- 0
  for (i in seq_along(fs)) {
    tryCatch(value(fs[[i]]), error = function(e) {
      errores <<- errores + 1
      msg <- sprintf("%s | %s %s %s %s | %s", format(Sys.time()), pend$brazo[i], pend$id[i],
                     pend$dato[i], pend$variante[i], conditionMessage(e))
      message("Falló: ", msg)
      cat(msg, "\n", file = file.path(DIR_RES, "errores.log"), append = TRUE)
    })
  }
  plan(sequential)
  cat("Tiempo total:", format(Sys.time() - t_ini), "; errores:", errores, "\n")
}

# ---- tablas ------------------------------------------------------------------
prop_ic <- function(x) {
  n <- length(x); p <- mean(x); z <- 1.96
  c(p = p, inf = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    sup = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), n = n)
}
escribir <- function(x, nombre) if (!is.null(x) && nrow(x))
  write.csv(x, sprintf("validation/tablas/fase5b_%s%s.csv", nombre, SUFIJO), row.names = FALSE)
O <- lapply(trabajos$archivo[file.exists(trabajos$archivo)], readRDS)
if (length(O)) {
  G <- do.call(rbind, lapply(O, function(o)
    cbind(o$trabajo, o$res, segundos = o$segundos,
          ESC[ESC$id == o$trabajo$id, c("z", "k", "sh", "hT", "E", "diseno", "bloque")])))
  G$deficit <- G$nll_prec_hat - G$nll_prec_true
  G$fallo_opt <- G$deficit > 1
  G$abierto <- G$abierto_bajo | G$abierto_alto
  G$err_z <- G$z_hat - G$z
  escribir(G, "G_simulacion")
  E <- do.call(rbind, lapply(O, function(o) {
    k <- merge(o$trabajo, o$rel); k$cubre <- o$res$cubre; k$fallo_opt <- (o$res$nll_prec_hat -
      o$res$nll_prec_true) > 1; k$err_z <- o$res$z_hat - ESC$z[ESC$id == o$trabajo$id]; k }))
  escribir(E, "E_chequeos")
  R <- do.call(rbind, lapply(split(G, list(G$brazo, G$z), drop = TRUE), function(q) {
    ic <- prop_ic(q$cubre)
    data.frame(brazo = q$brazo[1], z = q$z[1], n = nrow(q), cobertura = ic[["p"]],
               inf = ic[["inf"]], sup = ic[["sup"]], abiertos = mean(q$abierto),
               excluye_1 = mean(q$excluye_1), err_abs_med = median(abs(q$err_z)),
               fallo_opt = mean(q$fallo_opt), puntos_piso = mean(q$puntos_piso > 0),
               n_puntos = mean(q$n_puntos), horas = mean(q$segundos) / 3600)
  }))
  escribir(R, "G_resumen")
  cat("\nCobertura del IC de z por brazo y z verdadero:\n"); print(R, digits = 3)
  for (b in unique(G$brazo)) {
    ic <- prop_ic(G$cubre[G$brazo == b])
    cat(sprintf("Cobertura global %s: %.1f%% [%.1f, %.1f] (n = %d)\n", b, 100 * ic[["p"]],
                100 * ic[["inf"]], 100 * ic[["sup"]], ic[["n"]]))
  }
  # Comparación con la Fase 5 sobre datos idénticos (G1 frente a A)
  fa <- "validation/tablas/fase5_A_simulacion.csv"
  if (!RAPIDO && file.exists(fa) && any(G$brazo == "G1")) {
    A <- read.csv(fa, stringsAsFactors = FALSE)
    A <- A[A$bloque == "nucleo", ]
    g1 <- G[G$brazo == "G1", ]
    m <- merge(g1[, c("id", "dato", "cubre", "z_low", "z_high", "n_puntos", "segundos", "z_hat")],
               A[, c("id", "dato", "cubre", "z_low", "z_high", "n_puntos", "segundos", "z_hat")],
               by = c("id", "dato"), suffixes = c("_1.2rc", "_fase5"))
    cmp <- data.frame(n = nrow(m), cobertura_fase5 = mean(m$cubre_fase5),
                      cobertura_5b = mean(m$cubre_1.2rc),
                      ic_mas_angosto_5b = mean((m$z_high_1.2rc - m$z_low_1.2rc) <
                                               (m$z_high_fase5 - m$z_low_fase5), na.rm = TRUE),
                      puntos_fase5 = mean(m$n_puntos_fase5), puntos_5b = mean(m$n_puntos_1.2rc),
                      horas_fase5 = mean(m$segundos_fase5) / 3600, horas_5b = mean(m$segundos_1.2rc) / 3600)
    escribir(cmp, "G1_vs_fase5")
    cat("\nG1 frente a la Fase 5 (mismos datos):\n"); print(cmp, digits = 3)
  }
}
