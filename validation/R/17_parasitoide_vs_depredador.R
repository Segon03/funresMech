# ============================================================================
# 17_parasitoide_vs_depredador.R
#
# Comparación exploratoria del proceso parasitoide (el hospedador atacado sigue
# disponible y puede reencontrarse) frente al proceso depredador (el atacado se
# retira: equivale a una discriminación perfecta y sin costo de los hospedadores
# ya parasitados) en D2, D2 sin el ensayo atípico (fila 13) y D3.
# Preset exploratory, 2 semillas por combinación (12 ajustes). Solo para ver la
# tendencia del AIC: no es un resultado para informar con precisión.
#
# Uso (desde funresMech_dev, con el paquete 1.2.0-rc2 instalado):
#   rm(list = ls()); RAPIDO <- TRUE; source("validation/R/17_parasitoide_vs_depredador.R")  # ~5 min
#   rm(list = ls()); source("validation/R/17_parasitoide_vs_depredador.R")                   # ~2 h
# Reanudable. Resultados: validation/tablas/proceso_*.csv
# ============================================================================

library(funresMech)
library(future)
stopifnot(packageVersion("funresMech") >= "1.2.0")
if (!exists("RAPIDO")) RAPIDO <- FALSE
DIR <- file.path("validation", "results", if (RAPIDO) "proceso_rapido" else "proceso")
dir.create(DIR, recursive = TRUE, showWarnings = FALSE)
SUF <- if (RAPIDO) "_rapido" else ""
T_EXP <- 24
GRILLA <- if (RAPIDO) c(1, 2, 3) else seq(0.5, 5, by = 0.5)

leer <- function(nm) {
  d <- read.csv(file.path("tests", "testthat", "fixtures", paste0(sub("_sin_atipico", "", nm), ".csv")))
  if (grepl("sin_atipico", nm)) {
    stopifnot(d$dens[13] == 10, d$par[13] == 2)
    d <- d[-13, ]
  }
  d
}
trabajos <- expand.grid(datos = c("D2", "D2_sin_atipico", "D3"), proceso = c("parasitoid", "predator"),
                        semilla = 1:2, stringsAsFactors = FALSE)
if (RAPIDO) trabajos <- trabajos[trabajos$datos == "D3" & trabajos$semilla == 1, ]
trabajos$archivo <- file.path(DIR, sprintf("%s__%s__s%d.rds", trabajos$datos, trabajos$proceso, trabajos$semilla))
pend <- trabajos[!file.exists(trabajos$archivo), ]
cat(sprintf("Ajustes: %d en total, %d pendientes.\n", nrow(trabajos), nrow(pend)))

correr <- function(t) {
  library(funresMech)
  future::plan(future::sequential)
  d <- leer(t$datos)
  ctrl <- if (RAPIDO) funres_control("exploratory", n_sim = 300, itermax = 8, NP = 16, n_reeval = 2, bb = FALSE)
          else funres_control("exploratory", bb = FALSE)
  set.seed(1000 * t$semilla + sum(utf8ToInt(paste(t$datos, t$proceso))))
  t0 <- proc.time()[["elapsed"]]
  f <- suppressWarnings(fit_funres(d, T = T_EXP, model = funres_model(process = t$proceso),
                                   z_grid = GRILLA, control = ctrl, extend_z = FALSE, time_unit = "h"))
  saveRDS(list(trabajo = t, fit = f, segundos = proc.time()[["elapsed"]] - t0), t$archivo)
  TRUE
}
if (nrow(pend)) {
  plan(multisession, workers = min(nrow(pend), max(1L, parallel::detectCores() - 1L)))
  fs <- lapply(seq_len(nrow(pend)), function(i)
    future(correr(pend[i, ]), seed = TRUE, globals = list(correr = correr, pend = pend, i = i, leer = leer,
                                                         T_EXP = T_EXP, GRILLA = GRILLA, RAPIDO = RAPIDO)))
  for (i in seq_along(fs)) tryCatch(value(fs[[i]]), error = function(e)
    message("Falló ", pend$datos[i], " ", pend$proceso[i], ": ", conditionMessage(e)))
  plan(sequential)
}

# ---- tabla ---------------------------------------------------------------------
O <- lapply(trabajos$archivo[file.exists(trabajos$archivo)], readRDS)
R <- do.call(rbind, lapply(O, function(o) {
  f <- o$fit; p <- f$par
  data.frame(datos = o$trabajo$datos, proceso = o$trabajo$proceso, semilla = o$trabajo$semilla,
             ensayos = nrow(f$data), z_hat = p[["z"]], z_low = f$ci$z_low, z_high = f$ci$z_high,
             a = p[["a"]], h = p[["h"]], k = p[["k"]], s = p[["s"]], nll = f$nll, nll_se = f$nll_se,
             AIC = f$aic$AIC_full, fiabilidad = attr(reliability(f), "overall"), horas = o$segundos / 3600)
}))
write.csv(R, sprintf("validation/tablas/proceso_ajustes%s.csv", SUF), row.names = FALSE)
# ΔAIC = AIC(depredador) - AIC(parasitoide): negativo = el depredador ajusta mejor
w <- merge(R[R$proceso == "parasitoid", c("datos", "semilla", "AIC", "z_hat", "h", "k")],
           R[R$proceso == "predator", c("datos", "semilla", "AIC", "z_hat", "h", "k")],
           by = c("datos", "semilla"), suffixes = c("_par", "_dep"))
w$dAIC_dep_menos_par <- w$AIC_dep - w$AIC_par
write.csv(w, sprintf("validation/tablas/proceso_comparacion%s.csv", SUF), row.names = FALSE)
cat("\nΔAIC = AIC(depredador) - AIC(parasitoide); negativo = el depredador ajusta mejor\n")
print(w[, c("datos", "semilla", "dAIC_dep_menos_par", "z_hat_par", "z_hat_dep", "h_par", "h_dep", "k_par", "k_dep")],
      digits = 3, row.names = FALSE)
