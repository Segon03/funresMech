# ============================================================
# R/02b_ensayo1_confirmacion.R
#
# Diagnostico: los resultados de 02_ensayo1_simulacion.R (n_sim=300,
# NP=20) mostraron un sesgo sistematico de z_hat por debajo de 1,
# que empeora con rho (ver tablas/ensayo1_resumen_por_rho.csv:
# sesgo_z = -0.10, -0.24, -0.34 para rho = 0, 0.05, 0.10). Ademas
# hay mucha dispersion entre replicas incluso en rho=0 (SD=0.21),
# con casos como rho=0/rep=12 cayendo lejos de los parametros
# verdaderos (a_hat=0.04 vs a=0.3 real).
#
# Sospecha: NP=20 esta por debajo de lo que el propio DEoptim
# recomienda para 5 parametros (el aviso "NP deberia ser >= 10x los
# parametros" que aparece en cada corrida implica NP >= 50), y
# n_sim=300 es 10x menor que el n_sim=3000 del manuscrito. Con un
# optimizador asi de chico, es esperable que caiga en optimos
# locales distintos segun la replica, inflando sesgo y dispersion
# sin que eso sea una propiedad real del modelo.
#
# Este script vuelve a ajustar, con NP=50 y n_sim=1500 (a mitad de
# camino entre el barrido liviano y el n_sim=3000 del manuscrito -
# ver nota de costo mas abajo), EXACTAMENTE los mismos datasets
# simulados que ya se ajustaron en 02_ensayo1_simulacion.R para las
# replicas 1 y 2 de cada nivel de rho (6 ajustes en total - se
# reusan las mismas semillas, asi que el dataset es identico byte a
# byte y la comparacion aisla el efecto de NP/n_sim). Si z_hat se
# acerca a 1 con esta configuracion, confirma que el sesgo era un
# artefacto de la configuracion liviana, no del modelo. Si el sesgo
# persiste, es una senal real que hay que documentar.
#
# NOTA DE COSTO: con la formula de 00_benchmark.R (costo ~ NP x
# n_sim x (itermax+1)), y usando el tiempo real medido ahi
# (13.3 min con NP=20/n_sim=300), un ajuste con NP=50/n_sim=1500
# sale ~12.5x mas caro: ~166 min (2.8 h) por ajuste. Con 6 ajustes
# en paralelo (menos que los workers disponibles), el bloque
# completo deberia tardar ese mismo orden, ~2.5-3 h, no la suma.
#
# IMPORTANTE: este script NO reemplaza los resultados de
# 02_ensayo1_simulacion.R (se guardan aparte, en
# results/ensayo1_confirmacion/) y no recalcula el perfil local de
# z (solo el punto estimado) para mantener el costo acotado - el
# objetivo aca es diagnostico, no una nueva tanda completa de IC.
# ============================================================

library(funresMech)
library(future)
library(future.apply)
source(file.path("fun", "betabinom_model.R"))
source(file.path("fun", "simulate_ensayo1.R"))

dir.create(file.path("results", "ensayo1_confirmacion"), showWarnings = FALSE, recursive = TRUE)
dir.create("tablas", showWarnings = FALSE)

rho_grid     <- c(0, 0.05, 0.1)     # igual que 02_ensayo1_simulacion.R
confirm_reps <- c(1, 2)             # replicas a re-ajustar (mismas para las que ya
                                     # tenemos resultado con la config liviana)
T_exp        <- 1

NP_confirm      <- 50    # cumple la recomendacion de DEoptim (>=10x 5 parametros)
n_sim_confirm   <- 1500  # 5x mas que el barrido liviano (300), mitad del manuscrito (3000)
itermax_confirm <- 50    # sin cambios
reltol_confirm  <- 1e-2  # sin cambios

grid <- expand.grid(rho = rho_grid, rep_id = confirm_reps)
cat(sprintf("Confirmacion Ensayo 1: %d ajustes a correr (NP=%d, n_sim=%d)\n",
            nrow(grid), NP_confirm, n_sim_confirm))

n_cores   <- parallel::detectCores()
n_workers <- max(1, n_cores - 1)
cat(sprintf("Paralelizando en %d workers (de %d cores detectados).\n", n_workers, n_cores))
plan(multisession, workers = n_workers)

ajustar_confirmacion <- function(i) {
  library(funresMech)
  source(file.path("fun", "betabinom_model.R"))
  source(file.path("fun", "simulate_ensayo1.R"))

  rho    <- grid$rho[i]
  rep_id <- grid$rep_id[i]
  out_file <- file.path("results", "ensayo1_confirmacion",
                         sprintf("rho%.2f_rep%02d.rds", rho, rep_id))

  if (file.exists(out_file)) {
    return(sprintf("[%d/%d] rho=%.2f rep=%02d ya existe, se salta", i, nrow(grid), rho, rep_id))
  }

  # MISMA semilla y MISMA llamada que 02_ensayo1_simulacion.R -> mismo dataset
  # simulado exacto, para que la comparacion aisle el efecto de NP/n_sim.
  seed <- 1000 + rep_id + match(rho, rho_grid) * 1000
  set.seed(seed)
  df <- simulate_dataset_ensayo1(rho = rho, T = T_exp)

  t0 <- Sys.time()
  fit <- funresMech:::fit_full(
    data_spp = df, T_exp = T_exp,
    itermax = itermax_confirm, NP = NP_confirm, reltol = reltol_confirm,
    n_sim_profile = n_sim_confirm
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  saveRDS(
    list(rho = rho, rep_id = rep_id, seed = seed, fit = fit, elapsed_sec = elapsed,
         NP = NP_confirm, n_sim = n_sim_confirm),
    out_file
  )

  sprintf("[%d/%d] rho=%.2f rep=%02d -> z_hat=%.3f a_hat=%.3f h_hat=%.4f (%.1f min)",
          i, nrow(grid), rho, rep_id, fit$par["z"], fit$par["a"], fit$par["h"], elapsed / 60)
}

mensajes <- future_lapply(seq_len(nrow(grid)), ajustar_confirmacion, future.seed = TRUE)
invisible(lapply(mensajes, cat, "\n"))

# ---- Comparacion contra los resultados originales (config liviana) -----
cat("\nArmando tabla comparativa contra ensayo1_ajustes_individuales.csv...\n")

original <- read.csv(file.path("tablas", "ensayo1_ajustes_individuales.csv"))

filas <- lapply(seq_len(nrow(grid)), function(i) {
  rho    <- grid$rho[i]
  rep_id <- grid$rep_id[i]
  f <- file.path("results", "ensayo1_confirmacion",
                  sprintf("rho%.2f_rep%02d.rds", rho, rep_id))
  if (!file.exists(f)) return(NULL)
  r <- readRDS(f)

  orig <- original[original$rho == rho & original$rep_id == rep_id, ]

  data.frame(
    rho = rho, rep_id = rep_id,
    z_hat_liviano    = if (nrow(orig) == 1) orig$z_hat else NA,
    z_hat_confirm    = unname(r$fit$par["z"]),
    a_hat_liviano    = if (nrow(orig) == 1) orig$a_hat else NA,
    a_hat_confirm    = unname(r$fit$par["a"]),
    h_hat_liviano    = if (nrow(orig) == 1) orig$h_hat else NA,
    h_hat_confirm    = unname(r$fit$par["h"]),
    k_hat_liviano    = if (nrow(orig) == 1) orig$k_hat else NA,
    k_hat_confirm    = unname(r$fit$par["k"]),
    elapsed_min_confirm = r$elapsed_sec / 60
  )
})
comparacion <- do.call(rbind, filas[!vapply(filas, is.null, logical(1))])

write.csv(comparacion, file.path("tablas", "ensayo1_confirmacion_comparacion.csv"), row.names = FALSE)
cat("\nGuardado en tablas/ensayo1_confirmacion_comparacion.csv\n")
print(comparacion)

cat(sprintf(
  "\nsesgo_z liviano (estas %d replicas): %.3f | sesgo_z confirmacion: %.3f\n",
  nrow(comparacion),
  mean(comparacion$z_hat_liviano - 1, na.rm = TRUE),
  mean(comparacion$z_hat_confirm - 1, na.rm = TRUE)
))
