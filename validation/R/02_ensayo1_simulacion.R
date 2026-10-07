# ============================================================
# R/02_ensayo1_simulacion.R  (Secciones 3 y 5.2 del plan)
#
# Barrido de 60 ajustes (3 niveles de rho x 20 replicas) para
# evaluar la recuperacion de z=1 bajo distintos niveles de sobre-
# dispersion. Guarda CADA ajuste apenas termina (checkpointing):
# si el script se corta (luz, RStudio cerrado, PC colgada), correr
# este mismo script de nuevo retoma exactamente donde quedo (cada
# tarea, corra en el worker que corra, chequea su propio archivo
# antes de calcular nada).
#
# CAMBIOS respecto de la version original (ver 00_benchmark.R para
# el porque): n_sim y NP reducidos para el barrido completo, y el
# perfil LOCAL de z (para el IC95% de cada replica) ya NO se calcula
# para las 60 replicas sino solo para un subconjunto chico de control
# (n_ci_subset) - ver nota mas abajo, es la diferencia mas importante.
#
# IMPORTANTE: correr 00_benchmark.R antes de esto.
# ============================================================

library(funresMech)
library(future)
library(future.apply)
source(file.path("fun", "betabinom_model.R"))
source(file.path("fun", "simulate_ensayo1.R"))
source(file.path("fun", "profile_utils.R"))

dir.create(file.path("results", "ensayo1"), showWarnings = FALSE, recursive = TRUE)

rho_grid    <- c(0, 0.05, 0.1)
n_replicas  <- 20
n_sim_pilot <- 300     # antes 1000 - ver 00_benchmark.R para la justificacion
NP_fit      <- 20      # antes 40
itermax_fit <- 50
T_exp       <- 1       # ver nota de unidades, Seccion 3.4 del plan: NO usar horas

# ---- Perfil LOCAL de z por replica (para el IC95% de esa replica) -------
# NOTA IMPORTANTE (descubierta al implementar el computo paralelo): la
# version original de este script calculaba, para CADA una de las 60
# replicas, un perfil local de 9 puntos (ventana +-0.4 en pasos de 0.1),
# cada punto con su propio ajuste DEoptim completo. Eso multiplica por
# ~10 el costo total del Ensayo 1 (60 ajustes -> ~600 ajustes DEoptim en
# total), algo que NO estaba reflejado en la primera version de
# 00_benchmark.R (esa solo cronometraba el ajuste principal, no el
# perfil). Como el objetivo real del paper es mostrar sesgo/dispersion
# de z-hat entre replicas (eso SI se calcula para las 60), no
# necesariamente cobertura exacta del IC por replica, la cobertura se
# calcula solo sobre un subconjunto chico de "control" por nivel de rho.
# Para volver a calcularla para las 60 replicas (mucho mas caro, ver
# 00_benchmark.R), subir n_ci_subset a n_replicas.
n_ci_subset <- 3        # replicas con perfil local calculado, por nivel de rho
paso_perfil_local <- 0.2  # antes 0.1 (menos puntos por perfil local)

grid <- expand.grid(rho = rho_grid, rep_id = seq_len(n_replicas))

cat(sprintf("Ensayo 1: %d ajustes a correr (rho x replica)\n", nrow(grid)))
cat(sprintf("  De esos, %d replicas por nivel de rho (%d en total) tambien\n",
            n_ci_subset, n_ci_subset * length(rho_grid)))
cat("  calculan un perfil local de z para el IC95%.\n")

lower_prof <- c(a = 0.001, h = 0.001, k = 0.5, s = 0.001)
upper_prof <- c(a = 2.0,   h = 0.5,   k = 5.0, s = 0.5)

n_cores   <- parallel::detectCores()
n_workers <- max(1, n_cores - 1)
cat(sprintf("Paralelizando en %d workers (de %d cores detectados).\n", n_workers, n_cores))
plan(multisession, workers = n_workers)

ajustar_uno <- function(i) {
  # Esta funcion corre en un worker separado: tiene que hacer su propio
  # checkpointing (chequear/crear el archivo) sin depender de estado
  # compartido con el proceso principal.
  library(funresMech)
  source(file.path("fun", "betabinom_model.R"))
  source(file.path("fun", "simulate_ensayo1.R"))
  source(file.path("fun", "profile_utils.R"))

  rho    <- grid$rho[i]
  rep_id <- grid$rep_id[i]
  out_file <- file.path("results", "ensayo1",
                         sprintf("rho%.2f_rep%02d.rds", rho, rep_id))

  if (file.exists(out_file)) {
    return(sprintf("[%d/%d] rho=%.2f rep=%02d ya existe, se salta", i, nrow(grid), rho, rep_id))
  }

  seed <- 1000 + rep_id + match(rho, rho_grid) * 1000
  set.seed(seed)
  df <- simulate_dataset_ensayo1(rho = rho, T = T_exp)

  t0 <- Sys.time()
  fit <- funresMech:::fit_full(
    data_spp = df, T_exp = T_exp,
    itermax = itermax_fit, NP = NP_fit, reltol = 1e-2, n_sim_profile = n_sim_pilot
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  z_hat <- unname(fit$par["z"])
  ci_low <- NA_real_; ci_high <- NA_real_; z_grid <- NA; nll_local <- NA

  if (rep_id <= n_ci_subset) {
    # Perfil LOCAL de z alrededor del z estimado, solo para esta replica
    # "de control" (ver nota arriba sobre por que no se hace para las 60).
    z_grid <- seq(max(0.5, z_hat - 0.4), min(3.0, z_hat + 0.4), by = paso_perfil_local)
    nll_local <- numeric(length(z_grid))
    for (j in seq_along(z_grid)) {
      res <- DEoptim::DEoptim(
        fn = function(par) funresMech:::negloglik_fixed_z(par, z_grid[j], df, T_exp, n_sim_pilot),
        lower = lower_prof, upper = upper_prof,
        control = DEoptim::DEoptim.control(itermax = itermax_fit, NP = NP_fit, reltol = 1e-2, trace = FALSE)
      )
      nll_local[j] <- res$optim$bestval
    }
    ci <- ci_from_profile_okuyama(z_grid, nll_local)
    ci_low  <- ci$z_low
    ci_high <- ci$z_high
  }

  saveRDS(
    list(rho = rho, rep_id = rep_id, seed = seed, fit = fit,
         z_grid = z_grid, nll_local = nll_local,
         ci_low = ci_low, ci_high = ci_high, elapsed_sec = elapsed),
    out_file
  )

  sprintf("[%d/%d] rho=%.2f rep=%02d -> z_hat=%.3f IC=[%s,%s] (%.1f min)",
          i, nrow(grid), rho, rep_id, z_hat,
          ifelse(is.na(ci_low), "NA", sprintf("%.2f", ci_low)),
          ifelse(is.na(ci_high), "NA", sprintf("%.2f", ci_high)),
          elapsed / 60)
}

mensajes <- future_lapply(seq_len(nrow(grid)), ajustar_uno, future.seed = TRUE)
invisible(lapply(mensajes, cat, "\n"))

cat("\nEnsayo 1 completo (o retomado hasta donde estaba). Correr 03_ensayo1_resumen.R.\n")
