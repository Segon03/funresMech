# ============================================================
# R/04_ensayo2_datos_reales.R  (Seccion 4 del plan)
#
# Ajusta D2 y D3 con el modelo mecanistico (funresMech:::) y con
# el modelo beta-binomial propio (fun/betabinom_model.R). El
# AJUSTE PRINCIPAL (punto estimado de z) usa la configuracion real
# del manuscrito/app (n_sim=3000, itermax=50, NP=40) sin reducir -
# es el numero que se va a citar en el paper. El PERFIL de z (para
# el IC95%) usa una grilla mas gruesa y un n_sim mas chico (ver
# n_sim_profile abajo) solo para que sea viable en tiempo - ver
# 00_benchmark.R para la justificacion y una estimacion de tiempo.
#
# Cada ajuste principal y cada punto del perfil se guarda en su
# propio archivo dentro de results/ensayo2/ apenas termina
# (checkpointing): si se corta, correr de nuevo retoma solo lo que
# falta. Las tareas corren en paralelo (future.apply).
#
# IMPORTANTE: correr 00_benchmark.R y 01_extract_forage_d2_d3.R antes.
# ============================================================

library(funresMech)
library(future)
library(future.apply)
source(file.path("fun", "betabinom_model.R"))
source(file.path("fun", "profile_utils.R"))

dir.create(file.path("results", "ensayo2"), showWarnings = FALSE, recursive = TRUE)
dir.create("tablas", showWarnings = FALSE)

T_exp   <- 1          # ver nota de unidades (Seccion 4.4 del plan): NO usar 24
itermax <- 50
NP      <- 40         # sin cambios, en todo el Ensayo 2
reltol  <- 1e-2

n_sim_main    <- 3000   # ajuste PRINCIPAL (punto estimado) - config real del
                         # manuscrito, SIN reducir: es el numero citable.
n_sim_profile <- 1500   # SOLO para el perfil de z (antes 3000, igual que el
                         # principal). Si el tiempo de 00_benchmark.R sigue
                         # pareciendo alto, bajar esto (p.ej. a 700-800) es
                         # la palanca mas directa - no afecta el punto
                         # estimado principal, solo el ancho/forma del IC.

z_grid <- seq(0.5, 3.0, by = 0.2)   # antes by = 0.1 (26 puntos -> 13 puntos)

lower_prof <- c(a = 0.001, h = 0.001, k = 0.5, s = 0.001)
upper_prof <- c(a = 2.0,   h = 0.5,   k = 5.0, s = 0.5)

labs <- c("D2", "D3")
datasets <- setNames(lapply(labs, function(lab) {
  data_path <- file.path("data_clean", paste0(lab, ".rds"))
  if (!file.exists(data_path)) {
    stop("No se encontro ", data_path, ". Correr 01_extract_forage_d2_d3.R primero.")
  }
  readRDS(data_path)
}), labs)

# ---- Lista de tareas: 1 ajuste principal + length(z_grid) puntos de
#      perfil, por cada dataset. Cada tarea se identifica por su propio
#      archivo de checkpoint. ------------------------------------------
tareas <- list()
for (lab in labs) {
  tareas[[length(tareas) + 1]] <- list(type = "principal", lab = lab)
  for (z in z_grid) {
    tareas[[length(tareas) + 1]] <- list(type = "perfil", lab = lab, z = z)
  }
}
cat(sprintf("Ensayo 2: %d tareas a correr (2 ajustes principales + %d puntos de perfil)\n",
            length(tareas), length(tareas) - 2))

n_cores   <- parallel::detectCores()
n_workers <- max(1, n_cores - 1)
cat(sprintf("Paralelizando en %d workers (de %d cores detectados).\n", n_workers, n_cores))
plan(multisession, workers = n_workers)

correr_tarea <- function(tarea) {
  library(funresMech)
  source(file.path("fun", "betabinom_model.R"))
  source(file.path("fun", "profile_utils.R"))

  lab <- tarea$lab
  data_spp <- datasets[[lab]]

  if (tarea$type == "principal") {
    out_file <- file.path("results", "ensayo2", sprintf("principal_%s.rds", lab))
    if (file.exists(out_file)) return(sprintf("%s: ajuste principal ya existe, se salta", lab))
    fit <- funresMech:::fit_full(
      data_spp = data_spp, T_exp = T_exp,
      itermax = itermax, NP = NP, reltol = reltol, n_sim_profile = n_sim_main
    )
    saveRDS(fit, out_file)
    return(sprintf("%s: ajuste principal listo (z_hat=%.3f)", lab, fit$par["z"]))
  }

  # tarea$type == "perfil"
  out_file <- file.path("results", "ensayo2",
                         sprintf("perfil_%s_z%.2f.rds", lab, tarea$z))
  if (file.exists(out_file)) return(sprintf("%s z=%.2f: perfil ya existe, se salta", lab, tarea$z))
  res <- DEoptim::DEoptim(
    fn = function(par) funresMech:::negloglik_fixed_z(par, tarea$z, data_spp, T_exp, n_sim_profile),
    lower = lower_prof, upper = upper_prof,
    control = DEoptim::DEoptim.control(itermax = itermax, NP = NP, reltol = reltol, trace = FALSE)
  )
  saveRDS(list(z = tarea$z, nll = res$optim$bestval), out_file)
  sprintf("%s z=%.2f: nll=%.3f", lab, tarea$z, res$optim$bestval)
}

mensajes <- future_lapply(tareas, correr_tarea, future.seed = TRUE)
invisible(lapply(mensajes, cat, "\n"))

# ---- Reensamblar resultados desde los checkpoints ----------------------
resumen_final <- list()

for (lab in labs) {
  data_spp <- datasets[[lab]]

  fit_mech <- readRDS(file.path("results", "ensayo2", sprintf("principal_%s.rds", lab)))
  nll_profile <- vapply(z_grid, function(z) {
    r <- readRDS(file.path("results", "ensayo2", sprintf("perfil_%s_z%.2f.rds", lab, z)))
    r$nll
  }, numeric(1))

  ci_mech  <- ci_from_profile_okuyama(z_grid, nll_profile)
  aic_pkg  <- aic_full_vs_restricted(fit_mech$nll, z_grid, nll_profile)
  aic_mech <- 2 * fit_mech$nll + 2 * 5

  # Se guarda tambien con el nombre/forma originales (results/ensayo2_mech_*.rds)
  # para que 05_sintesis_y_figuras.R siga funcionando sin cambios.
  saveRDS(
    list(fit = fit_mech, z_grid = z_grid, nll_profile = nll_profile,
         ci = ci_mech, aic_pkg = aic_pkg, aic_mech = aic_mech),
    file.path("results", sprintf("ensayo2_mech_%s.rds", lab))
  )

  cat(sprintf("Ajustando modelo beta-binomial para %s...\n", lab))
  fit_bb <- fit_betabinom(data_spp, T = T_exp)
  aic_bb <- aic_value(fit_bb$value, k_params = 4)
  saveRDS(fit_bb, file.path("results", sprintf("ensayo2_bb_%s.rds", lab)))

  delta_aic_bb_vs_mech <- aic_bb - aic_mech   # > 0 favorece al mecanistico
  concluye_tipo_III    <- !is.na(ci_mech$z_low) && ci_mech$z_low > 1

  # FIX (revision post-hoc, ver notas de Claude en el chat): z_hat_mech
  # usaba ci_mech$z_hat, que es el argmin de la grilla GRUESA del perfil
  # (paso 0.2, n_sim=1500) - no el punto estimado continuo del ajuste
  # PRINCIPAL (n_sim=3000, DEoptim libre en 5 parametros), que es el
  # numero que el comentario de cabecera del script dice que se va a
  # citar en el paper. Para D2 esto daba 1.10 en vez de 1.486 (fit_mech);
  # para D3, 2.30 en vez de 2.096. El AIC ya usaba correctamente
  # fit_mech$nll, asi que solo el z_hat/etiqueta estaba desalineado
  # con el AIC reportado junto a el. Se agregan tambien aic_full/
  # aic_restricted/delta_aic_full_vs_restricted (ya calculados en
  # aic_pkg pero antes descartados) porque dan una senal Tipo II vs.
  # Tipo III mas directa que la comparacion contra el modelo beta-
  # binomial (que ademas tuvo un problema de convergencia en D3, ver
  # notas del chat y fun/betabinom_model.R).
  resumen_final[[lab]] <- data.frame(
    dataset = lab,
    z_hat_mech = unname(fit_mech$par["z"]),
    ci_low_mech = ci_mech$z_low, ci_high_mech = ci_mech$z_high,
    aic_mech = aic_mech, aic_bb = aic_bb, delta_aic_bb_vs_mech = delta_aic_bb_vs_mech,
    aic_restricted_z1 = aic_pkg$aic_restricted,
    delta_aic_full_vs_restricted = aic_pkg$delta_aic,
    favorece_mecanistico = delta_aic_bb_vs_mech > 0,
    favorece_tipoIII_vs_tipoII = !is.na(aic_pkg$delta_aic) && aic_pkg$delta_aic > 0,
    ci_excluye_z1 = concluye_tipo_III,
    reproduce_okuyama2026 = (delta_aic_bb_vs_mech > 0) && concluye_tipo_III
  )

  cat(sprintf(
    "%s: z_hat=%.2f IC=[%.2f,%.2f] | AIC mec=%.1f vs AIC bb=%.1f (delta=%.1f) | reproduce Okuyama 2026: %s\n",
    lab, ci_mech$z_hat, ci_mech$z_low, ci_mech$z_high, aic_mech, aic_bb,
    delta_aic_bb_vs_mech, resumen_final[[lab]]$reproduce_okuyama2026
  ))
}

tabla_final <- do.call(rbind, resumen_final)
write.csv(tabla_final, file.path("tablas", "ensayo2_D2_D3_resultado.csv"), row.names = FALSE)

cat("\nGuardado en tablas/ensayo2_D2_D3_resultado.csv\n")
print(tabla_final)
