# ============================================================
# R/04c_ensayo2_reensamblar.R
#
# Reensambla tablas/ensayo2_D2_D3_resultado.csv despues de aplicar
# las correcciones de 04b_ensayo2_correcciones.R:
#   - lee TODOS los perfil_{lab}_z*.rds presentes en results/ensayo2/
#     (dinamicamente, por si D3 quedo con mas puntos que D2 tras
#     extender su grilla) en vez de asumir el z_grid original de 13
#     puntos fijo;
#   - usa fit_mech$par["z"] (ajuste principal continuo) como
#     z_hat_mech, no el argmin de la grilla (fix aplicado tambien en
#     04_ensayo2_datos_reales.R para que una corrida nueva desde cero
#     ya salga bien);
#   - vuelve a ajustar el modelo beta-binomial con la version
#     multi-arranque de fit_betabinom() (fun/betabinom_model.R,
#     corregida) para los dos datasets, por consistencia aunque D2
#     ya hubiera convergido bien la primera vez;
#   - agrega aic_restricted_z1 / delta_aic_full_vs_restricted (ya se
#     calculaban, no se exportaban) como evidencia adicional Tipo II
#     vs. Tipo III, independiente de la comparacion con el modelo
#     beta-binomial.
#
# Correr DESPUES de 04_ensayo2_datos_reales.R y 04b_ensayo2_correcciones.R.
# No recalcula ningun ajuste DEoptim - solo lee checkpoints existentes
# y re-ajusta el beta-binomial (rapido, no usa DEoptim).
# ============================================================

library(funresMech)
source(file.path("fun", "betabinom_model.R"))
source(file.path("fun", "profile_utils.R"))

T_exp <- 1
labs  <- c("D2", "D3")

# backup de la tabla anterior, por si hace falta comparar
tabla_previa <- file.path("tablas", "ensayo2_D2_D3_resultado.csv")
if (file.exists(tabla_previa)) {
  file.copy(tabla_previa,
            file.path("tablas", "ensayo2_D2_D3_resultado_ANTES_DE_FIX.csv"),
            overwrite = TRUE)
  cat("Backup de la tabla anterior guardado como ensayo2_D2_D3_resultado_ANTES_DE_FIX.csv\n")
}

resumen_final <- list()

for (lab in labs) {
  data_spp <- readRDS(file.path("data_clean", paste0(lab, ".rds")))

  fit_mech <- readRDS(file.path("results", "ensayo2", sprintf("principal_%s.rds", lab)))

  archivos_perfil <- list.files(
    file.path("results", "ensayo2"),
    pattern = sprintf("^perfil_%s_z.*\\.rds$", lab),
    full.names = TRUE
  )
  puntos <- lapply(archivos_perfil, readRDS)
  z_grid      <- vapply(puntos, function(p) p$z,   numeric(1))
  nll_profile <- vapply(puntos, function(p) p$nll, numeric(1))
  ord <- order(z_grid)
  z_grid <- z_grid[ord]; nll_profile <- nll_profile[ord]

  cat(sprintf("\n%s: %d puntos de perfil (z de %.2f a %.2f)\n",
              lab, length(z_grid), min(z_grid), max(z_grid)))

  # aviso si sigue habiendo outliers evidentes en el perfil (para no
  # repetir a ciegas el problema que tuvo D2 la primera vez)
  mediana_local <- stats::median(nll_profile)
  outliers <- which(nll_profile > mediana_local + 100)
  # (excluye z=0.5, que es real y consistentemente alto en ambos datasets)
  outliers <- outliers[z_grid[outliers] > 0.6]
  if (length(outliers) > 0) {
    cat("ADVERTENCIA: posibles puntos no convergidos todavia en", lab, "en z =",
        paste(z_grid[outliers], collapse = ", "), "\n")
  }

  ci_mech  <- ci_from_profile_okuyama(z_grid, nll_profile)
  aic_pkg  <- aic_full_vs_restricted(fit_mech$nll, z_grid, nll_profile)
  aic_mech <- 2 * fit_mech$nll + 2 * 5

  saveRDS(
    list(fit = fit_mech, z_grid = z_grid, nll_profile = nll_profile,
         ci = ci_mech, aic_pkg = aic_pkg, aic_mech = aic_mech),
    file.path("results", sprintf("ensayo2_mech_%s.rds", lab))
  )

  cat(sprintf("Reajustando modelo beta-binomial para %s (multi-arranque)...\n", lab))
  fit_bb <- fit_betabinom(data_spp, T = T_exp)
  if (fit_bb$convergence != 0) {
    cat(sprintf("ADVERTENCIA: el ajuste beta-binomial de %s SIGUE sin converger (code=%d, %s)\n",
                lab, fit_bb$convergence, fit_bb$message))
  }
  aic_bb <- aic_value(fit_bb$value, k_params = 4)
  saveRDS(fit_bb, file.path("results", sprintf("ensayo2_bb_%s.rds", lab)))

  delta_aic_bb_vs_mech <- aic_bb - aic_mech
  concluye_tipo_III    <- !is.na(ci_mech$z_low) && ci_mech$z_low > 1

  resumen_final[[lab]] <- data.frame(
    dataset = lab,
    z_hat_mech = unname(fit_mech$par["z"]),
    ci_low_mech = ci_mech$z_low, ci_high_mech = ci_mech$z_high,
    z_max_grilla = max(z_grid),
    aic_mech = aic_mech, aic_bb = aic_bb, delta_aic_bb_vs_mech = delta_aic_bb_vs_mech,
    bb_convergio = fit_bb$convergence == 0,
    aic_restricted_z1 = aic_pkg$aic_restricted,
    delta_aic_full_vs_restricted = aic_pkg$delta_aic,
    favorece_mecanistico = delta_aic_bb_vs_mech > 0,
    favorece_tipoIII_vs_tipoII = !is.na(aic_pkg$delta_aic) && aic_pkg$delta_aic > 0,
    ci_excluye_z1 = concluye_tipo_III,
    reproduce_okuyama2026 = (delta_aic_bb_vs_mech > 0) && concluye_tipo_III
  )

  cat(sprintf(
    "%s: z_hat=%.3f IC=[%.3f,%s] | AIC mec=%.1f vs AIC bb=%.1f (delta=%.2f, bb convergio=%s) | AIC full vs restringido delta=%.2f | reproduce Okuyama 2026: %s\n",
    lab, fit_mech$par["z"], ci_mech$z_low,
    if (is.na(ci_mech$z_high)) "NA" else sprintf("%.3f", ci_mech$z_high),
    aic_mech, aic_bb, delta_aic_bb_vs_mech, fit_bb$convergence == 0,
    aic_pkg$delta_aic, resumen_final[[lab]]$reproduce_okuyama2026
  ))
}

tabla_final <- do.call(rbind, resumen_final)
write.csv(tabla_final, file.path("tablas", "ensayo2_D2_D3_resultado.csv"), row.names = FALSE)

cat("\nGuardado en tablas/ensayo2_D2_D3_resultado.csv (version corregida)\n")
print(tabla_final)
