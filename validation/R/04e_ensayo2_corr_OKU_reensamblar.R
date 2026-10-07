# ############################################################
# REEMPLAZADO por R/04g_ensayo2_corr_OKU_reensamblar_v2.R (25-sep-2026).
# NO usar para el manuscrito. Errores de este script: (1) toma z_hat
# del ajuste con z libre en vez del minimo del perfil (Okuyama usa el
# perfil); (2) IC con umbral +3.84 en vez de +1.92; (3) bb_convergio
# siempre FALSE (identical(0L, 0)); (4) cotas de h desincronizadas con
# 04d; (5) el encabezado dice n_sim=10000 pero 04d uso 3000.
# Se conserva solo como registro historico.
# ############################################################
# ============================================================
# R/04e_ensayo2_corr_OKU_reensamblar.R
#
# Reensambla los resultados de 04d_ensayo2_corr_OKU.R (T=24,
# itermax=200, NP=40, n_sim=10000 - parametros IDENTICOS a los
# de Okuyama 2026, Script S3/S4 de su material suplementario) y
# arma una tabla comparativa DIRECTA, lado a lado, con los
# valores PUBLICADOS por Okuyama en su Tabla S1 (material
# suplementario: papers/Okuyama_Supp_Mat/jen70148-sup-0009-
# tables1.xlsx), para D2 y D3, metodos "ML: beta-binom" y
# "ML: simulation".
#
# Requiere haber corrido 04d_ensayo2_corr_OKU.R primero (o al
# menos que existan los .rds correspondientes en
# results/ensayo2_corr_OKU/).
# ============================================================

source(file.path("fun", "profile_utils.R"))

labs    <- c("D2", "D3")
res_dir <- file.path("results", "ensayo2_corr_OKU")

# Mismas cotas que uso 04d_ensayo2_corr_OKU.R, solo para poder avisar
# si algun parametro convergio pegado a un borde (senal de que la
# cota quedo chica, no de que el ajuste este mal).
lower_ok <- c(a = 1e-6, h = 1e-6, k = 1e-6, s = 0)
upper_ok <- c(a = 5,    h = 5,    k = 5,    s = 3)
eps_borde <- 0.01  # 1% de margen

# ---- Valores publicados por Okuyama (2026), Tabla S1 --------------------
# Copiados a mano desde papers/Okuyama_Supp_Mat/jen70148-sup-0009-
# tables1.xlsx (hoja "estimates"), filas D2 y D3, metodos "ML: beta-
# binom" y "ML: simulation". Revisar contra el Excel si se vuelve a
# necesitar para otro dataset (D1, D4, D5).
okuyama_tabla_s1 <- data.frame(
  dataset = c("D2", "D2", "D3", "D3"),
  metodo  = c("beta-binom", "simulacion", "beta-binom", "simulacion"),
  a       = c(0.020, 0.001, 0.016, 0.012),
  h       = c(0.547, 0.606, 0.611, 0.626),
  z_hat   = c(1.864, 3.100, 2.343, 2.500),
  ci_low  = c(1.310, 1.100, 1.562, 1.400),
  ci_high = c(2.552, 4.600, 3.824, 4.800),
  rho_o_k = c(0.006, 0.037, 0.000, 0.902),
  s       = c(NA,    0.682, NA,    0.542),
  AIC     = c(346.916, 336.589, 273.507, 272.890),
  stringsAsFactors = FALSE
)

leer_rds <- function(f) if (file.exists(f)) readRDS(f) else NULL

avisar_si_en_borde <- function(lab, nombre, valor, lo, hi) {
  if (is.na(valor)) return(invisible())
  rango <- hi - lo
  if (valor <= lo + eps_borde * rango || valor >= hi - eps_borde * rango) {
    cat(sprintf("  AVISO %s: %s=%.4f quedo pegado a una cota de busqueda [%.4g, %.4g] - revisar si hace falta ampliarla.\n",
                lab, nombre, valor, lo, hi))
  }
}

resultados <- list()

for (lab in labs) {

  # --- Mejor de los intentos del ajuste principal (z libre) --------------
  archivos_principal <- list.files(res_dir,
    pattern = sprintf("^principal_%s_corrOKU_intento[0-9]+\\.rds$", lab), full.names = TRUE)
  ajustes_principal <- Filter(Negate(is.null), lapply(archivos_principal, leer_rds))
  if (length(ajustes_principal) == 0) {
    warning(sprintf("%s: todavia no hay ningun ajuste principal corr_OKU - saltando.", lab))
    next
  }
  nlls_principal <- vapply(ajustes_principal, function(f) f$nll, numeric(1))
  fit_mech <- ajustes_principal[[which.min(nlls_principal)]]
  cat(sprintf("%s: %d/%d intentos del ajuste principal listos. NLLs = %s (se usa el minimo).\n",
              lab, length(ajustes_principal),
              length(list.files(res_dir, pattern = sprintf("^principal_%s_corrOKU_intento", lab))),
              paste(sprintf("%.3f", sort(nlls_principal)), collapse = ", ")))

  avisar_si_en_borde(lab, "a", fit_mech$par["a"], lower_ok["a"], upper_ok["a"])
  avisar_si_en_borde(lab, "h", fit_mech$par["h"], lower_ok["h"], upper_ok["h"])
  avisar_si_en_borde(lab, "k", fit_mech$par["k"], lower_ok["k"], upper_ok["k"])
  avisar_si_en_borde(lab, "s", fit_mech$par["s"], lower_ok["s"], upper_ok["s"])

  # --- Perfil de z: todos los archivos presentes --------------------------
  archivos_perfil <- list.files(res_dir,
    pattern = sprintf("^perfil_%s_corrOKU_z.*\\.rds$", lab), full.names = TRUE)
  perfiles <- Filter(Negate(is.null), lapply(archivos_perfil, leer_rds))
  if (length(perfiles) < 3) {
    warning(sprintf("%s: muy pocos puntos de perfil corr_OKU todavia (%d) - saltando IC.", lab, length(perfiles)))
    next
  }
  z_grid_lab  <- vapply(perfiles, function(p) p$z, numeric(1))
  nll_profile <- vapply(perfiles, function(p) p$nll, numeric(1))
  ord <- order(z_grid_lab)
  z_grid_lab  <- z_grid_lab[ord]
  nll_profile <- nll_profile[ord]

  # --- IC95% por verosimilitud-perfil (misma funcion que server.R) -------
  ci_mech <- ci_from_profile_okuyama(z_grid_lab, nll_profile)

  # --- Beta-binomial (T=24) ------------------------------------------------
  fit_bb <- leer_rds(file.path(res_dir, sprintf("bb_%s_corrOKU.rds", lab)))
  if (is.null(fit_bb)) {
    warning(sprintf("%s: todavia no hay ajuste beta-binomial corr_OKU - saltando.", lab))
    next
  }

  aic_mech <- 2 * fit_mech$nll + 2 * 5
  aic_bb   <- 2 * fit_bb$value + 2 * 4

  resultados[[lab]] <- data.frame(
    dataset = lab,
    a_mech = unname(fit_mech$par["a"]), h_mech = unname(fit_mech$par["h"]),
    z_hat_mech = unname(fit_mech$par["z"]),
    ci_low_mech = ci_mech$z_low, ci_high_mech = ci_mech$z_high,
    k_mech = unname(fit_mech$par["k"]), s_mech = unname(fit_mech$par["s"]),
    nll_mech = fit_mech$nll, aic_mech = aic_mech,
    a_bb = unname(fit_bb$par["a"]), h_bb = unname(fit_bb$par["h"]),
    z_bb = unname(fit_bb$par["z"]), rho_bb = unname(fit_bb$par["rho"]),
    nll_bb = fit_bb$value, aic_bb = aic_bb,
    bb_convergio = identical(fit_bb$convergence, 0),
    delta_aic_bb_vs_mech = aic_bb - aic_mech,
    favorece_mecanistico = (aic_bb - aic_mech) > 0,
    ci_excluye_z1 = isTRUE(ci_mech$z_low > 1) || isTRUE(ci_mech$z_high < 1),
    stringsAsFactors = FALSE
  )
}

tabla_corr_OKU <- do.call(rbind, resultados)
write.csv(tabla_corr_OKU, file.path("tablas", "ensayo2_D2_D3_resultado_corr_OKU.csv"), row.names = FALSE)

cat("\n=== Tabla final: funresMech con parametros identicos a Okuyama (T=24, itermax=200, n_sim=10000) ===\n")
print(tabla_corr_OKU)

# --- Tabla comparativa DIRECTA, lado a lado, contra Okuyama (2026) --------
comparacion <- do.call(rbind, lapply(tabla_corr_OKU$dataset, function(lab) {
  fila   <- tabla_corr_OKU[tabla_corr_OKU$dataset == lab, ]
  ok_sim <- okuyama_tabla_s1[okuyama_tabla_s1$dataset == lab & okuyama_tabla_s1$metodo == "simulacion", ]
  ok_bb  <- okuyama_tabla_s1[okuyama_tabla_s1$dataset == lab & okuyama_tabla_s1$metodo == "beta-binom", ]
  data.frame(
    dataset               = lab,
    a_mech_funresMech     = fila$a_mech,           a_mech_okuyama   = ok_sim$a,
    h_mech_funresMech     = fila$h_mech,           h_mech_okuyama   = ok_sim$h,
    z_mech_funresMech     = fila$z_hat_mech,       z_mech_okuyama   = ok_sim$z_hat,
    ci_mech_funresMech    = sprintf("[%.2f, %.2f]", fila$ci_low_mech, fila$ci_high_mech),
    ci_mech_okuyama       = sprintf("[%.2f, %.2f]", ok_sim$ci_low, ok_sim$ci_high),
    aic_mech_funresMech   = fila$aic_mech,         aic_mech_okuyama = ok_sim$AIC,
    aic_bb_funresMech     = fila$aic_bb,           aic_bb_okuyama   = ok_bb$AIC,
    delta_aic_funresMech  = fila$delta_aic_bb_vs_mech,
    delta_aic_okuyama     = ok_bb$AIC - ok_sim$AIC,
    stringsAsFactors = FALSE
  )
}))

write.csv(comparacion, file.path("tablas", "ensayo2_comparacion_vs_okuyama2026.csv"), row.names = FALSE)

cat("\n=== Comparacion directa: funresMech (parametros corr_OKU) vs. Okuyama (2026) Tabla S1 ===\n")
print(comparacion)

cat("\nArchivos generados:\n",
    " - tablas/ensayo2_D2_D3_resultado_corr_OKU.csv\n",
    " - tablas/ensayo2_comparacion_vs_okuyama2026.csv\n",
    "\ndelta_aic > 0 favorece el modelo mecanistico (simulacion) sobre el\n",
    "beta-binomial, en ambas columnas (funresMech y Okuyama) - comparar el\n",
    "SIGNO ademas de la magnitud.\n")
