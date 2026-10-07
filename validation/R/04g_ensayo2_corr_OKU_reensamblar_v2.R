# ============================================================
# R/04g_ensayo2_corr_OKU_reensamblar_v2.R
#
# Reemplaza a 04e_ensayo2_corr_OKU_reensamblar.R. Arma la tabla
# final del Ensayo 2 (corr_OKU) siguiendo el PROCEDIMIENTO EXACTO
# de Okuyama (2026, Sec. 2.3) y la compara con su Tabla S1.
#
# DIFERENCIAS CON 04e (errores corregidos):
#
# 1) Estimador de z del modelo de simulacion. Okuyama NO ajusta z
#    libre: evalua una grilla de z, ajusta (a, h, k, s) condicional
#    a cada z, y toma como MLE el punto de la grilla con menor NLL
#    ("The parameter set associated with the smallest negative
#    log-likelihood over this grid was regarded as the MLE"). 04e
#    usaba el mejor de los 3 ajustes con z libre, que ademas
#    quedaron PEORES que el minimo del perfil (D2: 163.73 vs 162.93).
#    Aca: z_hat = argmin del perfil; AIC = 2*NLL_min + 2*5.
#    Los ajustes con z libre se siguen informando, solo como control.
#
# 2) Umbral del IC95%: NLL_min + qchisq(.95,1)/2 (= +1.92), no
#    + qchisq(.95,1) (= +3.84) como hace la app funresMech v1.0.4.
#    Ver fun/profile_utils.R. Este script VERIFICA el criterio
#    contra Okuyama con el modelo beta-binomial (sin ruido Monte
#    Carlo): el IC de D2 tiene que dar [1.310, 2.552].
#
# 3) IC tambien sobre la grilla discreta (sin interpolar), como
#    Okuyama. Se informan ambas lecturas (interior/exterior, ver
#    ci_grilla_okuyama()) y el interpolado.
#
# 4) Limites censurados: si el perfil no cruza el umbral dentro de
#    la grilla, se informa ">= max(grilla)" en vez de NA mudo.
#
# 5) bb_convergio: 04e usaba identical(convergence, 0), que da FALSE
#    cuando optim devuelve 0L (entero). Aca se usa == 0.
#
# 6) IC del beta-binomial: 04e no lo calculaba. Aca se calcula por
#    perfil (z fijo, a/h/rho optimizados en escala log/logit, sin
#    cotas que corten el perfil a z altos).
#
# 7) Parametro s: funresMech usa s = desvio en escala LOG del tiempo
#    de manipulacion (sdlog); Okuyama usa s = desvio en escala
#    NATURAL (rlnorm0_cpp en su Script S5). Es la misma familia
#    (lognormal con media h), asi que la verosimilitud maxima es la
#    misma, pero los valores de s no se comparan directo. Se agrega
#    s_natural = h * sqrt(exp(s^2) - 1) para poder comparar.
#
# Funciona con o sin la extension de 04f (si falta, avisa que el
# limite superior de D2 queda censurado).
# ============================================================

source(file.path("fun", "profile_utils.R"))
source(file.path("fun", "betabinom_model.R"))

labs    <- c("D2", "D3")
res_dir <- file.path("results", "ensayo2_corr_OKU")
T_exp   <- 24
grilla_esperada <- round(seq(0.5, 5.9, by = 0.2), 2)

# Cotas usadas en 04d/04f (para avisar si un parametro quedo pegado a una)
lower_ok <- c(a = 1e-6, h = 0.01, k = 1e-6, s = 0)
upper_ok <- c(a = 5,    h = 5,    k = 5,    s = 3)

# ---- Okuyama (2026), Tabla S1 (jen70148-sup-0009-tables1.xlsx) ----------
okuyama <- data.frame(
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
ok_de <- function(lab, met) okuyama[okuyama$dataset == lab & okuyama$metodo == met, ]

leer <- function(f) if (file.exists(f)) readRDS(f) else NULL
fmt_ci <- function(lo, hi, lo_c = FALSE, hi_c = FALSE, d = 2) {
  f <- paste0("%.", d, "f")
  sprintf("[%s, %s]",
          if (lo_c) paste0("<=", sprintf(f, lo)) else sprintf(f, lo),
          if (hi_c) paste0(">=", sprintf(f, hi)) else sprintf(f, hi))
}
avisar_borde <- function(lab, par) {
  for (p in intersect(names(par), names(lower_ok))) {
    # cerca del piso: a menos de 1e-4 del rango (a y k pueden ser
    # legitimamente chicos, p.ej. Okuyama D2: a=0.001, k=0.037);
    # cerca del techo: a menos del 1% del rango.
    r <- upper_ok[p] - lower_ok[p]
    if (par[p] <= lower_ok[p] + 1e-4 * r || par[p] >= upper_ok[p] - 0.01 * r)
      cat(sprintf("  AVISO %s: %s = %.4g quedo pegado a una cota [%g, %g].\n",
                  lab, p, par[p], lower_ok[p], upper_ok[p]))
  }
}

# ---- Perfil beta-binomial (z fijo) ---------------------------------------
# Optimiza (log a, log h, logit rho) sin cotas: a z alto el 'a' optimo
# es muy chico (a*x^z tiene que mantenerse), y una cota inferior fija
# (como el 1e-4 de fit_betabinom) cortaria el perfil.
perfil_bb_en_z <- function(z, d, T, inicio) {
  f <- function(q) negloglik_betabinom(c(exp(q[1]), exp(q[2]), z, plogis(q[3])), d, T)
  mejores <- lapply(inicio, function(st) {
    q0 <- c(log(st[1]), log(st[2]), qlogis(min(max(st[3], 1e-6), 0.9)))
    r1 <- optim(q0, f, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12))
    r2 <- tryCatch(optim(r1$par, f, method = "BFGS", control = list(maxit = 1000, reltol = 1e-12)),
                   error = function(e) r1)
    if (r2$value <= r1$value) r2 else r1
  })
  m <- mejores[[which.min(vapply(mejores, `[[`, numeric(1), "value"))]]
  list(nll = m$value, par = c(a = exp(m$par[1]), h = exp(m$par[2]), rho = plogis(m$par[3])))
}

perfil_bb <- function(d, T, fit_bb, z_grid) {
  a0 <- unname(fit_bb$par[1]); h0 <- unname(fit_bb$par[2]); z0 <- unname(fit_bb$par[3])
  xm <- median(d$dens)
  prev <- NULL
  out <- vector("list", length(z_grid))
  for (i in seq_along(z_grid)) {
    z <- z_grid[i]
    a_esc <- a0 * xm^(z0 - z)   # mantiene a*x^z en la escala del ajuste libre
    inicios <- list(c(a_esc, h0, 0.01), c(a_esc, h0, 0.2))
    if (!is.null(prev)) inicios <- c(inicios, list(prev))
    r <- perfil_bb_en_z(z, d, T, inicios)
    prev <- r$par
    out[[i]] <- r
  }
  data.frame(z = z_grid, nll = vapply(out, `[[`, numeric(1), "nll"))
}

# IC exacto del beta-binomial: busca el cruce en la grilla fina y lo
# refina con uniroot sobre la funcion perfil real (no interpolada).
ic_bb_exacto <- function(d, T, fit_bb, prof, umbral_delta = UMBRAL_IC95) {
  nll_min <- min(fit_bb$value, min(prof$nll))
  thr <- nll_min + umbral_delta
  a0 <- unname(fit_bb$par[1]); h0 <- unname(fit_bb$par[2]); z0 <- unname(fit_bb$par[3])
  xm <- median(d$dens)
  g <- function(z) perfil_bb_en_z(z, d, T, list(c(a0 * xm^(z0 - z), h0, 0.01),
                                                c(a0 * xm^(z0 - z), h0, 0.2)))$nll - thr
  i_min <- which.min(prof$nll)
  izq <- which(prof$nll[seq_len(i_min)] > thr)
  der <- which(prof$nll[i_min:nrow(prof)] > thr)
  lo <- if (length(izq)) uniroot(g, c(prof$z[max(izq)], prof$z[max(izq) + 1]), tol = 1e-6)$root else NA_real_
  hi <- if (length(der)) { j <- i_min - 1 + min(der)
                          uniroot(g, c(prof$z[j - 1], prof$z[j]), tol = 1e-6)$root } else NA_real_
  c(low = lo, high = hi)
}

# ============================================================
resultados <- list(); comparacion <- list(); perfiles_largo <- list()
verificacion_ok <- NA

for (lab in labs) {
  cat(sprintf("\n==================== %s ====================\n", lab))
  d <- readRDS(file.path("data_clean", paste0(lab, ".rds")))

  # --- Perfil del modelo de simulacion -----------------------------------
  archivos <- list.files(res_dir, pattern = sprintf("^perfil_%s_corrOKU_z.*\\.rds$", lab), full.names = TRUE)
  perf <- lapply(archivos, readRDS)
  prof <- data.frame(z = round(vapply(perf, function(p) as.numeric(p$z), 0), 2),
                     nll = vapply(perf, function(p) as.numeric(p$nll), 0))
  pars_perfil <- lapply(perf, function(p) p$par)   # NULL en los puntos de 04d
  o <- order(prof$z); prof <- prof[o, ]; pars_perfil <- pars_perfil[o]
  faltan <- setdiff(grilla_esperada, prof$z)
  cat(sprintf("Perfil simulacion: %d puntos (z = %.1f a %.1f).\n", nrow(prof), min(prof$z), max(prof$z)))
  if (length(faltan)) cat(sprintf("  Faltan puntos de la grilla 0.5-5.9: %s\n", paste(faltan, collapse = ", ")))
  if (any(diff(prof$z) > 0.21)) cat("  AVISO: hay huecos en la grilla (paso > 0.2).\n")

  i_min   <- which.min(prof$nll)
  z_hat   <- prof$z[i_min]
  nll_min <- prof$nll[i_min]
  aic_sim <- 2 * nll_min + 2 * 5

  ci_int  <- ci_from_profile_okuyama(prof$z, prof$nll)                               # +1.92
  ci_app  <- ci_from_profile_okuyama(prof$z, prof$nll, use_standard_factor = FALSE)  # +3.84 (app v1.0.4)
  ci_gr   <- ci_grilla_okuyama(prof$z, prof$nll)

  # --- Ajustes con z libre (solo control) ----------------------------------
  princ <- Filter(Negate(is.null), lapply(
    list.files(res_dir, pattern = sprintf("^principal_%s_corrOKU_intento[0-9]+\\.rds$", lab), full.names = TRUE), leer))
  nll_princ <- vapply(princ, function(f) as.numeric(f$nll), 0)
  best_pr <- princ[[which.min(nll_princ)]]
  cat(sprintf("Ajustes con z libre (control): NLL = %s | z = %s\n",
              paste(sprintf("%.2f", nll_princ), collapse = ", "),
              paste(sprintf("%.2f", vapply(princ, function(f) f$par[["z"]], 0)), collapse = ", ")))
  if (min(nll_princ) > nll_min)
    cat(sprintf("  El minimo del perfil (%.2f en z=%.1f) es MEJOR que el mejor ajuste con z libre (%.2f):\n",
                nll_min, z_hat, min(nll_princ)),
        "  esperable (Okuyama: optimizar 5 parametros es menos estable que 4). Se usa el perfil.\n")

  # --- Parametros en el minimo -------------------------------------------
  if (!is.null(pars_perfil[[i_min]])) {
    par_min <- c(pars_perfil[[i_min]][c("a", "h")], z = z_hat, pars_perfil[[i_min]][c("k", "s")])
    fuente_par <- sprintf("perfil en z=%.1f", z_hat)
  } else {
    par_min <- best_pr$par[c("a", "h", "z", "k", "s")]
    fuente_par <- sprintf("mejor ajuste con z libre (z=%.2f); el punto z=%.1f del perfil no guardo parametros",
                          best_pr$par[["z"]], z_hat)
  }
  avisar_borde(lab, par_min)
  s_nat <- unname(par_min["h"] * sqrt(exp(par_min["s"]^2) - 1))

  # --- Beta-binomial --------------------------------------------------------
  fit_bb <- readRDS(file.path(res_dir, sprintf("bb_%s_corrOKU.rds", lab)))
  aic_bb <- 2 * fit_bb$value + 2 * 4
  cat("Perfil beta-binomial (grilla fina 0.5-6.0, paso 0.01)... ")
  prof_bb <- perfil_bb(d, T_exp, fit_bb, round(seq(0.5, 6.0, by = 0.01), 2))
  if (min(prof_bb$nll) < fit_bb$value - 1e-3)
    cat(sprintf("\n  AVISO: el perfil bb encontro NLL %.4f < ajuste bb %.4f (el ajuste libre no llego al optimo).\n",
                min(prof_bb$nll), fit_bb$value))
  ic_bb     <- ic_bb_exacto(d, T_exp, fit_bb, prof_bb)
  ic_bb_app <- ic_bb_exacto(d, T_exp, fit_bb, prof_bb, umbral_delta = qchisq(0.95, 1))
  cat("listo.\n")

  # --- Verificacion del criterio del IC contra Okuyama (sin ruido MC) ------
  okb <- ok_de(lab, "beta-binom")
  cat(sprintf("IC beta-binomial:  +1.92 -> [%.3f, %.3f] | +3.84 (app) -> [%.3f, %.3f] | Okuyama [%.3f, %.3f]\n",
              ic_bb["low"], ic_bb["high"], ic_bb_app["low"], ic_bb_app["high"], okb$ci_low, okb$ci_high))
  if (lab == "D2") {
    verificacion_ok <- all(abs(ic_bb - c(okb$ci_low, okb$ci_high)) < 0.005)
    cat(if (verificacion_ok)
          "  VERIFICACION OK: con el umbral +1.92 el IC beta-binomial de D2 coincide con Okuyama (tolerancia 0.005).\n"
        else
          "  VERIFICACION FALLIDA: el IC beta-binomial de D2 no coincide con Okuyama. Revisar antes de seguir.\n")
  } else {
    cat(sprintf("  (En %s nuestro ajuste bb tiene AIC %.3f vs %.3f de Okuyama: su optim quedo algo lejos del\n",
                lab, aic_bb, okb$AIC),
        "   optimo, lo que desplaza levemente su umbral y su IC. Por eso la verificacion exacta se hace con D2.)\n")
  }

  # --- Resumen en consola -----------------------------------------------------
  oks <- ok_de(lab, "simulacion")
  cat(sprintf("SIMULACION  z_hat = %.1f (Okuyama %.1f) | NLL_min = %.3f | AIC = %.3f (Okuyama %.3f)\n",
              z_hat, oks$z_hat, nll_min, aic_sim, oks$AIC))
  cat(sprintf("  IC grilla interior %s | exterior %s | interpolado %s | Okuyama [%.1f, %.1f]\n",
              fmt_ci(ci_gr$interior_low, ci_gr$interior_high, ci_gr$low_censurado, ci_gr$high_censurado, 1),
              fmt_ci(ci_gr$exterior_low, ci_gr$exterior_high, d = 1),
              fmt_ci(ci_int$z_low, ci_int$z_high),
              oks$ci_low, oks$ci_high))
  cat(sprintf("  (con el umbral viejo de la app, +3.84: %s)\n", fmt_ci(ci_app$z_low, ci_app$z_high)))
  cat(sprintf("DELTA AIC (bb - simulacion) = %.2f (Okuyama %.2f)\n", aic_bb - aic_sim, okb$AIC - oks$AIC))

  ic_lo <- if (ci_int$low_censurado) min(prof$z) else ci_int$z_low
  resultados[[lab]] <- data.frame(
    dataset = lab,
    z_hat_sim = z_hat, nll_min_sim = nll_min, aic_sim = aic_sim,
    ic_sim_interp_low = ci_int$z_low, ic_sim_interp_high = ci_int$z_high,
    ic_sim_grilla_int_low = ci_gr$interior_low, ic_sim_grilla_int_high = ci_gr$interior_high,
    ic_sim_grilla_ext_low = ci_gr$exterior_low, ic_sim_grilla_ext_high = ci_gr$exterior_high,
    ic_sim_low_censurado = ci_gr$low_censurado, ic_sim_high_censurado = ci_gr$high_censurado,
    ic_sim_app_v104_low = ci_app$z_low, ic_sim_app_v104_high = ci_app$z_high,
    a_sim = unname(par_min["a"]), h_sim = unname(par_min["h"]),
    k_sim = unname(par_min["k"]), s_sim_sdlog = unname(par_min["s"]), s_sim_natural = s_nat,
    fuente_parametros = fuente_par,
    nll_mejor_z_libre = min(nll_princ), z_mejor_z_libre = best_pr$par[["z"]],
    a_bb = unname(fit_bb$par[1]), h_bb = unname(fit_bb$par[2]), z_bb = unname(fit_bb$par[3]),
    rho_bb = unname(fit_bb$par[4]), nll_bb = fit_bb$value, aic_bb = aic_bb,
    bb_convergio = isTRUE(fit_bb$convergence == 0),
    ic_bb_low = unname(ic_bb["low"]), ic_bb_high = unname(ic_bb["high"]),
    delta_aic_bb_menos_sim = aic_bb - aic_sim,
    ic_sim_excluye_z1 = ic_lo > 1,
    ic_bb_excluye_z1 = isTRUE(ic_bb["low"] > 1),
    stringsAsFactors = FALSE
  )

  comparacion[[lab]] <- data.frame(
    dataset = lab,
    metodo  = c("simulacion", "beta-binom"),
    z_hat_funresMech = c(z_hat, unname(fit_bb$par[3])),
    z_hat_okuyama    = c(oks$z_hat, okb$z_hat),
    ic_funresMech    = c(fmt_ci(ci_gr$interior_low, ci_gr$interior_high, ci_gr$low_censurado, ci_gr$high_censurado, 1),
                         fmt_ci(ic_bb["low"], ic_bb["high"], d = 3)),
    ic_funresMech_interpolado = c(fmt_ci(ci_int$z_low, ci_int$z_high, ci_int$low_censurado, ci_int$high_censurado),
                                  fmt_ci(ic_bb["low"], ic_bb["high"], d = 3)),
    ic_okuyama       = c(sprintf("[%.1f, %.1f]", oks$ci_low, oks$ci_high),
                         sprintf("[%.3f, %.3f]", okb$ci_low, okb$ci_high)),
    a_funresMech = c(unname(par_min["a"]), unname(fit_bb$par[1])), a_okuyama = c(oks$a, okb$a),
    h_funresMech = c(unname(par_min["h"]), unname(fit_bb$par[2])), h_okuyama = c(oks$h, okb$h),
    k_o_rho_funresMech = c(unname(par_min["k"]), unname(fit_bb$par[4])), k_o_rho_okuyama = c(oks$rho_o_k, okb$rho_o_k),
    s_natural_funresMech = c(s_nat, NA), s_okuyama = c(oks$s, NA),
    AIC_funresMech = c(aic_sim, aic_bb), AIC_okuyama = c(oks$AIC, okb$AIC),
    delta_aic_bb_menos_sim_funresMech = aic_bb - aic_sim,
    delta_aic_bb_menos_sim_okuyama    = okb$AIC - oks$AIC,
    stringsAsFactors = FALSE
  )

  perfiles_largo[[lab]] <- rbind(
    data.frame(dataset = lab, modelo = "simulacion", z = prof$z, nll = prof$nll),
    data.frame(dataset = lab, modelo = "beta-binom", z = prof_bb$z, nll = prof_bb$nll)
  )
}

tabla <- do.call(rbind, resultados)
comp  <- do.call(rbind, comparacion)
perfs <- do.call(rbind, perfiles_largo)
dir.create("tablas", showWarnings = FALSE)
write.csv(tabla, file.path("tablas", "ensayo2_D2_D3_resultado_corr_OKU_v2.csv"), row.names = FALSE)
write.csv(comp,  file.path("tablas", "ensayo2_comparacion_vs_okuyama2026_v2.csv"), row.names = FALSE)
write.csv(perfs, file.path("tablas", "ensayo2_perfiles_corr_OKU.csv"), row.names = FALSE)

cat("\n=== Comparacion directa funresMech vs. Okuyama (2026) Tabla S1 ===\n")
print(comp[, c("dataset", "metodo", "z_hat_funresMech", "z_hat_okuyama", "ic_funresMech",
               "ic_okuyama", "AIC_funresMech", "AIC_okuyama")], row.names = FALSE)

cat("\nArchivos generados:\n",
    " - tablas/ensayo2_D2_D3_resultado_corr_OKU_v2.csv   (todo el detalle)\n",
    " - tablas/ensayo2_comparacion_vs_okuyama2026_v2.csv (lado a lado con Okuyama)\n",
    " - tablas/ensayo2_perfiles_corr_OKU.csv             (perfiles completos, para figuras)\n", sep = "")
cat(sprintf("\nVerificacion del criterio del IC contra Okuyama (bb, D2): %s\n",
            if (isTRUE(verificacion_ok)) "OK" else "FALLIDA - revisar"))
cat("\nNota: el perfil de simulacion tiene ruido Monte Carlo (n_sim=3000: desvio de la NLL ~0.35-0.45).\n",
    "Diferencias de z_hat o de limites del IC de un paso de grilla respecto de Okuyama son esperables.\n", sep = "")
