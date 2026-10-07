# ============================================================
# R/04i_ensayo2_rcpp_reensamblar.R
#
# Arma los resultados del Ensayo 2 a partir de los perfiles de
# 04h_ensayo2_rcpp_log.R (motor Rcpp) y los compara con la Tabla S1
# de Okuyama (2026). Rapido (segundos); se puede correr cuantas veces
# se quiera, incluso con la corrida de 04h a medio terminar (usa los
# puntos que ya existan y avisa cuales faltan).
#
# Para cada dataset (D2, D3) y cada escala de busqueda ("log" =
# principal; "okuyama" = ablacion, si se corrio):
#   - z_hat = argmin del perfil (como Okuyama, Sec. 2.3)
#   - IC95% de z con el umbral correcto NLL_min + qchisq(.95,1)/2,
#     sobre la grilla (como Okuyama) y por interpolacion
#   - parametros (a, h, k, s) en z_hat, AIC = 2*NLL + 2*5
#   - lo mismo usando la NLL RE-EVALUADA (robustez: sin el sesgo del
#     minimo de DEoptim)
#   - comparacion contra Okuyama con tolerancias explicitas (abajo)
#   - diagnostico de suavidad del perfil (saltos entre z vecinos)
#   - Delta AIC beta-binomial - simulacion (conclusion de Okuyama)
#
# Salidas:
#   tablas/ensayo2_rcpp_perfiles.csv       todos los puntos (con parametros)
#   tablas/ensayo2_rcpp_resultado.csv      una fila por dataset x escala x NLL
#   tablas/ensayo2_rcpp_vs_okuyama.csv     lado a lado con la Tabla S1 + veredicto
#   figuras/fig_ensayo2_rcpp_perfiles.png  perfiles (Rcpp log / Rcpp okuyama /
#                                           motor R anterior) con IC y Okuyama
#
# REVISION (01-oct-2026): veredicto revisado + extension de D2
#   Se agrega un veredicto REVISADO (columna veredicto_revisado) SIN
#   tocar el original (columna veredicto, mismas reglas de antes),
#   porque el perfil de D2 es casi plano entre z ~1.6 y 6 y las reglas
#   originales lo leen mal:
#    - z_hat: con un perfil plano el argmin cae donde lo deja el ruido
#      Monte Carlo. Revisado: el z_hat de Okuyama es equivalente si
#      NLL(z_hat Okuyama) - NLL_min <= TOL_DNLL, con
#      TOL_DNLL = 2*sqrt(2)*SD de UNA evaluacion de la NLL (diferencia
#      entre dos puntos ruidosos, ~95%); la SD se estima del propio
#      perfil: mediana(nll_reeval_se * sqrt(n_reeval)).
#    - IC: si el perfil no cruza el umbral dentro de la grilla, ese
#      lado se informa "no acotado" y no como falla. Se informa ademas
#      NLL - NLL_min en los limites de Okuyama: si es < 1.92, nuestro
#      perfil deja ese limite DENTRO del IC.
#   Si existe results/ensayo2_rcpp/log_ext/ (R/04j), se agrega una
#   seccion con el perfil extendido de D2 (z hasta 20):
#     tablas/ensayo2_rcpp_extension_perfil.csv
#     tablas/ensayo2_rcpp_extension_resumen.csv
#     figuras/fig_ensayo2_rcpp_extension.png
#   La comparacion 1 a 1 con Okuyama sigue usando SOLO la grilla
#   0.5-6.0 (la suya); la extension no entra en esas tablas.
# ============================================================

source(file.path("fun", "profile_utils.R"))

MODO_PRUEBA <- FALSE   # TRUE para leer la corrida de humo de 04h
dir_res <- file.path("results", if (MODO_PRUEBA) "ensayo2_rcpp_PRUEBA" else "ensayo2_rcpp")
labs    <- c("D2", "D3")
T_exp   <- 24
Z_GRID  <- round(seq(0.5, 6.0, by = 0.1), 2)

# Tolerancias para declarar "reproduce a Okuyama". El perfil simulado
# tiene ruido Monte Carlo (NLL con n_sim = 10000: desvio ~0.2-0.3) y
# Okuyama informa z con 1 decimal y su IC "puede ser hasta 0.1 mas
# ancho de cada lado" (Sec. 2.3). Se acepta:
TOL_Z_HAT <- 0.3    # |z_hat - z_hat Okuyama|
TOL_IC    <- 0.4    # |limite IC - limite IC Okuyama|, cada lado
TOL_AIC   <- 3      # |AIC - AIC Okuyama| (~1.5 unidades de NLL)

okuyama <- data.frame(   # Tabla S1 (jen70148-sup-0009-tables1.xlsx), ML: simulation
  dataset = c("D2", "D3"), a = c(0.001, 0.012), h = c(0.606, 0.626),
  z_hat = c(3.1, 2.5), ci_low = c(1.1, 1.4), ci_high = c(4.6, 4.8),
  k = c(0.037, 0.902), s = c(0.682, 0.542), AIC = c(336.589, 272.890),
  AIC_bb = c(346.916, 273.507), stringsAsFactors = FALSE)

# ---- Leer perfiles ------------------------------------------------------------
escalas <- intersect(c("log", "okuyama"), list.dirs(dir_res, full.names = FALSE, recursive = FALSE))
if (!length(escalas)) stop("No hay resultados en ", dir_res, ". Correr 04h primero.")

leer_perfil <- function(f, esc) {
  r <- readRDS(f)
  data.frame(
    escala = esc, dataset = r$dataset, z = round(r$z, 2),
    nll = r$nll, nll_reeval = r$nll_reeval, nll_reeval_se = r$nll_reeval_se,
    a = r$par[["a"]], h = r$par[["h"]], k = r$par[["k"]], s = r$par[["s"]],
    en_cota = paste(names(r$en_cota)[r$en_cota], collapse = ","),
    minutos = r$segundos / 60, nsim = r$config$nsim, itermax = r$config$itermax,
    n_reeval = if (is.null(r$config$n_reeval)) NA_real_ else r$config$n_reeval,
    stringsAsFactors = FALSE)
}
filas <- list()
for (esc in escalas) for (f in list.files(file.path(dir_res, esc), pattern = "^perfil_.*\\.rds$", full.names = TRUE))
  filas[[length(filas) + 1]] <- leer_perfil(f, esc)
perf <- do.call(rbind, filas)
perf <- perf[order(perf$escala, perf$dataset, perf$z), ]
dir.create("tablas", showWarnings = FALSE); dir.create("figuras", showWarnings = FALSE)
write.csv(perf, file.path("tablas", "ensayo2_rcpp_perfiles.csv"), row.names = FALSE)

# Perfil del motor R anterior (04d/04f + 04g), solo para la figura
perf_viejo <- if (file.exists(file.path("tablas", "ensayo2_perfiles_corr_OKU.csv"))) {
  p <- read.csv(file.path("tablas", "ensayo2_perfiles_corr_OKU.csv"))
  p[p$modelo == "simulacion", ]
} else NULL

fmt_ci <- function(lo, hi, lo_c = FALSE, hi_c = FALSE, d = 1) {
  f <- paste0("%.", d, "f")
  sprintf("[%s, %s]", if (isTRUE(lo_c)) paste0("<=", sprintf(f, lo)) else sprintf(f, lo),
          if (isTRUE(hi_c)) paste0(">=", sprintf(f, hi)) else sprintf(f, hi))
}

# Ruido Monte Carlo de UNA evaluacion de la NLL, estimado del propio perfil
# (mediana: robusta a los pocos puntos con re-evaluaciones inestables)
sd_mc_perfil <- function(p) {
  sd1 <- p$nll_reeval_se * sqrt(p$n_reeval)
  median(sd1[is.finite(sd1)], na.rm = TRUE)
}
# Estado de un lado del IC frente a Okuyama: ok / difiere / no_acotado
estado_lado <- function(censurado, lim, lim_ok) {
  if (isTRUE(censurado)) "no_acotado" else if (abs(lim - lim_ok) <= TOL_IC + 1e-9) "ok" else "difiere"
}
ETIQ_LADO <- c(ic_low = "IC inf.", ic_high = "IC sup.")

# ---- Resultados por dataset x escala x tipo de NLL ----------------------------
res <- list(); comp <- list()
for (esc in escalas) for (lab in labs) {
  p <- perf[perf$escala == esc & perf$dataset == lab, ]
  if (nrow(p) < 3) { cat(sprintf("%s/%s: menos de 3 puntos, se omite.\n", esc, lab)); next }
  faltan <- setdiff(if (MODO_PRUEBA) numeric(0) else Z_GRID, p$z)
  ok <- okuyama[okuyama$dataset == lab, ]
  bb_f <- file.path("results", "ensayo2_corr_OKU", sprintf("bb_%s_corrOKU.rds", lab))
  aic_bb <- if (file.exists(bb_f)) 2 * readRDS(bb_f)$value + 2 * 4 else NA_real_

  cat(sprintf("\n========== %s | escala de busqueda: %s ==========\n", lab, esc))
  cat(sprintf("Puntos del perfil: %d%s\n", nrow(p),
              if (length(faltan)) sprintf(" (FALTAN %d: %s)", length(faltan),
                                          paste(head(faltan, 12), collapse = ", ")) else " (grilla completa)"))
  cotas <- p$en_cota[nzchar(p$en_cota)]
  if (length(cotas)) cat(sprintf("  Aviso: %d puntos con algun parametro en la cota de busqueda (%s)\n",
                                 length(cotas), paste(unique(unlist(strsplit(cotas, ","))), collapse = ", ")))

  for (tipo in c("nll", "nll_reeval")) {
    v <- p[[tipo]]
    i_min <- which.min(v)
    ci  <- ci_from_profile_okuyama(p$z, v)
    cig <- ci_grilla_okuyama(p$z, v)
    aic <- 2 * v[i_min] + 2 * 5
    # suavidad: segunda diferencia (curvatura local) muy negativa = "diente"
    saltos <- if (length(v) >= 3) max(abs(diff(v, differences = 2))) else NA
    etiqueta <- if (tipo == "nll") "NLL de DEoptim (como Okuyama)" else "NLL re-evaluada (robustez)"
    cat(sprintf("[%s]\n  z_hat = %.1f (Okuyama %.1f) | IC grilla %s | IC interpolado %s | Okuyama [%.1f, %.1f]\n",
                etiqueta, p$z[i_min], ok$z_hat,
                fmt_ci(cig$interior_low, cig$interior_high, cig$low_censurado, cig$high_censurado),
                fmt_ci(if (ci$low_censurado) min(p$z) else ci$z_low,
                       if (ci$high_censurado) max(p$z) else ci$z_high,
                       ci$low_censurado, ci$high_censurado, 2),
                ok$ci_low, ok$ci_high))
    cat(sprintf("  NLL_min = %.3f | AIC = %.3f (Okuyama %.3f) | a=%.4g h=%.3f k=%.3g s=%.3f (Okuyama a=%.3f h=%.3f k=%.3f s=%.3f)\n",
                v[i_min], aic, ok$AIC, p$a[i_min], p$h[i_min], p$k[i_min], p$s[i_min],
                ok$a, ok$h, ok$k, ok$s))
    cat(sprintf("  Delta AIC (beta-binom - simulacion) = %.2f (Okuyama %.2f) | max |2a diferencia| del perfil = %.2f\n",
                aic_bb - aic, ok$AIC_bb - ok$AIC, saltos))

    # IC para comparar: grilla interior (lo que informa Okuyama)
    lo <- cig$interior_low; hi <- cig$interior_high
    chk <- c(z_hat = abs(p$z[i_min] - ok$z_hat) <= TOL_Z_HAT + 1e-9,
             ic_low = !cig$low_censurado && abs(lo - ok$ci_low) <= TOL_IC + 1e-9,
             ic_high = !cig$high_censurado && abs(hi - ok$ci_high) <= TOL_IC + 1e-9,
             aic = abs(aic - ok$AIC) <= TOL_AIC,
             conclusion_z_mayor_1 = (lo > 1) == (ok$ci_low > 1),
             grilla_completa = length(faltan) == 0)
    veredicto <- if (all(chk)) "REPRODUCE" else paste("NO:", paste(names(chk)[!chk], collapse = ","))
    cat(sprintf("  Veredicto (tol. z %.1f, IC %.1f, AIC %.0f): %s\n", TOL_Z_HAT, TOL_IC, TOL_AIC, veredicto))

    # ---- Veredicto REVISADO (ver encabezado) ----
    tol_dnll       <- 2 * sqrt(2) * sd_mc_perfil(p)
    dnll_zhat_ok   <- nll_at_z(p$z, v, ok$z_hat)   - v[i_min]
    dnll_iclow_ok  <- nll_at_z(p$z, v, ok$ci_low)  - v[i_min]
    dnll_ichigh_ok <- nll_at_z(p$z, v, ok$ci_high) - v[i_min]
    est_low  <- estado_lado(cig$low_censurado,  lo, ok$ci_low)
    est_high <- estado_lado(cig$high_censurado, hi, ok$ci_high)
    chk2 <- c(z_hat_equivalente = isTRUE(dnll_zhat_ok <= tol_dnll),
              ic_low = est_low != "difiere", ic_high = est_high != "difiere",
              aic = chk[["aic"]], conclusion_z_mayor_1 = chk[["conclusion_z_mayor_1"]],
              grilla_completa = chk[["grilla_completa"]])
    no_acot <- c(ic_low = est_low, ic_high = est_high)
    no_acot <- names(no_acot)[no_acot == "no_acotado"]
    veredicto_rev <- if (!all(chk2)) {
      paste("NO:", paste(names(chk2)[!chk2], collapse = ","))
    } else if (length(no_acot)) {
      paste0("REPRODUCE salvo ", paste(ETIQ_LADO[no_acot], collapse = " e "), " no acotado en la grilla")
    } else "REPRODUCE"
    cat(sprintf(paste0("  Revisado: NLL-NLL_min en z_hat de Okuyama (%.1f) = %.2f (tol. ruido MC %.2f) |",
                       " en sus limites de IC: z=%.1f -> %.2f, z=%.1f -> %.2f (umbral %.2f)\n",
                       "  Veredicto revisado: %s\n"),
                ok$z_hat, dnll_zhat_ok, tol_dnll, ok$ci_low, dnll_iclow_ok,
                ok$ci_high, dnll_ichigh_ok, UMBRAL_IC95, veredicto_rev))

    res[[length(res) + 1]] <- data.frame(
      dataset = lab, escala = esc, tipo_nll = tipo, n_puntos = nrow(p), faltan = length(faltan),
      z_hat = p$z[i_min], nll_min = v[i_min], aic = aic,
      ic_grilla_low = lo, ic_grilla_high = hi,
      ic_grilla_ext_low = cig$exterior_low, ic_grilla_ext_high = cig$exterior_high,
      ic_interp_low = ci$z_low, ic_interp_high = ci$z_high,
      ic_low_censurado = cig$low_censurado, ic_high_censurado = cig$high_censurado,
      a = p$a[i_min], h = p$h[i_min], k = p$k[i_min], s = p$s[i_min],
      aic_bb = aic_bb, delta_aic_bb_menos_sim = aic_bb - aic,
      max_segunda_diferencia = saltos, veredicto = veredicto,
      veredicto_revisado = veredicto_rev, tol_dnll = tol_dnll,
      dnll_en_zhat_okuyama = dnll_zhat_ok, dnll_en_ic_low_okuyama = dnll_iclow_ok,
      dnll_en_ic_high_okuyama = dnll_ichigh_ok, estado_ic_low = est_low, estado_ic_high = est_high,
      stringsAsFactors = FALSE)
    comp[[length(comp) + 1]] <- data.frame(
      dataset = lab, escala = esc, tipo_nll = tipo,
      z_hat = p$z[i_min], z_hat_okuyama = ok$z_hat,
      ic = fmt_ci(lo, hi, cig$low_censurado, cig$high_censurado),
      ic_okuyama = sprintf("[%.1f, %.1f]", ok$ci_low, ok$ci_high),
      aic = round(aic, 3), aic_okuyama = ok$AIC,
      delta_aic_bb_sim = round(aic_bb - aic, 2), delta_aic_bb_sim_okuyama = round(ok$AIC_bb - ok$AIC, 2),
      a = signif(p$a[i_min], 3), a_okuyama = ok$a, h = round(p$h[i_min], 3), h_okuyama = ok$h,
      k = signif(p$k[i_min], 3), k_okuyama = ok$k, s = round(p$s[i_min], 3), s_okuyama = ok$s,
      veredicto = veredicto,
      dnll_en_zhat_okuyama = round(dnll_zhat_ok, 2), tol_dnll = round(tol_dnll, 2),
      dnll_en_ic_high_okuyama = round(dnll_ichigh_ok, 2),
      veredicto_revisado = veredicto_rev, stringsAsFactors = FALSE)
  }
}
tabla <- do.call(rbind, res); tcomp <- do.call(rbind, comp)
write.csv(tabla, file.path("tablas", "ensayo2_rcpp_resultado.csv"), row.names = FALSE)
write.csv(tcomp, file.path("tablas", "ensayo2_rcpp_vs_okuyama.csv"), row.names = FALSE)

# ---- Extension del perfil (R/04j), si existe ----------------------------------
dir_ext <- file.path(dir_res, "log_ext")
f_ext <- if (dir.exists(dir_ext)) list.files(dir_ext, pattern = "^perfil_.*\\.rds$", full.names = TRUE) else character(0)
if (length(f_ext) && "log" %in% escalas) {
  ext <- do.call(rbind, lapply(f_ext, leer_perfil, esc = "log"))
  res_ext <- list(); perf_ext <- list()
  for (lab in unique(ext$dataset)) {
    base <- perf[perf$escala == "log" & perf$dataset == lab, ]
    if (!nrow(base)) next
    e <- ext[ext$dataset == lab & !(ext$z %in% base$z), ]
    comb <- rbind(base, e); comb <- comb[order(comb$z), ]
    comb$origen <- ifelse(comb$z %in% e$z, "04j (extension)", "04h (grilla 0.5-6.0)")
    ok <- okuyama[okuyama$dataset == lab, ]
    z_fin <- max(base$z)

    cat(sprintf("\n========== %s | escala log | PERFIL EXTENDIDO (R/04j): z de %.1f a %.1f ==========\n",
                lab, min(comb$z), max(comb$z)))
    cat(sprintf("Puntos agregados: %s\n", paste(sprintf("%.1f", sort(e$z)), collapse = ", ")))
    ec <- e[nzchar(e$en_cota), ]
    if (nrow(ec)) cat(sprintf("  AVISO: parametro en la cota de busqueda en z = %s (%s): esos puntos pueden\n  estar subajustados (NLL inflada) y el cruce del umbral ahi no seria confiable.\n",
                              paste(sprintf("%.1f", ec$z), collapse = ", "), paste(unique(ec$en_cota), collapse = "; ")))
    cat("  z      a          h      k        s      dNLL(DEoptim)  dNLL(reeval)\n")
    for (i in which(comb$z >= z_fin))
      cat(sprintf("  %-5.1f  %-9.3g  %.3f  %-7.3g  %.3f  %6.2f         %6.2f\n", comb$z[i], comb$a[i], comb$h[i],
                  comb$k[i], comb$s[i], comb$nll[i] - min(comb$nll), comb$nll_reeval[i] - min(comb$nll_reeval)))

    comb$dnll <- comb$nll - min(comb$nll)
    comb$dnll_reeval <- comb$nll_reeval - min(comb$nll_reeval)
    perf_ext[[lab]] <- comb

    for (tipo in c("nll", "nll_reeval")) {
      v <- comb[[tipo]]; i_min <- which.min(v); d <- v - v[i_min]
      ci  <- ci_from_profile_okuyama(comb$z, v)
      cig <- ci_grilla_okuyama(comb$z, v)
      d_ext <- d[comb$z > z_fin]
      lectura <- if (!ci$high_censurado) {
        sprintf("el perfil cruza +%.2f en z ~ %.1f (entre %.1f y %.1f de la grilla): limite superior acotado",
                UMBRAL_IC95, ci$z_high, cig$interior_high, cig$exterior_high)
      } else if (max(d_ext) <= 1) {
        sprintf("perfil plano hasta z = %.0f (dNLL max %.2f << %.2f): z no identificable por arriba con estos datos",
                max(comb$z), max(d_ext), UMBRAL_IC95)
      } else {
        sprintf("el perfil sube (dNLL max %.2f) pero no llega a +%.2f hasta z = %.0f: limite superior > %.0f",
                max(d_ext), UMBRAL_IC95, max(comb$z), max(comb$z))
      }
      etiqueta <- if (tipo == "nll") "NLL de DEoptim" else "NLL re-evaluada"
      cat(sprintf("[%s]\n  z_hat = %.1f%s | IC interpolado [%s, %s] | dNLL en z = %.1f (Okuyama IC sup.) = %.2f\n  -> %s\n",
                  etiqueta, comb$z[i_min],
                  if (comb$z[i_min] > z_fin) "  (OJO: el minimo cae en la extension)" else "",
                  if (ci$low_censurado) sprintf("<=%.1f", min(comb$z)) else sprintf("%.2f", ci$z_low),
                  if (ci$high_censurado) sprintf(">=%.1f", max(comb$z)) else sprintf("%.2f", ci$z_high),
                  ok$ci_high, nll_at_z(comb$z, v, ok$ci_high) - v[i_min], lectura))

      res_ext[[length(res_ext) + 1]] <- data.frame(
        dataset = lab, escala = "log", tipo_nll = tipo, z_max_explorado = max(comb$z),
        z_hat = comb$z[i_min], nll_min = v[i_min], aic = 2 * v[i_min] + 2 * 5,
        ic_interp_low = ci$z_low, ic_interp_high = ci$z_high,
        ic_low_censurado = ci$low_censurado, ic_high_censurado = ci$high_censurado,
        ic_grilla_interior_high = cig$interior_high, ic_grilla_exterior_high = cig$exterior_high,
        dnll_max_extension = max(d_ext), dnll_en_z_max = d[length(d)],
        dnll_en_ic_high_okuyama = nll_at_z(comb$z, v, ok$ci_high) - v[i_min],
        tol_dnll = 2 * sqrt(2) * sd_mc_perfil(comb),
        puntos_en_cota = paste(sprintf("%.1f", ec$z), collapse = ";"),
        lectura = lectura, stringsAsFactors = FALSE)
    }
  }
  if (length(perf_ext)) {
    write.csv(do.call(rbind, perf_ext), file.path("tablas", "ensayo2_rcpp_extension_perfil.csv"), row.names = FALSE)
    write.csv(do.call(rbind, res_ext), file.path("tablas", "ensayo2_rcpp_extension_resumen.csv"), row.names = FALSE)

    png(file.path("figuras", "fig_ensayo2_rcpp_extension.png"), width = 1100 * length(perf_ext), height = 900, res = 170)
    op <- par(mfrow = c(1, length(perf_ext)), mar = c(4.2, 4.5, 2.5, 1))
    for (lab in names(perf_ext)) {
      cb <- perf_ext[[lab]]; ok <- okuyama[okuyama$dataset == lab, ]
      ymax <- max(3, min(15, max(c(cb$dnll, cb$dnll_reeval), na.rm = TRUE)))
      plot(NA, xlim = range(cb$z), ylim = c(0, ymax), xlab = "z", ylab = expression(NLL - NLL[min]),
           main = sprintf("%s, escala log: perfil extendido (Okuyama IC [%.1f, %.1f])", lab, ok$ci_low, ok$ci_high),
           cex.main = 0.9)
      rect(ok$ci_low, -1, ok$ci_high, ymax + 1, col = "#00000010", border = NA)
      abline(v = max(Z_GRID), col = "grey50", lty = 3)
      abline(h = UMBRAL_IC95, col = "firebrick", lty = 2)
      lines(cb$z, cb$dnll, col = "#1f5fa8", lwd = 2)
      lines(cb$z, cb$dnll_reeval, col = "#1f5fa8", lwd = 1, lty = 2)
      es_ext <- cb$origen != "04h (grilla 0.5-6.0)"
      points(cb$z[!es_ext], cb$dnll[!es_ext], pch = 16, cex = 0.5, col = "#1f5fa8")
      points(cb$z[es_ext], cb$dnll[es_ext], pch = 21, cex = 1.1, bg = "#d9822b", col = "black")
      if (any(nzchar(cb$en_cota))) points(cb$z[nzchar(cb$en_cota)], cb$dnll[nzchar(cb$en_cota)], pch = 4, cex = 1.6, col = "red")
      legend("topright", c("NLL DEoptim", "NLL re-evaluada", "puntos 04j (extension)", "parametro en cota",
                           "umbral IC95 (+1.92)", "fin grilla Okuyama (z = 6)", "IC Okuyama"),
             col = c("#1f5fa8", "#1f5fa8", "black", "red", "firebrick", "grey50", "#00000030"),
             pt.bg = c(NA, NA, "#d9822b", NA, NA, NA, NA), lty = c(1, 2, NA, NA, 2, 3, NA),
             pch = c(NA, NA, 21, 4, NA, NA, 15), lwd = c(2, 1, NA, NA, 1, 1, NA),
             pt.cex = c(1, 1, 1.1, 1.3, 1, 1, 2), bty = "o", bg = "white", box.col = "grey80", cex = 0.72)
    }
    par(op); dev.off()
  }
}

# ---- Figura -------------------------------------------------------------------
png(file.path("figuras", "fig_ensayo2_rcpp_perfiles.png"), width = 2000, height = 900, res = 170)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.5, 1))
for (lab in labs) {
  series <- list()
  for (esc in escalas) {
    p <- perf[perf$escala == esc & perf$dataset == lab, ]
    if (nrow(p)) series[[paste("Rcpp", esc)]] <- data.frame(z = p$z, d = p$nll - min(p$nll))
  }
  if (!is.null(perf_viejo)) {
    pv <- perf_viejo[perf_viejo$dataset == lab, ]
    if (nrow(pv)) series[["Motor R v1.0.4 (04d/04f)"]] <- data.frame(z = pv$z, d = pv$nll - min(pv$nll))
  }
  if (!length(series)) next
  cols <- c("Rcpp log" = "#1f5fa8", "Rcpp okuyama" = "#d9822b", "Motor R v1.0.4 (04d/04f)" = "grey55")
  ltys <- c("Rcpp log" = 1, "Rcpp okuyama" = 2, "Motor R v1.0.4 (04d/04f)" = 3)
  ymax <- min(15, max(unlist(lapply(series, function(s) s$d)), na.rm = TRUE))
  plot(NA, xlim = c(0.5, 6), ylim = c(0, ymax), xlab = "z", ylab = expression(NLL - NLL[min]),
       main = sprintf("%s (Okuyama: z = %.1f, IC [%.1f, %.1f])", lab,
                      okuyama$z_hat[okuyama$dataset == lab], okuyama$ci_low[okuyama$dataset == lab],
                      okuyama$ci_high[okuyama$dataset == lab]), cex.main = 0.95)
  ok <- okuyama[okuyama$dataset == lab, ]
  rect(ok$ci_low, -1, ok$ci_high, ymax + 1, col = "#00000010", border = NA)
  abline(v = ok$z_hat, col = "grey30", lty = 2)
  abline(h = UMBRAL_IC95, col = "firebrick", lty = 2)
  for (nm in names(series)) lines(series[[nm]]$z, series[[nm]]$d, col = cols[nm], lty = ltys[nm],
                                  lwd = 2, type = "o", pch = 16, cex = 0.5)
  legend("top", c(names(series), "umbral IC95 (+1.92)", "Okuyama z_hat / IC"),
         col = c(cols[names(series)], "firebrick", "grey30"), lty = c(ltys[names(series)], 2, 2),
         lwd = c(rep(2, length(series)), 1, 1), bty = "n", cex = 0.75)
}
par(op); dev.off()

cat("\n=== Resumen vs. Okuyama (2026) Tabla S1 ===\n")
print(tcomp[, c("dataset", "escala", "tipo_nll", "z_hat", "z_hat_okuyama", "ic", "ic_okuyama",
                "aic", "aic_okuyama", "veredicto", "veredicto_revisado")], row.names = FALSE)
cat("\nArchivos: tablas/ensayo2_rcpp_perfiles.csv, tablas/ensayo2_rcpp_resultado.csv,",
    "tablas/ensayo2_rcpp_vs_okuyama.csv, figuras/fig_ensayo2_rcpp_perfiles.png\n")
if (exists("perf_ext") && length(perf_ext)) {
  cat("\n=== Perfil extendido (R/04j) ===\n")
  print(do.call(rbind, res_ext)[, c("dataset", "tipo_nll", "z_max_explorado", "z_hat", "ic_interp_low",
                                    "ic_interp_high", "ic_high_censurado", "dnll_max_extension", "lectura")],
        row.names = FALSE)
  cat("Archivos: tablas/ensayo2_rcpp_extension_perfil.csv, tablas/ensayo2_rcpp_extension_resumen.csv,",
      "figuras/fig_ensayo2_rcpp_extension.png\n")
}
