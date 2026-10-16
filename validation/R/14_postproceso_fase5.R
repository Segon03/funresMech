# ============================================================================
# 14_postproceso_fase5.R
#
# Correcciones sobre los resultados GUARDADOS de la Fase 5 (no reajusta nada):
#   (1) re-evalúa, con 4x n_sim, los puntos del perfil que tocaron el piso de
#       probabilidad (H16), en los parámetros ya estimados, y recalcula el IC
#       de z con interpolación lineal y con spline monótona (H17);
#   (2) recalibra los umbrales de los chequeos de confiabilidad (H18) con los
#       valores de los ajustes "buenos" (IC cubre z, sin falla de optimización,
#       error de z <= 1), usando la cobertura corregida.
# Escribe validation/tablas/fase5_post_*.csv
# Uso (desde funresMech_dev):  source("validation/R/14_postproceso_fase5.R")
# ============================================================================

suppressMessages(library(funresMech))
TAB <- "validation/tablas"; RES <- "validation/results/fase5"
lee <- function(x) read.csv(file.path(TAB, paste0("fase5_", x, ".csv")), stringsAsFactors = FALSE)
escribe <- function(x, nombre) write.csv(x, file.path(TAB, paste0("fase5_post_", nombre, ".csv")),
                                         row.names = FALSE)
clave <- function(x) paste(x$id, x$dato)
prop_ic <- function(x) {
  n <- length(x); p <- mean(x); z <- 1.96
  c(p = p, inf = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    sup = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), n = n)
}
ci_spline <- function(z, nll) {
  f <- stats::splinefun(z, nll, method = "monoH.FC")
  zz <- seq(min(z), max(z), by = 0.01)
  funresMech:::.profile_ci(zz, f(zz))
}
cubre <- function(ci, zt) (is.na(ci$z_low) || ci$z_low <= zt) && (is.na(ci$z_high) || zt <= ci$z_high)

A <- lee("A_simulacion"); E <- lee("E_chequeos")

# ---- (1) re-evaluación de los puntos con piso ---------------------------------
set.seed(20261016)
reevalua <- function(o, i, factor = 4, rondas = 3) {
  fit <- o$fit; p <- fit$profile[i, ]
  par <- c(a = p$a, h = p$h, z = p$z, k = p$k, s = p$s)
  n_sim <- fit$control$n_sim
  for (r in seq_len(rondas)) {
    n_sim <- n_sim * factor
    re <- lapply(seq_len(fit$control$n_reeval), function(j)
      funresMech:::.nll_sim(par, fit$model, fit$data, fit$T, n_sim, fit$control$cut, detail = TRUE))
    nll <- mean(vapply(re, function(x) x$nll, 0)); nf <- max(vapply(re, function(x) x$n_floor, 0L))
    if (nf == 0) break
  }
  c(nll = nll, n_floor = nf, n_sim = n_sim)
}

fa <- list.files(RES, "^A__.*\\.rds$", full.names = TRUE)
filas <- lapply(fa, function(f) {
  o <- readRDS(f); p <- o$fit$profile; zt <- as.numeric(sub(".*_z([0-9.]+)_k.*", "\\1", o$trabajo$id))
  nllr <- p$nll; nfr <- p$n_floor
  for (i in which(p$n_floor > 0)) { r <- reevalua(o, i); nllr[i] <- r[["nll"]]; nfr[i] <- r[["n_floor"]] }
  ord <- order(p$z); z <- p$z[ord]
  q2 <- p[ord, ][p$n_floor[ord] == 0, ]
  data.frame(id = o$trabajo$id, dato = o$trabajo$dato, z = zt, tuvo_piso = any(p$n_floor > 0),
             piso_persiste = any(nfr > 0),
             actual = cubre(funresMech:::.profile_ci(z, p$nll[ord]), zt),
             sin_piso = cubre(funresMech:::.profile_ci(q2$z, q2$nll), zt),
             reeval_lineal = cubre(funresMech:::.profile_ci(z, nllr[ord]), zt),
             reeval_spline = cubre(ci_spline(z, nllr[ord]), zt),
             zl_reeval = ci_spline(z, nllr[ord])$z_low, zh_reeval = ci_spline(z, nllr[ord])$z_high,
             zhat_reeval = z[which.min(nllr[ord])])
})
R <- do.call(rbind, filas)
escribe(R, "reevaluacion_por_ajuste")
res <- rbind(aggregate(cbind(actual, sin_piso, reeval_lineal, reeval_spline) ~ z, R, mean),
             data.frame(z = NA, t(colMeans(R[, c("actual", "sin_piso", "reeval_lineal", "reeval_spline")]))))
res$n_con_piso <- c(tapply(R$tuvo_piso, R$z, sum), sum(R$tuvo_piso))
res$piso_persiste <- c(tapply(R$piso_persiste, R$z, sum), sum(R$piso_persiste))
escribe(res, "cobertura")
print(res)
cat("\nCobertura global corregida (re-evaluación + spline):\n")
print(prop_ic(R$reeval_spline))

# ---- (2) recalibración de los umbrales ----------------------------------------
m <- match(clave(E), clave(A))
A$cubre_corr <- R$reeval_spline[match(clave(A), clave(R))]
A$bueno <- A$cubre_corr & A$deficit <= 1 & abs(A$z_hat - A$z) <= 1
E$bueno <- A$bueno[m]
auc <- function(x, y) { ok <- !is.na(x); x <- x[ok]; y <- y[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(x); n1 <- sum(y); n0 <- sum(!y); (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
th <- funresMech:::.rel_thresholds()
cal <- do.call(rbind, lapply(split(E, E$check), function(q) {
  data.frame(check = q$check[1], n_buenos = sum(q$bueno), n_malos = sum(!q$bueno),
             alerta_actual = mean(q$status != "ok"),
             alerta_en_buenos = mean(q$status[q$bueno] != "ok"),
             alerta_en_malos = mean(q$status[!q$bueno] != "ok"),
             auc_malo = auc(q$value, !q$bueno),
             p50_buenos = unname(quantile(q$value[q$bueno], 0.5, na.rm = TRUE)),
             p90_buenos = unname(quantile(q$value[q$bueno], 0.9, na.rm = TRUE)),
             p90_malos = unname(quantile(q$value[!q$bueno], 0.9, na.rm = TRUE)))
}))
escribe(cal, "calibracion_chequeos")
print(cal, digits = 3); cat("\nUmbrales actuales:\n"); str(th)
