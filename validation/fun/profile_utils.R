# ============================================================
# fun/profile_utils.R
#
# IC95% de z a partir de un perfil de verosimilitud, y AIC_full
# vs. AIC_restricted (z=1).
#
# ------------------------------------------------------------
# CORRECCION (25-sep-2026) - UMBRAL DEL IC
# ------------------------------------------------------------
# La version anterior replicaba tal cual el calculo de R/server.R
# de funresMech v1.0.4 (lineas 225-226):
#     threshold <- min_nll + qchisq(0.95, 1)          # = +3.84
# Eso esta MAL: el criterio de razon de verosimilitudes (Hilborn
# & Mangel 1997; Bolker 2008 - las mismas referencias que cita
# Okuyama 2026, Sec. 2.3) es
#     2 * (NLL(z) - NLL_min) <= qchisq(0.95, 1)
#     <=>  NLL(z) <= NLL_min + qchisq(0.95, 1) / 2      # = +1.92
# Sin el "/2" el intervalo resultante es en realidad un IC de
# ~99.5%, no de 95%.
#
# VERIFICACION CONTRA OKUYAMA (hecha antes de este cambio): se
# recalculo el perfil beta-binomial de D2 (T=24, modelo cerrado,
# sin ruido Monte Carlo) y se busco el IC con ambos umbrales:
#     umbral +1.92  ->  [1.310, 2.552]   <- Okuyama Tabla S1: [1.310, 2.552]
#     umbral +3.84  ->  [1.098, 2.912]
# Coincidencia exacta a 3 decimales con +1.92. El script
# 04g_ensayo2_corr_OKU_reensamblar_v2.R repite esta verificacion
# automaticamente dentro de R cada vez que se corre.
#
# El argumento use_standard_factor se conserva por compatibilidad,
# pero ahora el default es TRUE (criterio correcto). Pasar FALSE
# solo para documentar/reproducir lo que reportaba la app v1.0.4.
#
# NOVEDADES:
#  - ci_from_profile_okuyama() devuelve ademas low_censurado /
#    high_censurado: TRUE si el perfil nunca cruza el umbral de
#    ese lado dentro de la grilla (antes solo devolvia NA sin
#    explicar por que). En ese caso el limite real es <= min(grilla)
#    o >= max(grilla).
#  - ci_grilla_okuyama(): IC sobre la grilla discreta, SIN
#    interpolar, como describe Okuyama (2026, Sec. 2.3).
# ============================================================

UMBRAL_IC95 <- qchisq(0.95, 1) / 2   # = 1.920729

#' IC95% de z por interpolacion lineal entre puntos del perfil
#' (mismo patron de busqueda que R/server.R, lineas ~225-255, pero
#' con el umbral corregido).
#'
#' @param z_grid vector de valores de z evaluados (ordenado)
#' @param nll_vals NLL minimizado en cada z_grid (misma longitud)
#' @param use_standard_factor TRUE (default): umbral correcto
#'   NLL_min + qchisq(.95,1)/2. FALSE: umbral de la app v1.0.4
#'   (sin /2), solo para comparacion.
ci_from_profile_okuyama <- function(z_grid, nll_vals, use_standard_factor = TRUE) {
  stopifnot(length(z_grid) == length(nll_vals), !is.unsorted(z_grid))
  min_nll   <- min(nll_vals)
  idx_min   <- which.min(nll_vals)
  divisor   <- if (use_standard_factor) 2 else 1
  threshold <- min_nll + qchisq(0.95, 1) / divisor

  # Lado izquierdo: desde el minimo hacia z chicos, primer punto que
  # cruza el umbral; se interpola entre ese punto y el siguiente.
  z_low <- NA_real_
  for (i in seq(idx_min, 1, by = -1)) {
    if (nll_vals[i] >= threshold) {
      z1 <- z_grid[i];     n1 <- nll_vals[i]
      z2 <- z_grid[i + 1]; n2 <- nll_vals[i + 1]
      z_low <- z1 + (z2 - z1) * (threshold - n1) / (n2 - n1)
      break
    }
  }

  # Lado derecho: idem hacia z grandes.
  z_high <- NA_real_
  for (i in seq(idx_min, length(z_grid), by = 1)) {
    if (nll_vals[i] >= threshold) {
      z1 <- z_grid[i - 1]; n1 <- nll_vals[i - 1]
      z2 <- z_grid[i];     n2 <- nll_vals[i]
      z_high <- z1 + (z2 - z1) * (threshold - n1) / (n2 - n1)
      break
    }
  }

  list(z_hat = z_grid[idx_min], z_low = z_low, z_high = z_high,
       low_censurado  = is.na(z_low),    # limite real <= min(z_grid)
       high_censurado = is.na(z_high),   # limite real >= max(z_grid)
       threshold = threshold, min_nll = min_nll)
}

#' IC95% de z sobre la grilla discreta, sin interpolar (Okuyama 2026,
#' Sec. 2.3: "the smallest range of z values satisfying the likelihood
#' ratio criterion"). Se toma el tramo CONTIGUO de puntos alrededor del
#' minimo que cumplen NLL <= NLL_min + qchisq(.95,1)/2.
#'
#' La frase de Okuyama admite dos lecturas; se devuelven ambas:
#'  - interior: primer/ultimo punto de la grilla que SI cumple el
#'    criterio.
#'  - exterior: primer punto que YA NO cumple, a cada lado (coincide
#'    con su aclaracion de que el IC "may be up to 0.1 units wider on
#'    each side than the exact interval").
#' El IC exacto (interpolado) siempre queda entre ambas.
ci_grilla_okuyama <- function(z_grid, nll_vals) {
  stopifnot(length(z_grid) == length(nll_vals), !is.unsorted(z_grid))
  threshold <- min(nll_vals) + UMBRAL_IC95
  idx_min <- which.min(nll_vals)
  lo <- idx_min
  while (lo > 1 && nll_vals[lo - 1] <= threshold) lo <- lo - 1
  hi <- idx_min
  while (hi < length(z_grid) && nll_vals[hi + 1] <= threshold) hi <- hi + 1
  list(
    z_hat = z_grid[idx_min],
    interior_low  = z_grid[lo],
    interior_high = z_grid[hi],
    exterior_low  = if (lo > 1) z_grid[lo - 1] else NA_real_,
    exterior_high = if (hi < length(z_grid)) z_grid[hi + 1] else NA_real_,
    low_censurado  = lo == 1,
    high_censurado = hi == length(z_grid),
    threshold = threshold
  )
}

#' Interpola linealmente el NLL del perfil en un valor z0 dado
#' (replica el patron de server.R lineas ~186-207, usado ahi
#' para z0 = 1)
nll_at_z <- function(z_grid, nll_vals, z0) {
  if (z0 < min(z_grid) || z0 > max(z_grid)) return(NA_real_)
  idx <- which.min(abs(z_grid - z0))
  if (abs(z_grid[idx] - z0) < 0.01) return(nll_vals[idx])
  idx_low  <- max(which(z_grid <= z0))
  idx_high <- min(which(z_grid >= z0))
  z_low <- z_grid[idx_low];   n_low  <- nll_vals[idx_low]
  z_high <- z_grid[idx_high]; n_high <- nll_vals[idx_high]
  n_low + (n_high - n_low) * (z0 - z_low) / (z_high - z_low)
}

#' AIC_full (z libre, 5 parametros) vs. AIC_restricted (z=1, 4
#' parametros), tal como los calcula server.R. delta_aic > 0
#' favorece el modelo con z libre (tipo III).
aic_full_vs_restricted <- function(nll_full, z_grid, nll_profile) {
  aic_full <- 2 * nll_full + 2 * 5
  nll_1 <- nll_at_z(z_grid, nll_profile, 1)
  if (is.na(nll_1)) return(list(aic_full = aic_full, aic_restricted = NA, delta_aic = NA))
  aic_restricted <- 2 * nll_1 + 2 * 4
  list(aic_full = aic_full, aic_restricted = aic_restricted,
       delta_aic = aic_restricted - aic_full)
}
