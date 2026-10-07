# ============================================================
# fun/diagnostico_datos.R
#
# Diagnostico de observaciones atipicas e influyentes para el modelo
# mecanistico de Okuyama. Dos herramientas complementarias:
#
# (1) ANTES del ajuste (rapido, segundos): cribado_atipicos()
#     Ajusta una Beta-Binomial cuya media es la respuesta funcional de
#     Okuyama en forma cerrada (fun/betabinom_model.R; la misma que usa
#     Okuyama 2026 como modelo de comparacion) y, para cada ensayo,
#     calcula la probabilidad predictiva de su resultado con el modelo
#     ajustado SIN ese ensayo (leave-one-out). Un ensayo con
#     probabilidad muy baja (cola < alpha/n, Bonferroni) es atipico.
#     Ademas mide cuanto cambia la sobredispersion (rho) al quitarlo:
#     en el modelo mecanistico la sobredispersion entre ensayos la
#     absorbe k (forma de la Gamma del tiempo de busqueda), y un k muy
#     chico vuelve a z casi no identificable por arriba (Ensayo 1,
#     escenario C; Ensayo 2, D2 y R/04k).
#
# (2) DESPUES del ajuste (segundos con el motor C++): residuos_pit()
#     Residuos cuantilicos aleatorizados (Dunn & Smyth 1996) a partir de
#     la distribucion predictiva SIMULADA con los parametros ajustados
#     (motor_okuyama_cpp). Bajo el modelo correcto son ~N(0,1); |r| > 3
#     senala ensayos mal explicados.
#
# Ninguna de las dos quita datos: senalan ensayos para revisar
# (errores de carga, hembras defectuosas, etc.) y para un analisis de
# sensibilidad (re-ajustar sin ellos, como R/04k).
# ============================================================

# Necesita: fun/betabinom_model.R (fit_betabinom, mu_okuyama) y, para
# residuos_pit(), fun/motor_rcpp.R con el motor cargado.

#' Probabilidades Beta-Binomiales de 0..n
.pmf_bb <- function(n, pr, rho) {
  pr <- min(max(pr, 1e-9), 1 - 1e-9); rho <- max(rho, 1e-6)
  al <- pr * (1 - rho) / rho; be <- (1 - pr) * (1 - rho) / rho
  extraDistr::dbbinom(0:n, size = n, alpha = al, beta = be)
}

#' Cribado de ensayos atipicos/influyentes antes del ajuste mecanistico
#'
#' @param datos data.frame con columnas dens y par
#' @param T duracion del ensayo (misma unidad que en el ajuste)
#' @param alpha nivel global; cada ensayo se marca si su p (dos colas,
#'   leave-one-out) < alpha / n
#' @return data.frame ordenado por p, con p_loo, el cambio relativo en
#'   la sobredispersion al quitar el ensayo y las marcas
cribado_atipicos <- function(datos, T, alpha = 0.05) {
  n <- nrow(datos)
  ajuste <- fit_betabinom(datos, T)
  rho_todos <- ajuste$par[["rho"]]
  filas <- lapply(seq_len(n), function(i) {
    f <- fit_betabinom(datos[-i, ], T)
    p <- f$par
    N <- datos$dens[i]; y <- datos$par[i]
    pmf <- .pmf_bb(N, mu_okuyama(N, p[["a"]], p[["h"]], p[["z"]], T) / N, p[["rho"]])
    p_baja <- sum(pmf[seq_len(y + 1)])          # P(Y <= y)
    p_alta <- sum(pmf[(y + 1):(N + 1)])         # P(Y >= y)
    data.frame(fila = i, dens = N, par = y,
               esperado_loo = sum((0:N) * pmf),
               p_loo = min(1, 2 * min(p_baja, p_alta)),
               rho_sin = p[["rho"]],
               cambio_rho = p[["rho"]] / rho_todos - 1,
               nll_sin = f$value)
  })
  r <- do.call(rbind, filas)
  r$atipico <- r$p_loo < alpha / n
  r$influye_dispersion <- r$cambio_rho < -0.5     # sin el ensayo, rho cae mas de la mitad
  r <- r[order(r$p_loo), ]
  attr(r, "rho_todos") <- rho_todos
  attr(r, "umbral_p") <- alpha / n
  r
}

#' Residuos cuantilicos aleatorizados con la predictiva simulada
#'
#' @param par c(a, h, k, s) en la parametrizacion del motor (s = desvio natural)
#' @param z valor de z (p. ej. z_hat del perfil)
#' @param nsim ensayos simulados por densidad
residuos_pit <- function(datos, par, z, T, nsim = 1e5, semilla = NULL) {
  if (!is.null(semilla)) set.seed(semilla)
  r <- numeric(nrow(datos)); pb <- pa <- numeric(nrow(datos))
  for (N in unique(datos$dens)) {
    sim <- motor_okuyama_cpp(as.integer(N), par[[1]], par[[2]], z, par[[3]], par[[4]],
                             T, as.integer(nsim), 0L, TRUE)
    Fcum <- cumsum(tabulate(sim + 1L, nbins = N + 1L)) / nsim
    idx <- which(datos$dens == N)
    for (i in idx) {
      y <- datos$par[i]
      lo <- if (y == 0) 0 else Fcum[y]; hi <- Fcum[y + 1]
      u <- stats::runif(1, lo, hi)
      r[i] <- stats::qnorm(min(max(u, 1e-12), 1 - 1e-12))
      pb[i] <- hi; pa[i] <- 1 - lo
    }
  }
  data.frame(fila = seq_len(nrow(datos)), dens = datos$dens, par = datos$par,
             p_baja = pb, p_alta = pa, residuo = r, marca = abs(r) > 3)
}
