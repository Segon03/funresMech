# ============================================================
# fun/simulate_ensayo1.R
#
# Generador de datos sinteticos para el Ensayo 1 (recuperacion de
# parametros, Seccion 3 del plan de validacion). Usa la misma
# expectativa mecanistica (mu_okuyama, en betabinom_model.R) que
# el ajuste beta-binomial del Ensayo 2, y muestrea los conteos
# observados de una Beta-Binomial con dispersion rho (rho=0 ->
# Binomial), tal como describe la Seccion 2.4 original del
# manuscrito.
#
# NOTA: este archivo no estaba en la lista original de fun/ (que
# solo preveia betabinom_model.R); se agrego para no duplicar la
# logica de simulacion entre 00_benchmark.R y 02_ensayo1_simulacion.R.
# Ver README.md, seccion "Notas sobre la estructura".
# ============================================================

# mu_okuyama() vive en betabinom_model.R - quien use esta funcion
# debe haber corrido antes: source("fun/betabinom_model.R")

#' Simula un dataset del Ensayo 1 para un nivel de sobre-dispersion rho
#'
#' @param rho dispersion beta-binomial (0 = binomial, sin sobre-dispersion)
#' @param a,h,z parametros verdaderos (por defecto los del manuscrito y
#'   de Okuyama 2026, Fig. 1: a=0.3, h=0.02, z=1)
#' @param T duracion del ensayo (1 unidad de tiempo - ver nota de
#'   unidades, Seccion 3.4 del plan)
#' @param x_levels niveles de densidad de hospedadores (por defecto
#'   5,10,...,100 como en el manuscrito)
#' @param n_rep replicas por nivel de densidad (por defecto 10)
#'
#' @return data.frame con columnas dens, par (listo para
#'   funresMech:::fit_full())
simulate_dataset_ensayo1 <- function(rho, a = 0.3, h = 0.02, z = 1,
                                      T = 1, x_levels = seq(5, 100, by = 5),
                                      n_rep = 10) {
  do.call(rbind, lapply(x_levels, function(x) {
    mu <- mu_okuyama(x, a, h, z, T)
    pr <- min(max(mu / x, 1e-9), 1 - 1e-9)

    if (rho <= 0) {
      y <- rbinom(n_rep, size = x, prob = pr)
    } else {
      alpha <- pr * (1 - rho) / rho
      beta  <- (1 - pr) * (1 - rho) / rho
      y <- extraDistr::rbbinom(n_rep, size = x, alpha = alpha, beta = beta)
    }
    data.frame(dens = x, par = y)
  }))
}
