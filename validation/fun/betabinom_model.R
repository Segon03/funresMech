# ============================================================
# fun/betabinom_model.R
#
# Unica pieza de MODELO ESTADISTICO que no existe en el paquete
# funresMech: se verifico contra el codigo fuente real
# (github.com/Segon03/funresMech, v1.0.4 - DESCRIPTION e Imports,
# y toda la carpeta R/) que el paquete NO incluye VGAM, extraDistr,
# gamlss.dist ni ninguna implementacion de beta-binomial propia.
# Necesaria para el Ensayo 2 (comparacion mecanistico vs. beta-
# binomial, Seccion 4.3 del plan de validacion).
# ============================================================

#' Valor esperado de hosts parasitados bajo el modelo de Okuyama
#' (ecuaciones 1-3 del manuscrito / Okuyama 2026, Eq. 1-3)
#'
#' @param x densidad de hospedadores (dens)
#' @param a tasa de ataque
#' @param h tiempo de manipulacion medio
#' @param z exponente de escalado por densidad
#' @param T duracion del ensayo (misma unidad que a y h; ver nota de
#'   unidades en el plan de validacion, Seccion 3.4/4.4: usar T=1,
#'   no las 24 horas literales del protocolo)
mu_okuyama <- function(x, a, h, z, T = 1) {
  f <- a * x^z / (1 + a * h * x^z)
  x * (1 - exp(-f * T / x))
}

#' Log-verosimilitud negativa beta-binomial para el modelo de
#' Okuyama con dispersion rho (rho -> 0 equivale al binomial)
#'
#' @param par vector c(a, h, z, rho)
#' @param data_spp data.frame con columnas dens y par (mismos
#'   nombres que espera funresMech:::fit_full())
negloglik_betabinom <- function(par, data_spp, T = 1) {
  a <- par[1]; h <- par[2]; z <- par[3]; rho <- par[4]
  if (a <= 0 || h <= 0 || z <= 0) return(1e10)

  mu <- mu_okuyama(data_spp$dens, a, h, z, T)
  pr <- mu / data_spp$dens
  pr <- pmin(pmax(pr, 1e-9), 1 - 1e-9)   # evita 0/1 exactos (borde de la Beta)
  rho <- max(rho, 1e-6)

  alpha <- pr * (1 - rho) / rho
  beta  <- (1 - pr) * (1 - rho) / rho

  ll <- extraDistr::dbbinom(data_spp$par, size = data_spp$dens,
                             alpha = alpha, beta = beta, log = TRUE)
  if (any(!is.finite(ll))) return(1e10)
  -sum(ll)
}

#' Ajuste beta-binomial por maxima verosimilitud (L-BFGS-B)
#'
#' FIX (revision post-hoc, ver notas de Claude en el chat): el ajuste
#' de D3 con un unico arranque c(a=0.1,h=0.05,z=1,rho=0.01) terminaba
#' con convergence=52 ("ERROR: ABNORMAL_TERMINATION_IN_LNSRCH") - una
#' falla tipica de L-BFGS-B cerca del borde inferior de rho (rho->1e-6,
#' equivalente al binomial sin sobredispersion) donde la superficie de
#' verosimilitud queda mal condicionada (ver el Hessian de
#' ensayo2_bb_D3.rds: elementos de h en el orden de 1e6, sintoma de mal
#' escalado/casi-colinealidad ahi). El $value de esa corrida NO es un
#' minimo confiable, y de el depende aic_bb / delta_aic_bb_vs_mech /
#' favorece_mecanistico / reproduce_okuyama2026 para D3 - la conclusion
#' final del ensayo para ese dataset estaba construida sobre un ajuste
#' que el propio optim() reporto como fallido.
#'
#' Ahora se prueban varios arranques con L-BFGS-B y se toma el mejor
#' que SI converja (convergence==0); si ninguno converge, se cae a
#' Nelder-Mead (sin gradiente, mas robusto en superficies mal
#' condicionadas o cerca de bordes) y se pule el resultado con un
#' L-BFGS-B final partiendo de ahi.
#'
#' Devuelve un objeto tipo optim(): $par, $value (nll minimo),
#' $convergence, $hessian.
fit_betabinom <- function(data_spp, T = 1) {
  lower <- c(1e-4, 1e-4, 0.1, 1e-6)
  upper <- c(5, 5, 5, 0.999)

  starts <- list(
    c(a = 0.1, h = 0.05, z = 1,   rho = 0.01),
    c(a = 0.5, h = 0.02, z = 1.5, rho = 0.05),
    c(a = 1.0, h = 0.01, z = 2,   rho = 0.10),
    c(a = 0.3, h = 0.10, z = 1,   rho = 0.20),
    c(a = 0.05, h = 0.03, z = 0.8, rho = 0.001)
  )

  intentos <- lapply(starts, function(st) {
    tryCatch(
      optim(par = st, fn = negloglik_betabinom, data_spp = data_spp, T = T,
            method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE),
      error = function(e) NULL
    )
  })

  ok <- Filter(function(f) !is.null(f) && f$convergence == 0, intentos)

  if (length(ok) > 0) {
    return(ok[[which.min(vapply(ok, function(f) f$value, numeric(1)))]])
  }

  # Ningun arranque convergio con L-BFGS-B: respaldo con Nelder-Mead
  # (sin restricciones explicitas de caja, por eso se parte del mejor
  # intento aunque no haya convergido) y se pule con L-BFGS-B despues.
  mejor_intento <- intentos[[which.min(vapply(intentos, function(f) {
    if (is.null(f)) Inf else f$value
  }, numeric(1)))]]

  fit_nm <- optim(par = mejor_intento$par, fn = negloglik_betabinom,
                   data_spp = data_spp, T = T, method = "Nelder-Mead",
                   control = list(maxit = 5000))

  fit_final <- tryCatch(
    optim(par = pmin(pmax(fit_nm$par, lower), upper),
          fn = negloglik_betabinom, data_spp = data_spp, T = T,
          method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE),
    error = function(e) fit_nm
  )
  fit_final
}

#' AIC generico: AIC = 2*k + 2*NLL
aic_value <- function(negloglik, k_params) 2 * k_params + 2 * negloglik
