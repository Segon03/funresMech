# Tests de las funciones internas de calculo del modelo de Okuyama.
#
# CRAN pidio tests porque la interfaz 'shiny' no se puede verificar de
# forma automatica: estos tests cubren las funciones no exportadas que
# hacen el trabajo real (motor C++, verosimilitud simulada, perfil), usando
# funresMech::: para acceder a ellas.
#
# Reglas seguidas aca: la semilla se fija en el test (nunca dentro de la
# funcion testeada), itermax/NP/n_sim se mantienen chicos para no exceder
# el tiempo de check de CRAN, y no hay escritura a disco ni library().
# Todos los tests se pensaron para correr en pocos segundos; los dos mas
# costosos (screening, perfil completo) se saltean en CRAN.

datos_chicos <- data.frame(
  dens = rep(c(5, 10), each = 5),
  par  = c(1, 2, 1, 3, 2, 4, 5, 3, 6, 4)
)

# ---------------------------------------------------------------------------
# Motor C++ vs. implementacion de referencia en R
# ---------------------------------------------------------------------------

test_that("simulate_trial() devuelve un entero valido dentro de rango", {
  set.seed(42)
  x <- 10
  resultado <- funresMech:::simulate_trial(
    x = x, T = 1, a = 0.5, h = 0.1, z = 1, k = 1, s = 0
  )
  expect_type(resultado, "integer")
  expect_length(resultado, 1)
  expect_gte(resultado, 0)
  expect_lte(resultado, x)
})

test_that("simulate_trial() con densidad 0 devuelve 0", {
  set.seed(42)
  resultado <- funresMech:::simulate_trial(
    x = 0, T = 1, a = 0.5, h = 0.1, z = 1, k = 1, s = 0
  )
  expect_equal(resultado, 0)
})

test_that("simulate_trial() es estocastico: semillas distintas dan resultados distintos", {
  args <- list(x = 20, T = 5, a = 0.3, h = 0.2, z = 1, k = 1, s = 0.3)
  set.seed(1)
  r1 <- do.call(funresMech:::simulate_trial, args)
  set.seed(2)
  r2 <- do.call(funresMech:::simulate_trial, args)
  set.seed(3)
  r3 <- do.call(funresMech:::simulate_trial, args)
  expect_true(length(unique(c(r1, r2, r3))) > 1)
})

test_that("el motor C++ reproduce exactamente a simulate_trial() con la misma semilla", {
  # Sin corte por saturacion consumen los mismos numeros aleatorios, en el
  # mismo orden: los resultados deben ser identicos (con s = 0 y con s > 0)
  casos <- list(
    list(x = 30, T = 2,  a = 0.5,   h = 0.2, z = 1.5, k = 0.5,  s = 0),
    list(x = 8,  T = 24, a = 0.01,  h = 0.6, z = 3,   k = 0.05, s = 0.8),
    list(x = 50, T = 24, a = 0.001, h = 0.5, z = 2.5, k = 0.04, s = 0.6)
  )
  for (cs in casos) {
    n <- 100L
    set.seed(7)
    motor <- funresMech:::motor_okuyama_cpp(
      as.integer(cs$x), cs$a, cs$h, cs$z, cs$k, cs$s, cs$T, n,
      corte_saturacion = FALSE
    )
    set.seed(7)
    ref <- vapply(seq_len(n), function(i) {
      funresMech:::simulate_trial(cs$x, cs$T, cs$a, cs$h, cs$z, cs$k, cs$s)
    }, integer(1))
    expect_identical(as.integer(motor), as.integer(ref))
  }
})

test_that("el corte por saturacion no altera la distribucion simulada", {
  # Con parametros que saturan a menudo (todos los huespedes parasitados)
  n <- 20000L
  set.seed(11)
  con <- funresMech:::motor_okuyama_cpp(6L, 2, 0.05, 1, 1, 0.05, 5, n,
                                        corte_saturacion = TRUE)
  set.seed(12)
  sin <- funresMech:::motor_okuyama_cpp(6L, 2, 0.05, 1, 1, 0.05, 5, n,
                                        corte_saturacion = FALSE)
  tab <- rbind(tabulate(con + 1L, nbins = 7), tabulate(sin + 1L, nbins = 7))
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  expect_gt(suppressWarnings(chisq.test(tab)$p.value), 0.001)
  expect_gt(mean(con == 6L), 0.1)   # el escenario efectivamente satura
})

test_that("el corte por saturacion evita la trampa computacional (h ~ 0, k ~ 0)", {
  # Sin el corte, 10 ensayos a H = 50 con h = k = 1e-6 requieren decenas de
  # millones de vueltas del bucle; con el corte, unas pocas miles
  vueltas <- funresMech:::contar_vueltas_cpp(50L, 1, 1e-6, 1, 1e-6, 0, 24, 10L,
                                             corte_saturacion = TRUE)
  expect_lt(vueltas, 1e5)
  set.seed(1)
  r <- funresMech:::motor_okuyama_cpp(50L, 1, 1e-6, 1, 1e-6, 0, 24, 10L, TRUE)
  expect_true(all(r == 50L))
})

test_that("el motor rechaza parametros invalidos con un error claro", {
  expect_error(funresMech:::motor_okuyama_cpp(0L, 1, 1, 1, 1, 0, 1, 10L), "H must be")
  expect_error(funresMech:::motor_okuyama_cpp(5L, -1, 1, 1, 1, 0, 1, 10L), "invalid parameters")
  expect_error(funresMech:::motor_okuyama_cpp(5L, 1, 0, 1, 1, 0.5, 1, 10L), "h must be")
})

# ---------------------------------------------------------------------------
# Verosimilitud simulada
# ---------------------------------------------------------------------------

test_that("simulate_distribution() devuelve probabilidades validas", {
  set.seed(3)
  p <- funresMech:::simulate_distribution(x = 10, T = 1, a = 0.5, h = 0.1,
                                          z = 1, k = 1, s = 0.05, n_sim = 500)
  expect_length(p, 11)
  expect_named(p, as.character(0:10))
  expect_true(all(p > 0))
  expect_equal(sum(p), 1, tolerance = 1e-6)
})

test_that("negloglik_fixed_z() devuelve un escalar numerico finito", {
  set.seed(42)
  nll <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 0.3, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = datos_chicos, T = 1, n_sim = 50
  )
  expect_type(nll, "double")
  expect_length(nll, 1)
  expect_true(is.finite(nll))
})

test_that("negloglik_fixed_z() devuelve un valor grande fuera del dominio valido", {
  data_spp <- data.frame(dens = c(5, 10), par = c(1, 2))
  nll <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = -1, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = data_spp, T = 1, n_sim = 50
  )
  expect_true(is.infinite(nll) || nll >= 1e9)
  # parasitismo mayor que la densidad: dato imposible
  nll2 <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 0.3, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = data.frame(dens = 5, par = 7), T = 1, n_sim = 50
  )
  expect_gte(nll2, 1e9)
})

test_that("negloglik_fixed_z() es sensible a los parametros", {
  set.seed(42)
  nll_1 <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 0.2, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = datos_chicos, T = 1, n_sim = 50
  )
  nll_2 <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 1.5, h = 0.4, k = 3, s = 0.4),
    z_fixed = 2, data_spp = datos_chicos, T = 1, n_sim = 50
  )
  expect_false(isTRUE(all.equal(nll_1, nll_2)))
})

test_that("negloglik_fixed_z() es reproducible con la semilla del usuario", {
  f <- function() funresMech:::negloglik_fixed_z(
    c(a = 0.3, h = 0.1, k = 1, s = 0.1), z_fixed = 1,
    data_spp = datos_chicos, T = 1, n_sim = 200)
  set.seed(99); n1 <- f()
  set.seed(99); n2 <- f()
  expect_identical(n1, n2)
})

# ---------------------------------------------------------------------------
# Reparametrizacion log
# ---------------------------------------------------------------------------

test_that("la reparametrizacion log es de ida y vuelta exacta", {
  H_ref <- funresMech:::.h_ref(data.frame(dens = c(2, 10, 50)))
  expect_equal(H_ref, exp(mean(log(c(2, 10, 50)))))
  par_nat <- c(a = 3.9e-5, h = 0.626, k = 0.026, s = 0.669)
  for (z in c(0.5, 2.5, 4.5, 20)) {
    theta <- funresMech:::.nat_to_theta(par_nat, z, H_ref)
    back  <- funresMech:::.theta_to_nat(theta, z, H_ref)
    expect_equal(unname(back), unname(par_nat), tolerance = 1e-12)
  }
})

test_that("las cotas del optimizador cubren los valores de D2 y D3 de Okuyama", {
  b <- funresMech:::.log_bounds()
  H_ref <- 15.4
  # a ~ 1e-5 a z = 6, h ~ 0.6, k ~ 0.05, s ~ 0.9 deben quedar dentro
  for (z in c(1, 4.5, 6, 20)) {
    th <- funresMech:::.nat_to_theta(c(a = 1e-3 / H_ref^(z - 1), h = 0.6, k = 0.05, s = 0.9),
                                     z, H_ref)
    expect_true(all(th >= b$lower & th <= b$upper))
  }
})

# ---------------------------------------------------------------------------
# Intervalo de verosimilitud
# ---------------------------------------------------------------------------

test_that("el umbral del IC95 es qchisq(0.95, 1) / 2 y no qchisq(0.95, 1)", {
  expect_equal(funresMech:::.umbral_ic95(), 1.920729, tolerance = 1e-6)
})

test_that(".profile_ci() recupera el intervalo exacto de un perfil cuadratico", {
  # NLL(z) = 0.5 * ((z - 3) / 0.5)^2: cruza el umbral en 3 +- 0.5 * 1.95996
  z <- seq(0.5, 6, by = 0.01)
  nll <- 100 + 0.5 * ((z - 3) / 0.5)^2
  ci <- funresMech:::.profile_ci(z, nll)
  expect_equal(ci$z_hat, 3)
  expect_equal(ci$z_low,  3 - 0.5 * qnorm(0.975), tolerance = 1e-3)
  expect_equal(ci$z_high, 3 + 0.5 * qnorm(0.975), tolerance = 1e-3)
  expect_false(ci$low_censored)
  expect_false(ci$high_censored)
})

test_that(".profile_ci() informa un limite abierto cuando el perfil no cruza el umbral", {
  z <- seq(0.5, 6, by = 0.5)
  nll <- 100 + 0.5 * pmin(z - 2, 0)^2 * 4 + 0.02 * pmax(z - 2, 0)   # plano a la derecha
  ci <- funresMech:::.profile_ci(z, nll)
  expect_true(ci$high_censored)
  expect_true(is.na(ci$z_high))
  expect_false(ci$low_censored)
  expect_true(is.finite(ci$z_low))
})

test_that("un limite superior abierto solo se llama 'no identificado' si la grilla llego al tope", {
  expect_match(funresMech:::.high_open_text(5), "not reached within the grid")
  expect_match(funresMech:::.high_open_text(20), "not identified")
})

test_that(".refine_points() agrega puntos solo dentro de un intervalo ancho", {
  # grilla 1..5 en pasos de 0.5 + salto a 8: el limite superior (5.12) cae en (5, 8)
  z <- c(1, 1.5, 2, 2.5, 3, 4, 5, 8)
  nuevos <- funresMech:::.refine_points(z, z_low = 1.53, z_high = 5.12, step0 = 0.5)
  expect_equal(nuevos, c(5.5, 6, 6.5, 7, 7.5))
  # un limite en un intervalo de ancho normal (2, 2.5) no agrega nada
  expect_length(funresMech:::.refine_points(z, z_low = 1.53, z_high = 2.2, step0 = 0.5), 0)
  # limites abiertos (NA): nada que refinar
  expect_length(funresMech:::.refine_points(z, NA_real_, NA_real_, 0.5), 0)
})

test_that(".nll_at_z() interpola linealmente y devuelve NA fuera de la grilla", {
  z <- c(0.5, 1.5, 2.5)
  nll <- c(10, 6, 8)
  expect_equal(funresMech:::.nll_at_z(z, nll, 1.5), 6)
  expect_equal(funresMech:::.nll_at_z(z, nll, 1), 8)
  expect_true(is.na(funresMech:::.nll_at_z(z, nll, 3)))
})

# ---------------------------------------------------------------------------
# Ajuste
# ---------------------------------------------------------------------------

test_that("fit_z_fixed() devuelve la estructura esperada", {
  # itermax y NP chicos: CRAN limita el tiempo de check. DEoptim avisa que
  # NP deberia ser >= 10x la cantidad de parametros; es solo una
  # recomendacion de convergencia y se ignora aca para mantener el test rapido.
  set.seed(42)
  ajuste <- suppressWarnings(funresMech:::fit_z_fixed(
    z = 1, data_spp = datos_chicos, T_exp = 1,
    n_sim = 50, itermax = 10, NP = 20, reltol = 1e-2, n_reeval = 2
  ))
  expect_type(ajuste, "list")
  expect_named(ajuste$par, c("a", "h", "k", "s"))
  expect_true(all(ajuste$par[c("a", "h", "k")] > 0))
  expect_gte(ajuste$par[["s"]], 0)
  expect_true(is.finite(ajuste$nll))
  expect_true(is.finite(ajuste$nll_reeval))
  expect_length(ajuste$at_bound, 4)
})

test_that("fit_profile() devuelve el perfil, z_hat, IC y AIC", {
  skip_on_cran()
  set.seed(42)
  ajuste <- suppressWarnings(funresMech:::fit_profile(
    data_spp = datos_chicos, T_exp = 1, z_grid = c(0.5, 1, 2),
    n_sim = 50, itermax = 10, NP = 20, reltol = 1e-2, n_reeval = 2,
    extend_z = FALSE
  ))
  expect_named(ajuste, c("profile", "par", "nll", "aic", "ci", "extended", "refined", "notes"))
  expect_named(ajuste$par, c("a", "h", "z", "k", "s"))
  expect_equal(nrow(ajuste$profile), 3)
  expect_true(all(is.finite(ajuste$profile$nll)))
  # z_hat es el minimo del perfil y la NLL informada es la del minimo
  expect_equal(ajuste$par[["z"]], ajuste$profile$z[which.min(ajuste$profile$nll)])
  expect_equal(ajuste$nll, min(ajuste$profile$nll))
  expect_equal(ajuste$aic$AIC_full, 2 * ajuste$nll + 10)
  expect_false(is.na(ajuste$aic$AIC_restricted))   # z = 1 esta en la grilla
})

test_that("fit_profile() agrega z = 1 a la grilla y valida los datos", {
  skip_on_cran()
  set.seed(1)
  ajuste <- suppressWarnings(funresMech:::fit_profile(
    data_spp = datos_chicos, T_exp = 1, z_grid = c(0.5, 1.5),
    n_sim = 30, itermax = 5, NP = 12, n_reeval = 1, extend_z = FALSE))
  expect_true(1 %in% ajuste$profile$z)
  expect_error(funresMech:::fit_profile(
    data_spp = data.frame(dens = 5, par = 9), T_exp = 1, z_grid = 1),
    "between 0 and the host density")
})

# ---------------------------------------------------------------------------
# Screening de ensayos atipicos
# ---------------------------------------------------------------------------

test_that("screen_outliers() marca un ensayo atipico y no marca datos coherentes", {
  skip_on_cran()
  # D2 (Kishani Farahani & Goldansaz 2013, FoRAGE): el unico ensayo atipico
  # es 2 parasitados de 10 hospedadores (vecinos: 85-95 %)
  d2 <- data.frame(
    dens = c(50, 50, 2, 2, 2, 4, 4, 6, 6, 6, 8, 8, 10, 10, 10, 15, 15, 15, 15,
             20, 20, 20, 20, 20, 25, 25, 25, 25, 25, 25, 30, 30, 30, 30, 30, 30,
             30, 30, 30, 30, 35, 35, 35, 35, 35, 35, 35, 35, 40, 40, 40, 40, 40,
             40, 40, 40, 40, 45, 45, 45, 45, 45, 45, 45, 45, 45, 50, 50, 50, 50,
             50, 50, 50),
    par  = c(31, 32, 0, 1, 2, 3, 4, 4, 5, 6, 7, 8, 2, 8, 10, 12, 13, 14, 15,
             15, 16, 17, 18, 19, 15, 18, 19, 20, 21, 22, 16, 17, 18, 19, 20, 21,
             22, 23, 24, 20, 21, 22, 23, 24, 25, 26, 27, 29, 21, 23, 24, 25, 26,
             27, 28, 30, 31, 20, 22, 23, 24, 25, 27, 28, 29, 30, 25, 26, 27, 28,
             29, 30, 33)
  )
  expect_equal(nrow(d2), 73)
  scr <- funresMech:::screen_outliers(d2, T_exp = 24)
  expect_equal(sum(scr$flagged), 1)
  top <- scr[1, ]
  expect_equal(c(top$dens, top$par), c(10, 2))
  expect_true(top$influences_dispersion)

  # datos generados sin atipicos
  set.seed(5)
  N <- rep(c(5, 10, 20, 40), each = 6)
  y <- rbinom(length(N), N, 0.5)
  scr0 <- funresMech:::screen_outliers(data.frame(dens = N, par = y), T_exp = 24)
  expect_equal(sum(scr0$flagged), 0)
})

# ---------------------------------------------------------------------------
# Plantilla del reporte
# ---------------------------------------------------------------------------

test_that("la plantilla del reporte esta instalada y accesible", {
  ruta <- system.file("app", "report_template.Rmd", package = "funresMech")
  expect_true(nzchar(ruta))
  expect_true(file.exists(ruta))
})

test_that(".workers_ok() is TRUE under a sequential plan", {
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  expect_true(funresMech:::.workers_ok())
})
