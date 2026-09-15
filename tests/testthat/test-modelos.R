# Tests de las funciones internas de calculo del modelo de Okuyama.
#
# CRAN pidio tests porque la interfaz 'shiny' no se puede verificar de
# forma automatica: estos tests cubren las funciones no exportadas que
# hacen el trabajo real (simulate_trial, negloglik_fixed_z, fit_full),
# usando funresMech::: para acceder a ellas.
#
# Reglas seguidas aca: la semilla se fija en el test (nunca dentro de la
# funcion testeada), itermax/NP/n_sim se mantienen chicos para no exceder
# el tiempo de check de CRAN, y no hay escritura a disco ni library().

test_that("simulate_trial() devuelve un entero valido dentro de rango", {
  # El resultado nunca puede ser negativo ni superar la densidad ofrecida (x)
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
  # Con x = 0 no hay hospedadores para atacar: el resultado debe ser 0
  set.seed(42)
  resultado <- funresMech:::simulate_trial(
    x = 0, T = 1, a = 0.5, h = 0.1, z = 1, k = 1, s = 0
  )
  expect_equal(resultado, 0)
})

test_that("simulate_trial() es estocastico: semillas distintas dan resultados distintos", {
  # Al ser un proceso estocastico (Gamma + Lognormal), semillas distintas
  # deben producir al menos un resultado diferente entre si
  args <- list(x = 20, T = 5, a = 0.3, h = 0.2, z = 1, k = 1, s = 0.3)
  set.seed(1)
  r1 <- do.call(funresMech:::simulate_trial, args)
  set.seed(2)
  r2 <- do.call(funresMech:::simulate_trial, args)
  set.seed(3)
  r3 <- do.call(funresMech:::simulate_trial, args)
  expect_true(length(unique(c(r1, r2, r3))) > 1)
})

test_that("negloglik_fixed_z() devuelve un escalar numerico finito", {
  # Con parametros y datos razonables, la log-verosimilitud negativa debe
  # ser un unico numero finito
  set.seed(42)
  data_spp <- data.frame(
    dens = rep(c(5, 10), each = 5),
    par  = c(1, 2, 1, 3, 2, 4, 5, 3, 6, 4)
  )
  nll <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 0.3, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = data_spp, T = 1, n_sim = 50
  )
  expect_type(nll, "double")
  expect_length(nll, 1)
  expect_true(is.finite(nll))
})

test_that("negloglik_fixed_z() devuelve un valor grande fuera del dominio valido", {
  # Con a <= 0 (fuera del dominio del modelo) la funcion no debe fallar,
  # sino devolver un valor centinela grande (o Inf)
  data_spp <- data.frame(dens = c(5, 10), par = c(1, 2))
  nll <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = -1, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = data_spp, T = 1, n_sim = 50
  )
  expect_true(is.infinite(nll) || nll >= 1e9)
})

test_that("negloglik_fixed_z() es sensible a los parametros", {
  # Dos juegos de parametros bien distintos deben dar NLL distintos
  set.seed(42)
  data_spp <- data.frame(
    dens = rep(c(5, 10), each = 5),
    par  = c(1, 2, 1, 3, 2, 4, 5, 3, 6, 4)
  )
  nll_1 <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 0.2, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = data_spp, T = 1, n_sim = 50
  )
  nll_2 <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 1.5, h = 0.4, k = 3, s = 0.4),
    z_fixed = 2, data_spp = data_spp, T = 1, n_sim = 50
  )
  expect_false(isTRUE(all.equal(nll_1, nll_2)))
})

test_that("fit_full() devuelve la estructura esperada", {
  # itermax y NP chicos: CRAN limita el tiempo de check. Se corre igual
  # (no se saltea) porque con estos valores el ajuste es rapido.
  # DEoptim avisa (via warning()) que NP deberia ser >= 10x la cantidad de
  # parametros; es solo una recomendacion de convergencia, no un error, y
  # a proposito la ignoramos aca para mantener el test rapido.
  set.seed(42)
  data_spp <- data.frame(
    dens = rep(c(5, 10), each = 5),
    par  = c(1, 2, 1, 3, 2, 4, 5, 3, 6, 4)
  )
  ajuste <- suppressWarnings(funresMech:::fit_full(
    data_spp = data_spp, T_exp = 1,
    itermax = 20, NP = 20, reltol = 1e-2, n_sim_profile = 50
  ))

  expect_type(ajuste, "list")
  expect_named(ajuste, c("par", "nll"))
  expect_named(ajuste$par, c("a", "h", "z", "k", "s"))
  expect_true(is.finite(ajuste$nll))

  # Los estimados deben caer dentro de los limites del optimizador
  lower <- c(a = 0.001, h = 0.001, z = 0.5, k = 0.5, s = 0.001)
  upper <- c(a = 2.0,   h = 0.5,   z = 3.0, k = 5.0, s = 0.5)
  expect_true(all(ajuste$par >= lower[names(ajuste$par)]))
  expect_true(all(ajuste$par <= upper[names(ajuste$par)]))
})

test_that("la plantilla del reporte esta instalada y accesible", {
  # system.file() debe encontrar report_template.Rmd dentro del paquete
  # instalado, sin depender de getwd() ni de rutas relativas
  ruta <- system.file("app", "report_template.Rmd", package = "funresMech")
  expect_true(nzchar(ruta))
  expect_true(file.exists(ruta))
})
