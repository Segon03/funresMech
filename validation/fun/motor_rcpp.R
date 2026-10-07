# ============================================================
# fun/motor_rcpp.R
#
# Capa R del motor Rcpp (fun/motor_okuyama.cpp). Reemplaza, SOLO
# dentro de este proyecto de validacion, al motor R puro de
# funresMech v1.0.4 (simulate_trial / simulate_distribution /
# negloglik_fixed_z). El paquete funresMech NO se modifica hasta que
# la validacion este cerrada (ver docs/MOTOR_RCPP_cambios_y_validacion.md).
#
# Contenido:
#   cargar_motor_rcpp()        compila (una vez, con cache) y carga el C++
#   nll_gamma_lnorm_S4()       COPIA LITERAL de Okuyama (2026) Script S4
#                              (referencia; usa el simulador literal S5)
#   nll_motor()                NLL rapida con el motor (B), parametros
#                              naturales (a, h, k, s)
#   negloglik_fixed_z_rcpp()   misma firma que funresMech:::negloglik_fixed_z
#                              (reemplazo directo para el paquete)
#   theta_a_nat(), nat_a_theta(), COTAS_LOG
#                              reparametrizacion en escala log
#   ajustar_z_fijo()           un punto del perfil de z: DEoptim +
#                              re-evaluacion de la NLL en el optimo
#
# Convenciones de parametros (IMPORTANTE para comparar con Okuyama):
#   a, h, z, k  : igual que Okuyama y funresMech.
#   s           : en modo "okuyama" (el que se usa en toda la
#                 validacion) es el DESVIO ESTANDAR del tiempo de
#                 manipulacion en escala NATURAL (horas), igual que la
#                 columna s de la Tabla S1 de Okuyama. En funresMech
#                 v1.0.4 s es el desvio en escala LOG (sdlog). Son la
#                 misma familia de distribuciones con otra
#                 parametrizacion: sdlog = sqrt(log(1 + (s/h)^2)).
# ============================================================

MOTOR_RCPP_VERSION <- "1.0 (2026-09-28)"

# ---- Compilacion / carga --------------------------------------------------

#' Compila fun/motor_okuyama.cpp (tarda ~20-60 s) y deja las funciones
#' C++ disponibles en el entorno global. Llamarla tambien DENTRO de cada
#' worker.
#'
#' CORRECCION (04-oct-2026): la cache de compilacion ahora es PROPIA DE
#' CADA PROCESO de R (dentro de tempdir()), no compartida en fun/.cache_rcpp.
#' Rcpp::sourceCpp() lee y reescribe un indice (cache.rds) en la carpeta de
#' cache en CADA llamada; con varios workers usando la misma carpeta, uno
#' podia leer el indice justo mientras otro lo reescribia y encontrarlo
#' vacio -> "Error en load(file = index_file): archivo de entrada vacio
#' (cero bytes)", que cortaba toda la corrida (30-sep con 04h y 03-oct con
#' 02c). Con una cache por proceso esa carrera no puede ocurrir. Costo:
#' cada worker compila una vez al arrancar (~20-60 s). El codigo C++ y los
#' resultados no cambian.
cargar_motor_rcpp <- function(archivo = file.path("fun", "motor_okuyama.cpp"),
                              cache   = file.path(tempdir(), "cache_rcpp_motor"),
                              verbose = FALSE) {
  # Si ya esta cargado Y funciona (un puntero de otra sesion, p.ej.
  # copiado a un worker de future, existe pero no funciona), no hacer nada.
  ya_ok <- tryCatch({
    f <- get("motor_okuyama_cpp", envir = globalenv(), inherits = FALSE)
    is.integer(f(2L, 1, 1, 1, 1, 0, 1, 1L, 0L, TRUE))
  }, error = function(e) FALSE)
  if (isTRUE(ya_ok)) return(invisible(TRUE))
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Falta el paquete Rcpp: install.packages('Rcpp')")
  }
  if (.Platform$OS.type == "windows" && !nzchar(Sys.which("make"))) {
    stop("No se encontro Rtools (no hay 'make' en el PATH). En Windows hace falta ",
         "instalar Rtools de la MISMA version que R (https://cran.r-project.org/bin/windows/Rtools/) ",
         "y reiniciar RStudio. Verificar con: Sys.which('make')")
  }
  if (!file.exists(archivo)) stop("No se encontro ", archivo,
                                  ". Correr desde la carpeta del proyecto (.Rproj).")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(archivo, cacheDir = cache, env = globalenv(),
                  verbose = verbose, rebuild = FALSE)
  invisible(TRUE)
}

# ---- (A) Referencia: COPIA LITERAL de Okuyama (2026) Script S4 ----------
# Cambios: el nombre (para no chocar con otros objetos) y el argumento
# 'simulador', que permite elegir el simulador LITERAL de Okuyama
# (parasitism_gamma_lnorm_cpp, Script S5; default) o su version con la
# correccion del RNGScope anidado (parasitism_gamma_lnorm_cpp_corr).
nll_gamma_lnorm_S4 <- function(p, H0, H1, T, z = NA, nsim = 10000, minp = NA, s.sqrt = FALSE,
                               simulador = NULL) {
  if (!is.null(simulador)) parasitism_gamma_lnorm_cpp <- simulador
  a <- p[1]
  h <- p[2]
  k <- p[3]
  if(s.sqrt==FALSE) s <- p[4]
  if(s.sqrt==TRUE)  s <- p[4]*p[4]
  if (is.na(z)) z <- p[5]
  if (is.na(minp)) minp <- .Machine$double.xmin

  if (any(p < 0)) return(NA)
  if(k==0) return(NA)

  unique.H0 <- unique(H0)
  log_likelihood <- 0

  for (h0_value in unique.H0) {
     observed_H1 <- H1[H0==h0_value]
    simulated_H1 <- parasitism_gamma_lnorm_cpp(h0_value, a, h, z, k, s, T, nsim)

    prob <- pmax(rowMeans(outer(observed_H1, simulated_H1, "==")), minp)

    log_likelihood <- log_likelihood + sum(log(prob))
  }
  return(-log_likelihood)
}

# ---- (B) NLL rapida con el motor del proyecto ----------------------------

.MODOS <- c(okuyama = 0L, funresMech = 1L)

#' NLL simulada (Okuyama 2026, ec. de la verosimilitud por simulacion).
#' @param par_nat c(a, h, k, s) en escala natural (s segun 'modo').
#' @param z valor fijo de z.
#' @param datos data.frame con columnas dens (H0) y par (H1).
#' @param modo "okuyama" (Script S5) o "funresMech" (motor v1.0.4).
#' @param corte_saturacion TRUE: corte exacto al saturar (ver .cpp).
#' Devuelve 1e10 si los parametros no son validos (convencion de
#' funresMech, para que DEoptim los descarte sin romperse).
nll_motor <- function(par_nat, z, datos, T, nsim = 10000,
                      modo = "okuyama", corte_saturacion = TRUE,
                      minp = .Machine$double.xmin) {
  a <- par_nat[[1]]; h <- par_nat[[2]]; k <- par_nat[[3]]; s <- par_nat[[4]]
  if (!is.finite(a) || !is.finite(h) || !is.finite(k) || !is.finite(s) ||
      a <= 0 || k <= 0 || h < 0 || s < 0 || (s > 0 && h <= 0)) return(1e10)
  m <- .MODOS[[modo]]
  # Okuyama recorre unique(H0) en el orden de los datos; funresMech en
  # orden creciente. Solo cambia el orden de consumo de la semilla.
  niveles <- if (m == 0L) unique(datos$dens) else sort(unique(datos$dens))
  nll <- 0
  for (x in niveles) {
    y   <- datos$par[datos$dens == x]
    sim <- motor_okuyama_cpp(as.integer(x), a, h, z, k, s, T, as.integer(nsim),
                             m, corte_saturacion)
    frec <- tabulate(sim + 1L, nbins = x + 1L) / nsim
    nll  <- nll - sum(log(pmax(frec[y + 1L], minp)))
  }
  nll
}

#' Reemplazo directo de funresMech:::negloglik_fixed_z() (misma firma,
#' mismos parametros, mismo resultado con la misma semilla - prueba 2
#' de R/00b), pero en C++. Pensado para la proxima version del paquete.
negloglik_fixed_z_rcpp <- function(par_vec, z_fixed, data_spp, T, n_sim = 1500,
                                   corte_saturacion = TRUE) {
  a <- par_vec[1]; h <- par_vec[2]; k <- par_vec[3]; s <- par_vec[4]
  if (a <= 0 || h <= 0 || k <= 0 || s < 0) return(1e10)
  nll_motor(c(a, h, k, s), z_fixed, data_spp, T, n_sim,
            modo = "funresMech", corte_saturacion = corte_saturacion)
}

# ---- Reparametrizacion en escala log --------------------------------------
#
# DEoptim muta y cruza SUMANDO diferencias entre vectores de la
# poblacion en la escala en que se le dan las cotas. En escala lineal
# [0, 5] un valor como a = 1e-5 (necesario a z alto) o k = 0.01 ocupa
# una fraccion infima del rango y practicamente nunca se explora bien.
# Por eso se busca en:
#
#   theta1 = log(lambda_ref),  lambda_ref = a * H_ref^z
#            (tasa de encuentro a la densidad de referencia H_ref:
#             hace que el rango de busqueda NO dependa de z; a z alto
#             'a' tiene que achicarse ~H_ref^z veces y esta
#             transformacion lo absorbe)
#   theta2 = log(h)
#   theta3 = log(k)
#   theta4 = u,  s = u^2   (igual que Okuyama, s.sqrt = TRUE; se deja
#            lineal porque s = 0 - manipulacion fija - tiene que ser
#            alcanzable)
#
# H_ref = media geometrica de las densidades distintas del dataset.
# Las cotas cubren y exceden las de Okuyama (0-5) en todo el rango
# plausible; si un ajuste termina pegado a una cota, ajustar_z_fijo()
# lo marca.

COTAS_LOG <- list(
  lower = c(log_lambda_ref = log(1e-3), log_h = log(1e-3), log_k = log(1e-5), u_s = 0),
  upper = c(log_lambda_ref = log(1e3),  log_h = log(5),    log_k = log(5),    u_s = 5)
)
# Cotas de Okuyama (Script S3): p en [0,5]^4 con s = p4^2
COTAS_OKUYAMA <- list(lower = c(a = 0, h = 0, k = 0, u_s = 0),
                      upper = c(a = 5, h = 5, k = 5, u_s = 5))

h_ref_de <- function(datos) exp(mean(log(unique(datos$dens))))

theta_a_nat <- function(theta, z, H_ref) {
  c(a = exp(theta[[1]]) / H_ref^z, h = exp(theta[[2]]),
    k = exp(theta[[3]]), s = theta[[4]]^2)
}
nat_a_theta <- function(par_nat, z, H_ref) {
  c(log_lambda_ref = log(par_nat[["a"]]) + z * log(H_ref),
    log_h = log(par_nat[["h"]]), log_k = log(par_nat[["k"]]),
    u_s = sqrt(par_nat[["s"]]))
}

# ---- Un punto del perfil de z ---------------------------------------------

#' Ajusta (a, h, k, s) con z fijo usando DEoptim, como Okuyama (Script
#' S3), y re-evalua la NLL en el optimo con semillas nuevas.
#'
#' @param escala "log" (reparametrizacion de arriba, recomendada) u
#'   "okuyama" (exactamente la busqueda de Okuyama: lineal, [0,5]^4,
#'   s = p4^2).
#' @param n_reeval, nsim_reeval la NLL que devuelve DEoptim (bestval) es
#'   el MINIMO de muchas evaluaciones ruidosas y por eso esta sesgada
#'   hacia abajo; se re-evalua n_reeval veces en el optimo con
#'   nsim_reeval simulaciones y se informa la media (nll_reeval) y su
#'   error estandar. Ambas se guardan; ver 04i.
ajustar_z_fijo <- function(z, datos, T, nsim = 10000,
                           escala = c("log", "okuyama"),
                           itermax = 200, NP = 40,
                           n_reeval = 5, nsim_reeval = nsim,
                           modo = "okuyama", corte_saturacion = TRUE,
                           cotas = NULL, trace = FALSE) {
  escala <- match.arg(escala)
  H_ref  <- h_ref_de(datos)
  t0 <- proc.time()[["elapsed"]]

  if (escala == "log") {
    if (is.null(cotas)) cotas <- COTAS_LOG
    fn <- function(th) nll_motor(theta_a_nat(th, z, H_ref), z, datos, T, nsim,
                                 modo, corte_saturacion)
  } else {
    if (is.null(cotas)) cotas <- COTAS_OKUYAMA
    fn <- function(p) {
      if (any(p < 0) || p[3] == 0 || p[1] == 0) return(1e10)
      nll_motor(c(p[1], p[2], p[3], p[4]^2), z, datos, T, nsim, modo, corte_saturacion)
    }
  }

  res <- DEoptim::DEoptim(fn, lower = cotas$lower, upper = cotas$upper,
                          control = DEoptim::DEoptim.control(itermax = itermax, NP = NP,
                                                             trace = trace))
  mejor <- res$optim$bestmem
  par_nat <- if (escala == "log") theta_a_nat(mejor, z, H_ref) else
    c(a = mejor[[1]], h = mejor[[2]], k = mejor[[3]], s = mejor[[4]]^2)

  reevals <- vapply(seq_len(n_reeval), function(i)
    nll_motor(par_nat, z, datos, T, nsim_reeval, modo, corte_saturacion), numeric(1))

  rango <- cotas$upper - cotas$lower
  en_cota <- (mejor <= cotas$lower + 1e-3 * rango) | (mejor >= cotas$upper - 1e-3 * rango)
  names(en_cota) <- names(cotas$lower)

  list(z = z, par = par_nat, theta = setNames(as.numeric(mejor), names(cotas$lower)),
       nll = res$optim$bestval,              # comparable con Okuyama (bestval de DEoptim)
       nll_reeval = mean(reevals),
       nll_reeval_se = if (n_reeval > 1) sd(reevals) / sqrt(n_reeval) else NA_real_,
       en_cota = en_cota, H_ref = H_ref,
       iter = res$optim$iter, nfeval = res$optim$nfeval,
       segundos = proc.time()[["elapsed"]] - t0,
       config = list(escala = escala, nsim = nsim, itermax = itermax, NP = NP, T = T,
                     modo = modo, corte_saturacion = corte_saturacion,
                     n_reeval = n_reeval, nsim_reeval = nsim_reeval,
                     motor = MOTOR_RCPP_VERSION))
}
