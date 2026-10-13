# ============================================================================
# 10_piloto_costos.R
#
# Piloto de costos de cómputo (Plan de respuesta a Okuyama, tarea 5.1 / D6).
# Mide en ESTA PC:
#   1. el tiempo de una evaluación de la verosimilitud simulada según n_sim y
#      el diseño experimental (D2, D3 y un diseño grande);
#   2. el tiempo de una optimización completa en un punto del perfil
#      (preset "exploratory": 2000 simulaciones, 50 generaciones) con cada
#      método; el costo del preset "standard" (5000, 100) se proyecta x 5;
#   3. la proyección del tiempo de la comparación de optimizadores (Fase 4) y
#      de la validación ampliada (Fase 5) según el número de núcleos.
#
# Uso (R en Windows, desde la carpeta funresMech_dev):
#   devtools::install()            # instala la versión 1.2.0 de desarrollo
#   devtools::test()               # todas las pruebas deben pasar
#   source("validation/R/10_piloto_costos.R")
#
# Duración aproximada: 10 a 20 minutos (un solo núcleo).
# Resultado: validation/tablas/piloto_costos.csv y un resumen en pantalla.
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.2.0")
for (p in c("cmaes", "nloptr"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

dir.create("validation/tablas", showWarnings = FALSE, recursive = TRUE)
d2 <- read.csv("tests/testthat/fixtures/D2.csv")
d3 <- read.csv("tests/testthat/fixtures/D3.csv")
grande <- data.frame(dens = rep(c(5, 10, 20, 30, 40, 50, 60, 80, 100, 120), each = 10))
grande$par <- 0L   # solo se mide el tiempo; el valor de la NLL no importa

m <- funres_model()
p_d3 <- c(a = 0.012, h = 0.626, z = 2.5, k = 0.902, s = 0.542)
p_d2 <- c(a = 0.001, h = 0.606, z = 3.1, k = 0.037, s = 0.682)
T_exp <- 24

# ---- 1. Evaluación de la NLL ------------------------------------------------
cat("\n1. Tiempo por evaluación de la NLL\n")
eval_time <- function(datos, par, n_sim, rep = 5) {
  set.seed(1)
  t <- system.time(for (i in seq_len(rep))
    funresMech:::.nll_sim(par, m, datos, T_exp, n_sim))[["elapsed"]]
  t / rep
}
filas <- list()
for (dis in c("D2", "D3", "grande")) {
  datos <- switch(dis, D2 = d2, D3 = d3, grande = grande)
  par <- if (dis == "D2") p_d2 else p_d3
  for (ns in c(1000, 2000, 5000, 10000)) {
    s <- eval_time(datos, par, ns)
    filas[[length(filas) + 1]] <- data.frame(medida = "nll_eval", diseno = dis,
                                             metodo = NA, n_sim = ns, segundos = s)
    cat(sprintf("  %-7s n_sim = %5d: %.3f s\n", dis, ns, s))
  }
}

# ---- 2. Optimización en un punto (preset exploratory) -----------------------
cat("\n2. Optimización en un punto del perfil (D3, z = 2.5, preset exploratory)\n")
for (meth in c("DEoptim", "CMAES", "CRS2")) {
  ctrl <- funres_control("exploratory", method = meth)
  set.seed(2)
  t <- system.time(
    pt <- funresMech:::.fit_point(2.5, d3, T_exp, m, ctrl,
                                  funresMech:::.make_bounds(list(), T_exp),
                                  funresMech:::.h_ref(d3)))[["elapsed"]]
  filas[[length(filas) + 1]] <- data.frame(medida = "punto_exploratory", diseno = "D3",
                                           metodo = meth, n_sim = ctrl$n_sim, segundos = t)
  cat(sprintf("  %-8s %.0f s  (NLL reevaluada %.2f)\n", meth, t, pt$nll_reeval))
}
tab <- do.call(rbind, filas)
tab$nucleos_pc <- parallel::detectCores()
tab$fecha <- format(Sys.time(), "%Y-%m-%d %H:%M")
write.csv(tab, "validation/tablas/piloto_costos.csv", row.names = FALSE)

# ---- 3. Proyección ------------------------------------------------------------
nuc <- max(1, parallel::detectCores() - 1)
t_pt <- 5 * mean(tab$segundos[tab$medida == "punto_exploratory"])   # standard ~ 5 x exploratory
cat(sprintf("\n3. Proyección con %d núcleos (punto standard ~ %.0f s)\n", nuc, t_pt))
proy <- function(lbl, n_puntos) cat(sprintf("  %-62s %7.1f h\n", lbl, n_puntos * t_pt / nuc / 3600))
proy("Fase 4: 9 escenarios x 10 datos x 3 métodos x 3 puntos de z", 9 * 10 * 3 * 3)
proy("Fase 4: método conjunto vs perfil, 9 x 10 x (1 + ~12 puntos)", 9 * 10 * 13)
proy("Fase 5 nivel 1: 150 combinaciones x 10 datos x ~13 puntos (conjunto)", 150 * 10 * 13)
proy("Fase 5 nivel 2: 30 combinaciones x 10 datos x 46 puntos (perfil)", 30 * 10 * 46)
cat("\nEnviar a Claude el archivo validation/tablas/piloto_costos.csv.\n")
