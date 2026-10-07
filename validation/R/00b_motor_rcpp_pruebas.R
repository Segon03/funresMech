# ============================================================
# R/00b_motor_rcpp_pruebas.R
#
# Pruebas del motor Rcpp (fun/motor_okuyama.cpp + fun/motor_rcpp.R)
# ANTES de usarlo en la validacion. Correr desde la carpeta del
# proyecto (.Rproj). Tarda ~5-10 min. No hace falta paralelizar.
#
# Que demuestra cada prueba (resumen para el paper / el paquete):
#   1. El motor del proyecto reproduce BIT A BIT al simulador de
#      Okuyama (Script S5) con la misma semilla (1a, 1b), y
#      DIAGNOSTICO (1c): con las versiones actuales de Rcpp el Script
#      S5 literal reutiliza numeros aleatorios cuando s > 0 (RNGScope
#      anidado en rlnorm0_cpp); el motor del proyecto no tiene ese
#      problema y coincide con una implementacion independiente en R.
#   2. El mismo motor, en modo "funresMech", reproduce BIT A BIT al
#      motor R puro de funresMech v1.0.4 (simulate_trial /
#      negloglik_fixed_z). => el C++ es un reemplazo exacto.
#   3. El corte por saturacion NO cambia la distribucion simulada
#      (solo evita vueltas inutiles del while).
#   4. La "trampa computacional" (h ~ 0, k ~ 0) existe en el
#      simulador original y desaparece con el corte.
#   5. Velocidad: R puro vs C++.
#   6. Anclaje: la NLL en los parametros publicados por Okuyama
#      (Tabla S1) coincide con la NLL implicita en su AIC.
#   7. La reparametrizacion log ida y vuelta es exacta.
#   8. La NLL literal de Okuyama (Script S4) y nll_motor() coinciden.
#
# Salida: tablas/motor_rcpp_pruebas.csv  (una fila por prueba)
# ============================================================

source(file.path("fun", "motor_rcpp.R"))
cat("Compilando/cargando el motor C++ (la primera vez tarda ~1 min)...\n")
cargar_motor_rcpp()
dir.create("tablas", showWarnings = FALSE)

d2 <- as.data.frame(readRDS(file.path("data_clean", "D2.rds")))
d3 <- as.data.frame(readRDS(file.path("data_clean", "D3.rds")))
T_exp <- 24

resultados <- list()
registrar <- function(id, prueba, pasa, detalle) {
  resultados[[length(resultados) + 1]] <<- data.frame(
    id = id, prueba = prueba, resultado = if (isTRUE(pasa)) "PASA" else "FALLA",
    detalle = detalle, stringsAsFactors = FALSE)
  cat(sprintf("[%s] %s: %s\n    %s\n", id, if (isTRUE(pasa)) "PASA " else "FALLA", prueba, detalle))
}

# Parametros de prueba: cubren s = 0 (manipulacion fija) y s > 0,
# densidades chicas y grandes, k chico y grande.
casos <- list(
  list(H = 2,  a = 0.012, h = 0.626, z = 2.5, k = 0.902, s = 0.542),
  list(H = 10, a = 0.02,  h = 0.5,   z = 1.0, k = 1.5,   s = 0),
  list(H = 50, a = 0.001, h = 0.606, z = 3.1, k = 0.037, s = 0.682),
  list(H = 35, a = 0.3,   h = 0.05,  z = 1.2, k = 4.0,   s = 0.2)
)

# ---- 1. Motor (B) == Script S5 de Okuyama, bit a bit ---------------------
# 1a: contra S5 corregido (A', sin RNGScope anidado), todos los casos.
# 1b: contra S5 LITERAL (A), en los casos con s = 0 (donde el RNGScope
#     anidado no interviene).
comparar <- function(simulador, p, n = 5000) {
  set.seed(2026); ref <- simulador(p$H, p$a, p$h, p$z, p$k, p$s, T_exp, n)
  set.seed(2026); nue <- motor_okuyama_cpp(p$H, p$a, p$h, p$z, p$k, p$s, T_exp, n,
                                           modo = 0L, corte_saturacion = FALSE)
  identical(ref, nue)
}
ig_corr <- vapply(casos, function(p) comparar(parasitism_gamma_lnorm_cpp_corr, p), logical(1))
casos_s0 <- Filter(function(p) p$s == 0, casos)
ig_lit  <- vapply(casos_s0, function(p) comparar(parasitism_gamma_lnorm_cpp, p), logical(1))
registrar("1a", "Motor C++ (modo okuyama, sin corte) identico a Okuyama S5 (con RNGScope corregido)",
          all(ig_corr), sprintf("%d/%d casos identicos (5000 ensayos c/u, misma semilla)",
                                sum(ig_corr), length(ig_corr)))
registrar("1b", "Motor C++ identico a Okuyama S5 LITERAL cuando s = 0",
          all(ig_lit), sprintf("%d/%d casos con s = 0 identicos", sum(ig_lit), length(ig_lit)))

# 1c: DIAGNOSTICO del S5 literal cuando s > 0 (RNGScope anidado).
#   (i) demostracion minima: un RNGScope anidado repite el numero anterior
#   (ii) distribucion: referencia en R puro (independiente de Rcpp) vs
#        motor del proyecto vs S5 literal
sim_R_puro <- function(H, a, h, z, k, s, T) {   # Okuyama (2026), ec. del modelo, en R
  t <- 0; hit <- rep(FALSE, H)
  phi <- sqrt(log(1 + (s / h)^2)); mu <- log(h) - phi^2 / 2
  while (t < T) {
    t <- t + rgamma(1, shape = k, scale = 1 / (k * a * H^z))
    if (t > T) break
    hit[floor(runif(1, 0, H)) + 1] <- TRUE
    t <- t + if (s == 0) h else exp(rnorm(1, mu, phi))
  }
  sum(hit)
}
set.seed(1); demo <- demo_rngscope_anidado(5)
repite <- all(demo[, 1] == demo[, 2])
pc <- casos[[1]]   # H = 2, s = 0.542 (parametros de Okuyama para D3)
nR <- 40000
set.seed(101); r_ref <- replicate(nR, sim_R_puro(pc$H, pc$a, pc$h, pc$z, pc$k, pc$s, T_exp))
set.seed(102); r_mot <- motor_okuyama_cpp(pc$H, pc$a, pc$h, pc$z, pc$k, pc$s, T_exp, nR, 0L, TRUE)
set.seed(103); r_lit <- parasitism_gamma_lnorm_cpp(pc$H, pc$a, pc$h, pc$z, pc$k, pc$s, T_exp, nR)
tb <- function(x) tabulate(x + 1L, pc$H + 1L)
p_mot <- chisq.test(rbind(tb(r_ref), tb(r_mot)))$p.value
p_lit <- chisq.test(rbind(tb(r_ref), tb(r_lit)))$p.value
fr <- function(x) paste(sprintf("%.3f", tb(x) / length(x)), collapse = "/")
nll_lit <- replicate(3, nll_gamma_lnorm_S4(c(0.012, 0.626, 0.902, 0.542), d3$dens, d3$par,
                                           T = T_exp, z = 2.5, nsim = 10000))
nll_cor <- replicate(3, nll_gamma_lnorm_S4(c(0.012, 0.626, 0.902, 0.542), d3$dens, d3$par,
                                           T = T_exp, z = 2.5, nsim = 10000,
                                           simulador = parasitism_gamma_lnorm_cpp_corr))
registrar("1c", "DIAGNOSTICO: S5 literal con Rcpp actual reutiliza numeros aleatorios si s > 0; el motor del proyecto no",
          repite && p_mot > 0.001 && p_lit < 1e-6,
          sprintf(paste0("RNGScope anidado repite el numero previo: %s (Rcpp %s). ",
                         "P(0/1/2 parasitados) H=2, params D3 de Okuyama: R puro %s | motor %s (chi2 p=%.2f) | S5 literal %s (chi2 p=%.1e). ",
                         "NLL D3 en los parametros publicados (z=2.5; Okuyama: %.2f): S5 literal %s | S5 corregido %s"),
                  if (repite) "SI" else "NO", as.character(packageVersion("Rcpp")),
                  fr(r_ref), fr(r_mot), p_mot, fr(r_lit), p_lit, (272.890 - 10) / 2,
                  paste(sprintf("%.1f", nll_lit), collapse = ", "),
                  paste(sprintf("%.1f", nll_cor), collapse = ", ")))

# ---- 2. Modo funresMech == motor R puro de funresMech v1.0.4 -------------
if (requireNamespace("funresMech", quietly = TRUE)) {
  iguales2 <- vapply(casos, function(p) {
    n <- 300
    set.seed(7); ref <- vapply(seq_len(n), function(i)
      funresMech:::simulate_trial(p$H, T_exp, p$a, p$h, p$z, p$k, p$s), numeric(1))
    set.seed(7); nue <- motor_okuyama_cpp(p$H, p$a, p$h, p$z, p$k, p$s, T_exp, n,
                                          modo = 1L, corte_saturacion = FALSE)
    identical(as.integer(ref), nue)
  }, logical(1))
  set.seed(11); nll_R <- funresMech:::negloglik_fixed_z(c(0.02, 0.5, 1.2, 0.3), z_fixed = 1.5,
                                                         data_spp = d2, T = T_exp, n_sim = 200)
  set.seed(11); nll_C <- negloglik_fixed_z_rcpp(c(0.02, 0.5, 1.2, 0.3), z_fixed = 1.5,
                                                data_spp = d2, T = T_exp, n_sim = 200,
                                                corte_saturacion = FALSE)
  ok2 <- all(iguales2) && isTRUE(all.equal(unname(nll_R), nll_C, tolerance = 1e-10))
  registrar(2, "Motor C++ (modo funresMech) identico al motor R de funresMech v1.0.4", ok2,
            sprintf("simulate_trial: %d/%d casos identicos; negloglik_fixed_z D2: R=%.6f C++=%.6f",
                    sum(iguales2), length(iguales2), nll_R, nll_C))
} else {
  registrar(2, "Motor C++ (modo funresMech) identico al motor R de funresMech v1.0.4", NA,
            "funresMech no esta instalado: prueba omitida")
}

# ---- 3. El corte por saturacion no cambia la distribucion ---------------
# Casos donde la saturacion es frecuente (H chico o a alto).
casos_sat <- list(list(H = 2,  a = 0.012, h = 0.626, z = 2.5, k = 0.902, s = 0.542),
                  list(H = 10, a = 0.3,   h = 0.05,  z = 1.2, k = 0.2,   s = 0.2),
                  list(H = 6,  a = 0.02,  h = 0.3,   z = 3.0, k = 0.05,  s = 0.6))
n3 <- 200000
p3 <- vapply(casos_sat, function(p) {
  set.seed(31); sin_corte <- motor_okuyama_cpp(p$H, p$a, p$h, p$z, p$k, p$s, T_exp, n3, 0L, FALSE)
  set.seed(32); con_corte <- motor_okuyama_cpp(p$H, p$a, p$h, p$z, p$k, p$s, T_exp, n3, 0L, TRUE)
  tab <- rbind(tabulate(sin_corte + 1L, p$H + 1L), tabulate(con_corte + 1L, p$H + 1L))
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  if (ncol(tab) < 2) return(1)
  suppressWarnings(chisq.test(tab)$p.value)
}, numeric(1))
# y sobre la NLL completa de D2 (misma escala que usa el ajuste)
nll_sin <- replicate(10, nll_motor(c(0.001, 0.606, 0.037, 0.682), 3.1, d2, T_exp, 10000,
                                   corte_saturacion = FALSE))
nll_con <- replicate(10, nll_motor(c(0.001, 0.606, 0.037, 0.682), 3.1, d2, T_exp, 10000,
                                   corte_saturacion = TRUE))
tt <- t.test(nll_sin, nll_con)
registrar(3, "Corte por saturacion no altera la distribucion simulada",
          all(p3 > 0.001) && tt$p.value > 0.001,
          sprintf("chi2 p-valores (2e5 ensayos c/u): %s; NLL D2 media sin/con corte: %.3f / %.3f (t-test p=%.2f)",
                  paste(sprintf("%.3f", p3), collapse = ", "), mean(nll_sin), mean(nll_con), tt$p.value))

# ---- 4. Trampa computacional ----------------------------------------------
# Punto tipico que DEoptim visita a z alto: h y k casi 0 (piso de 04d/04f).
# Se cuentan las vueltas del while en 10 ensayos a H = 50 (una
# evaluacion de la NLL necesita 10000 ensayos x 12 densidades).
trampa <- list(a = 1, h = 1e-6, z = 3, k = 1e-6, s = 0)
set.seed(4); v_sin <- contar_vueltas_cpp(50, trampa$a, trampa$h, trampa$z, trampa$k, trampa$s,
                                         T_exp, 10, FALSE, max_vueltas = 5e8)
set.seed(4); v_con <- contar_vueltas_cpp(50, trampa$a, trampa$h, trampa$z, trampa$k, trampa$s,
                                         T_exp, 10, TRUE)
t_con <- system.time(nll_motor(unlist(trampa[c("a", "h", "k", "s")]), trampa$z, d2, T_exp, 10000))[["elapsed"]]
registrar(4, "Trampa h~0/k~0: el simulador original se 'cuelga'; con corte, no",
          (is.infinite(v_sin) || v_sin > 1e7) && v_con < 1e4 && t_con < 60,
          sprintf("vueltas del while en 10 ensayos (H=50): sin corte %s; con corte %.0f. NLL completa D2 (10000 ensayos) con corte: %.2f s",
                  if (is.infinite(v_sin)) ">5e8 (se corto la cuenta)" else format(v_sin, big.mark = ","),
                  v_con, t_con))

# ---- 5. Velocidad -----------------------------------------------------------
par_tip <- c(a = 0.02, h = 0.4, k = 1, s = 0.3)
t_cpp <- system.time(for (i in 1:3) nll_motor(par_tip, 1.5, d2, T_exp, 10000))[["elapsed"]] / 3
if (requireNamespace("funresMech", quietly = TRUE)) {
  n_r <- 300
  t_R <- system.time(funresMech:::negloglik_fixed_z(c(0.02, 0.4, 1, 0.3), 1.5, d2, T_exp, n_r))[["elapsed"]]
  t_R10k <- t_R * 10000 / n_r
  det5 <- sprintf("1 NLL de D2 con 10000 ensayos: C++ %.2f s; R puro ~%.0f s (extrapolado de %d ensayos) -> %.0fx mas rapido",
                  t_cpp, t_R10k, n_r, t_R10k / t_cpp)
} else det5 <- sprintf("1 NLL de D2 con 10000 ensayos: C++ %.2f s", t_cpp)
horas_tarea <- t_cpp * 40 * 201 / 3600
registrar(5, "Velocidad del motor C++", t_cpp < 30,
          paste0(det5, sprintf(". => 1 punto del perfil (NP=40, itermax=200) ~%.1f h por nucleo", horas_tarea)))

# ---- 6. Anclaje contra la Tabla S1 de Okuyama -----------------------------
# NLL implicita en el AIC publicado (5 parametros): NLL = (AIC - 10)/2.
# 'a' esta redondeado en la tabla (D2: 0.001 -> 1 cifra significativa),
# asi que ademas se busca el mejor 'a' cerca del publicado dejando fijos
# h, k, s (mismas semillas para todos los 'a': numeros aleatorios comunes).
anclas <- list(
  D2 = list(d = d2, a = 0.001, h = 0.606, z = 3.1, k = 0.037, s = 0.682, AIC = 336.589),
  D3 = list(d = d3, a = 0.012, h = 0.626, z = 2.5, k = 0.902, s = 0.542, AIC = 272.890))
det6 <- character(0); ok6 <- TRUE
for (lab in names(anclas)) {
  A <- anclas[[lab]]
  nll_pub <- (A$AIC - 10) / 2
  nll_en_pub <- mean(replicate(5, nll_motor(c(A$a, A$h, A$k, A$s), A$z, A$d, T_exp, 10000)))
  f_a <- function(la) { set.seed(66); nll_motor(c(exp(la), A$h, A$k, A$s), A$z, A$d, T_exp, 10000) }
  opt <- optimize(f_a, log(A$a) + c(-log(2), log(2)))    # a entre a/2 y 2a
  a_opt <- exp(opt$minimum)
  nll_a_opt <- mean(replicate(5, nll_motor(c(a_opt, A$h, A$k, A$s), A$z, A$d, T_exp, 10000)))
  ok6 <- ok6 && abs(nll_a_opt - nll_pub) < 1.5
  det6 <- c(det6, sprintf("%s: NLL Okuyama (del AIC) = %.2f; motor en parametros publicados = %.2f; con 'a' re-ajustado (a=%.5f, publicado %.3f) = %.2f",
                          lab, nll_pub, nll_en_pub, a_opt, A$a, nll_a_opt))
}
registrar(6, "Anclaje: NLL en los parametros de la Tabla S1 ~ NLL implicita en el AIC de Okuyama",
          ok6, paste(det6, collapse = " | "))

# ---- 7. Reparametrizacion log: ida y vuelta --------------------------------
pn <- c(a = 0.00123, h = 0.606, k = 0.037, s = 0.682)
H_ref <- h_ref_de(d2)
vuelta <- theta_a_nat(nat_a_theta(pn, 3.1, H_ref), 3.1, H_ref)
registrar(7, "Reparametrizacion log: ida y vuelta exacta",
          isTRUE(all.equal(vuelta, pn, tolerance = 1e-12)),
          sprintf("H_ref(D2) = %.3f; max error relativo = %.1e", H_ref, max(abs(vuelta / pn - 1))))

# ---- 8. NLL de Okuyama (Script S4) == nll_motor() ------------------------
# S4 literal, con el simulador S5 corregido (A'), s.sqrt = TRUE como en S3.
p8 <- c(0.012, 0.626, 0.902, sqrt(0.542))
set.seed(88); nll_S4 <- nll_gamma_lnorm_S4(p8, d3$dens, d3$par, T = T_exp, z = 2.5, nsim = 3000,
                                           s.sqrt = TRUE, simulador = parasitism_gamma_lnorm_cpp_corr)
set.seed(88); nll_B  <- nll_motor(c(p8[1], p8[2], p8[3], p8[4]^2), 2.5, d3, T_exp, 3000,
                                  modo = "okuyama", corte_saturacion = FALSE)
registrar(8, "NLL de Okuyama (Script S4, simulador S5 corregido) == nll_motor() con la misma semilla",
          isTRUE(all.equal(nll_S4, nll_B, tolerance = 1e-10)),
          sprintf("D3, z=2.5: S4 = %.8f; nll_motor = %.8f", nll_S4, nll_B))

# ---- Resumen ---------------------------------------------------------------
tabla <- do.call(rbind, resultados)
tabla$motor <- MOTOR_RCPP_VERSION
tabla$R_version <- R.version.string
tabla$fecha <- format(Sys.time(), "%Y-%m-%d %H:%M")
write.csv(tabla, file.path("tablas", "motor_rcpp_pruebas.csv"), row.names = FALSE)
cat("\n==== RESUMEN ====\n")
print(tabla[, c("id", "resultado", "prueba")], row.names = FALSE)
if (all(tabla$resultado == "PASA")) {
  cat("\nTODAS LAS PRUEBAS PASAN. Se puede correr R/04h_ensayo2_rcpp_log.R\n")
} else {
  cat("\nHAY PRUEBAS QUE FALLAN: revisar antes de correr 04h.\n")
}
