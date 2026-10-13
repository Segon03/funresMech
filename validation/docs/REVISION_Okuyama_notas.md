# Revisión en respuesta a Okuyama (2026-10-13): notas de trabajo

Plan completo: `PLAN_respuesta_Okuyama.md` (carpeta `Okuyama_rev` y proyecto).

## Flujo de trabajo

* Línea base congelada: funresMech 1.1.0 tal como se publicó (CRAN 2026-10-06,
  GitHub commit 52a779d). Copia intacta en `Okuyama_rev/funresMech_1.1.0_OR`.
* Copia de trabajo con historial git: `Okuyama_rev/funresMech_dev`,
  rama `revision-okuyama`, etiqueta `baseline-1.1.0` en el commit base.
* El código se edita y se prueba en el entorno de Claude (R 4.3.3, Linux) y
  cada paso se aplica a `funresMech_dev` como commit. Las corridas largas
  (validación, banco de optimizadores) se ejecutan en la PC del autor (Windows,
  23 núcleos). Nada se sube a GitHub ni a CRAN hasta la Fase 9.

## Puerta de regresión (Fase 0.3)

`tests/testthat/test-regression.R` + `tests/testthat/fixtures/golden_1.1.0.rds`.

* Ancla contra el código de Okuyama (2026): D3, z = 2.5, par = (0.012, 0.626,
  0.902, 0.542), 3000 ensayos simulados, semilla 88, sin corte por saturación
  → NLL = 132.04994217 (reproducido con el paquete 1.1.0 instalado).
* Conteos simulados exactos del motor (5 casos × con/sin corte) y NLL exactas
  (5 casos D2/D3) de la versión 1.1.0.
* El modelo por defecto (búsqueda Gamma, manipulación Lognormal) debe
  reproducir estos valores bit a bit en todas las versiones siguientes.
* El valor 7944.851530 (D2) de la validación compara con la implementación en
  R puro de 1.0.4, cuya `s` era el desvío en escala log; no es reproducible con
  la definición actual de `s` y por eso no forma parte de la puerta.

## Inconsistencias registradas (a resolver en las fases indicadas)

| ID | Inconsistencia | Fase |
|---|---|---|
| H3 | El IC de z y z_hat usan `nll` (mínimo ruidoso del optimizador), no `nll_reeval`. | 3.1 |
| H4 | `n_reeval` = 3 por defecto en el código; el S1 Text dice que el NLL se reevaluó cinco veces. | 3.1 |
| H6 | El paquete compara AIC de z libre vs z = 1; el manuscrito dice que compara con beta-binomial. | 3.11 |
| H8 | Grilla de z por defecto en la app: 0.5–5 paso 0.25; manuscrito: 0.5–6.0 paso 0.1. | 3.8 / 8 |
| H1–H2 | Cotas absolutas (h ≤ 5, λ_ref en [1e-3, 1e3]) y T rotulado en horas: dependencia de la unidad de tiempo. | 2 |
| H14 | `reltol` de DEoptim no tenía efecto en 1.1.0: DEoptim lo evalúa cada `steptol` generaciones y `steptol` vale `itermax` por defecto, así que nunca había parada temprana. | resuelto en `funres_control()` (reltol y steptol explícitos; por defecto sin parada temprana) |
| H15 | Las cotas del modelo beta-binomial del screening (a, h <= 5) también son absolutas y dependen de la unidad de T. | 2 |

## Decisiones del autor (2026-10-13)

* D1: se exporta una API pública además de `run_app()`.
* D2: distribuciones en tres niveles: familias en C++, muestreadores en R y
  muestreadores en C++ (puntero externo).
* D3: modo depredador: verificar primero. Verificado: el Apéndice S1 y el
  Script S2 de Okuyama (2026) definen el proceso (presa consumida retirada,
  dH/dt = -P a H^z / (1 + a h H^z)); el propio apéndice indica que la elección
  depende del sistema (un parasitoide que no revisita hospedadores
  parasitados equivale al proceso depredador). Se implementa como
  `process = "predator"`.
* D4: versión 1.2.0.
* D5 (optimizadores) y D6 (presupuesto de cómputo): pendientes del piloto.

## Fase 1 (estado)

* Motor generalizado (`src/motor.cpp`): familias de búsqueda Gamma, Lognormal,
  Weibull, R y C++ del usuario (media 1/lambda, CV^2 = 1/k); manipulación
  Lognormal, Gamma, Weibull, constante, R y C++ (media h, desvío s); proceso
  parasitoide o depredador. Ruta rápida del modelo por defecto idéntica bit
  a bit a 1.1.0.
* API: `funres_model()`, `funres_sampler()`, `test_sampler()`,
  `simulate_response()`, `funres_control()`, `fit_funres()`.
* Parametrización adimensional (tiempo en unidades de T) y cotas por defecto
  independientes de la unidad (adelanto de la Fase 2).
* Perfil e IC sobre la media de las reevaluaciones (`nll`), con su error
  estándar (`nll_se`); el mínimo del optimizador queda como `nll_opt`
  (adelanto de la Fase 3.1).
