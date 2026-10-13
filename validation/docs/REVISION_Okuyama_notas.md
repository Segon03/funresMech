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

## Decisiones pendientes del autor

D1 API pública · D2 alcance de distribuciones propias · D3 modo depredador ·
D4 número de versión (propuesta 1.2.0) · D5 alcance del banco de optimizadores ·
D6 presupuesto de cómputo (tras piloto).
