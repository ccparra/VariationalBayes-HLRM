# VariationalBayes-HLRM

Implementación en R de modelos bayesianos de regresión lineal:

- `LRM`: regresión lineal para una población.
- `HLRM`: regresión jerárquica con coeficientes y varianza residual por grupo.
- `CHLRM`: mezcla finita jerárquica para agrupar unidades con perfiles de regresión compartidos.

Los tres modelos incluyen MCMC, CAVI y SVI. SVI usa minibatches de
observaciones en LRM y de grupos completos en HLRM/CHLRM.

Este repositorio contiene exclusivamente el paquete estadístico reusable. La
tesis, los datos aplicados y la evaluación E1--E5 viven en el repositorio
privado de investigación `VariationalBayes-HLRM-research`, que consume esta API
sin copiar sus kernels.

## Estructura

```text
R/             implementación reusable de LRM, HLRM y CHLRM
man/           referencia de la API
tests/         pruebas automáticas del núcleo estadístico
dev/           validación para desarrollo
docs/          contrato de API y roadmap vigente
```

El contrato estable de métodos y resultados está en
[contrato de la API](docs/API.md) y el trabajo pendiente se mantiene en un
único [roadmap del paquete](docs/ROADMAP.md).

## Requisitos

El código fue validado con R 4.6.1 y utiliza los paquetes:

```r
c("R6", "mvtnorm", "LaplacesDemon", "gtools", "posterior", "loo",
  "testthat", "clue")
```

Instalación básica:

```sh
Rscript -e 'install.packages(c("R6","mvtnorm","LaplacesDemon","gtools","posterior","loo","testthat","clue"), repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

## Uso

Ejecute los comandos desde la raíz del repositorio.

```sh
Rscript dev/run-validation.R
R CMD build .
R CMD check --no-manual VariationalBayesHLRM_0.1.0.tar.gz
```

Cada ejecución escribe en su propio subdirectorio `output/`, ignorado por Git.
`dev/run-validation.R` conserva localmente `tests/output/latest-summary.txt`,
`results.csv` y `session-info.txt`; no guarda draws ni objetos de ajuste.

`tests/` no es material legado: contiene regresiones automatizadas de las
identidades algebraicas, las actualizaciones CAVI/SVI, el relabeling y los
diagnósticos. Debe permanecer separado de los experimentos porque sus pruebas
son rápidas, deterministas y se ejecutan después de cada cambio del núcleo.

Ejemplo mínimo:

```r
library(VariationalBayesHLRM)

X <- cbind(intercept = 1, predictor = seq(-2, 2, length.out = 40))
data <- list(X = X, y = drop(X %*% c(1, 2)) + rnorm(40, sd = .5))
prior <- lrm_prior(
  coefficient_mean = c(0, 0), coefficient_covariance = 10 * diag(2),
  residual_degrees_of_freedom = 4, residual_variance_scale = 1)

fit <- fit_model(
  lrm_model(2), data, prior, method = "cavi",
  initialization = initialization_ols(),
  control = cavi_control(max_iterations = 2000), seed = 10,
  verbose = FALSE)
summary <- summary(fit)
coef(fit)
nobs(fit)
draws <- posterior::as_draws(fit, n_draws = 500, seed = 11)
```

Los resultados públicos tienen contratos comunes entre modelos y métodos:

- `fit_model()` devuelve `fit-result-2`, con `parameters`, `objective`,
  `diagnostics`, `timing` y `metadata`; CHLRM añade `latent_variables`;
- `summarize_posterior()` devuelve `posterior-summary-1`;
- `predict()` y `posterior_predict()` devuelven `posterior-prediction-2`;
- `evaluate_in_sample_fit()` devuelve `in-sample-fit-1`;
- `posterior_predictive_check()` devuelve `posterior-predictive-check-2`;
- `predictive_criteria()` devuelve `predictive-criteria-2`.

Los nombres de cada cantidad son canónicos, pero cada modelo devuelve sólo lo
que matemáticamente posee. LRM no expone componentes ni variables latentes;
HLRM no expone parámetros de mezcla; CHLRM sí añade asignaciones o
responsabilidades. En el resultado de ajuste, MCMC guarda draws y VI/SVI los
parámetros de cada factor; ELBO, log-likelihood y log-joint viven por separado
en `fit$objective`.

El ajuste implementa `print()`, `summary()`, `coef()`, `fitted()`,
`residuals()` y `nobs()`. CHLRM añade, sólo donde corresponden,
`components()`, `responsibilities()`, `classify_groups()`,
`classification_uncertainty()` y `co_clustering()`.
HLRM ofrece `population_parameters()`, `group_parameters()` y
`group_deviations()`; esta última devuelve la media posterior de
$\beta_j-\beta$, sin fabricar intervalos a partir de marginales. Para tablas
y gráficos, `tidy_parameters()` produce una fila por elemento de parámetro.

Antes de ajustar, `estimate_fit_resources()` aplica el mismo cálculo preventivo
que usa `fit_model()`. El límite `max_memory_mb` está disponible en
`mcmc_control()`, `cavi_control()` y `svi_control()`, incluso para una sola
cadena o inicio:

```r
control <- mcmc_control(
  warmup = 1000, post_warmup_draws = 2000, chains = 4,
  stopping = mcmc_stopping(), max_memory_mb = 2048)
resources <- estimate_fit_resources(hlrm_model(2), grouped_data, control)
print(resources)
```

La estimación separa el pico esperado durante el ajuste, el estado retenido y
el almacenamiento potencial de una evaluación posterior. Incluye holgura para
objetos de R, pero no es una medición ni un límite duro del sistema operativo.

### MCMC multichain

`fit_model()` orquesta varias cadenas con el mismo contrato para LRM, HLRM y
CHLRM. Las ejecuta secuencialmente para no sobresuscribir CPU/RAM, exige
`thinning = 1`, registra semillas y tiempos, y usa `posterior` para R-hat
rank-normalizado, bulk ESS, tail ESS y MCSE. En CHLRM aprende una referencia
común y alinea todas las cadenas antes del diagnóstico.

```r
fit_mc <- fit_model(
  hlrm_model(2), grouped_data, prior, method = "mcmc",
  initialization = initialization_from_prior(),
  control = mcmc_control(
    warmup = 1000, post_warmup_draws = 2000, chains = 4,
    stopping = mcmc_stopping(), max_memory_mb = 1024),
  seed = 2026, verbose = FALSE
)
fit_mc$diagnostics$metrics$variables
draws <- posterior_draws(fit_mc)
```

Los umbrales por defecto son R-hat <= 1.01 y bulk/tail ESS >= 1000. El
presupuesto de memoria es una estimación conservadora previa, no un límite duro
del proceso.

Internamente se extiende por bloques y se recalculan los diagnósticos hasta
alcanzar sus umbrales o el máximo declarado. Un ajuste MCMC completo puede
continuarse sin mutar el original y guardarse como estado plano:

```r
extended <- resume_mcmc(fit_mc, additional_samples = 1000)
save_mcmc_checkpoint(extended, "fit-mcmc-v1.rds")
restored <- read_mcmc_checkpoint("fit-mcmc-v1.rds")
```

En varias cadenas se conserva el estado terminal y RNG de cada cadena, y tras
continuar se recalculan el relabeling CHLRM y los diagnósticos conjuntos. El
guardado rehúsa sobrescribir una ruta existente para favorecer checkpoints
versionados y recuperables.

### VI/SVI multistart

`fit_model()` separa la especificación de inicialización de la aleatoriedad del
algoritmo. CAVI y SVI usan los mismos nombres, reglas de parada y contrato de
salida; con `starts > 1` conservan todos los intentos y seleccionan por ELBO
finito bajo la regla de convergencia declarada.

```r
fit_svi <- fit_model(
  chlrm_model(3, 2), grouped_data, prior, method = "svi",
  initialization = initialization_from_prior(),
  control = svi_control(
    max_iterations = 50000, minibatch_size = 32,
    learning_rate = power_schedule(offset = 100, exponent = .7),
    starts = 5, require_converged = TRUE, max_memory_mb = 2048),
  seed = 2029, verbose = FALSE
)
fit_svi$metadata$start_summary
```

Se conservan parámetros, ELBO, iteraciones, convergencia, tiempo, semilla y
error de cada intento. Por defecto sólo son elegibles los intentos convergentes
con ELBO terminal finito; entre ellos se selecciona el mayor ELBO. Si ninguno
cumple, el estado es `no_valid_start` y no se elige silenciosamente un ajuste.
LRM admite multistart CAVI/SVI; su SVI muestrea observaciones, mientras
HLRM/CHLRM muestrean grupos completos.

## Decisiones estadísticas importantes

### Contratos numéricos y de prior

Los constructores de LRM, HLRM y CHLRM validan datos, dimensiones y
hiperparámetros antes de iniciar MCMC, VI o SVI. Las covarianzas, precisiones y
escalas deben ser simétricas definidas positivas. `solve_spd()` y
`logdet_spd()` usan Cholesky; no se reparan matrices indefinidas. La función
`symmetrize_numeric()` sólo elimina asimetrías del orden del error numérico y
falla si la discrepancia puede representar un error de fórmula.

Para una prior inverse-Wishart de dimensión `p` se exige `n0 > p - 1`, que es
la condición de distribución propia. Esto no implica que exista su media; si
se reporta ese momento deben imponerse además los grados de libertad que exige
la parametrización concreta. Los generadores empirical-Bayes regularizan los
ajustes OLS singulares o perfectos para producir escalas SPD, pero siguen
reutilizando los datos y no sustituyen una prior sustantivamente justificada.

### Ejemplo controlado de CHLRM

Este ejemplo pequeño comprueba la conexión entre la simulación genérica y una
likelihood marginal. La geometría y la prior concretas del review no son
defaults del paquete:

```r
library(VariationalBayesHLRM)
model <- chlrm_model(component_count = 2, coefficient_count = 3)
d <- simulate_data(
  conditional_dgp(model, list(
    component_coefficients = rbind(c(0, 1, -.5), c(2, -1, .5)),
    component_variances = c(.7, 1.2),
    mixture_weights = c(.4, .6))),
  simulation_design("independent_normal"),
  simulation_sample(group_count = 12, observations_per_group = 20,
                    test_observations_per_group = 10),
  seed = 80, max_output_mb = 64)

# Tras ajustar con fit_model(...), la operación pública calcula WAIC y,
# opcionalmente, conserva la log-verosimilitud pointwise usada:
criteria <- predictive_criteria(fit, return_log_likelihood = TRUE,
                                max_output_mb = 64)
ll <- criteria$log_likelihood$group
```

`log_likelihood()` devuelve el contrato versionado
`log-likelihood-1` con los valores pointwise y el estimando explícito. Las
funciones especializadas `loglik_*` permanecen como kernels numéricos;
`predictive_criteria()` y WAIC usan el contrato común para no confundir observaciones, grupos ni
marginalización de variables latentes.

Los valores exactos de una prior experimental, geometría de simulación,
escenarios y número de réplicas pertenecen al proyecto que consume el paquete;
no son defaults generales recomendados para sus usuarios.

### Evaluación: responsabilidades separadas

Después de `fit_model()`, sin depender de `summarize_posterior()`, el ajuste ofrece:

```r
fit_description <- evaluate_in_sample_fit(fit)
ppc <- posterior_predictive_check(
  fit,
  n_draws = 100, seed = 93, return_yrep = FALSE
)
criteria <- predictive_criteria(
  fit,
  n_draws = 100, seed = 94, predictive_target = "new_groups",
  return_log_likelihood = TRUE
)
```

- `evaluate_in_sample_fit()` calcula MSE/R² descriptivos sobre los datos usados
  para ajustar; no afirma capacidad fuera de muestra.
- `posterior_predictive_check()` evalúa adecuación generativa mediante réplicas
  posteriores. Solo materializa `y_rep` cuando se solicita; informa MCSE IID y
  MCSE corregido por ESS cuando los draws provienen de MCMC.
- `predictive_criteria()` calcula WAIC bajo un estimando explícito. Con
  `return_log_likelihood = TRUE` conserva también la
  matriz pointwise `log-likelihood-1` usada en el cálculo. No diagnostica
  convergencia.
- `predict()` y `posterior_predict()` producen predicciones para diseños dados;
  la segunda incluye ruido observacional.

| Modelo | Objetivo público | Likelihood elegida internamente |
|---|---|---|
| LRM | `observed_data`; unidad `observation` o `group` | condicional en coeficientes y varianza; agrupar sólo cambia la unidad |
| HLRM | `observed_groups` | condicional en efectos y varianzas locales |
| HLRM | `new_groups` | marginal, integrando efectos y varianzas locales |
| CHLRM | `new_groups` | marginal por grupo, integrando la asignación de mezcla |

En CHLRM, por ejemplo:

```r
criteria <- predictive_criteria(
  fit,
  n_draws = 100, seed = 93, return_log_likelihood = TRUE,
  max_output_mb = 64
)
criteria$criteria$waic_group_marginal
log_lik_draw_unit <- as_log_lik_matrix(criteria$log_likelihood$group)
```

DIC fue retirado de la API, los tests y los reportes: no se devuelve como campo
vacío ni se conserva un argumento de compatibilidad. HLRM marginal usa
integración gaussiana y cuadratura determinista, con tolerancia y presupuesto
configurables. En LRM, una unidad agrupada se
reporta como `waic_group`, nunca como marginal, porque no existe un efecto de
grupo que integrar.

La matriz pointwise no se guarda automáticamente en el objeto ajustado: puede
ser grande y no forma parte del estado de inferencia. Los runners activos la
solicitan y persisten en `evaluation-with-log-likelihood.rds`. En el ajuste,
MCMC conserva además su traza escalar de log-likelihood total; VI/SVI conserva
ELBO. Esas trazas cumplen funciones algorítmicas distintas y no deben
compararse como si fueran el mismo criterio.

**Predictiva posterior:** después de `fit_model()`, sin requerir
`summarize_posterior()`, los nueve motores comparten esta API:

```r
# HLRM/CHLRM: 1 identifica un grupo usado en el ajuste; NA crea uno nuevo.
mean_prediction <- predict(
  fit,
  X_new = list(X_known, X_new_group), group_ids = c(1, NA),
  n_draws = 500, level = 0.95, seed = 10
)
response_prediction <- posterior_predict(
  fit,
  X_new = list(X_known, X_new_group), group_ids = c(1, NA),
  n_draws = 500, level = 0.95, seed = 10, return_draws = TRUE
)
ppc <- pp_check_from_prediction(y_new, response_prediction)
```

En LRM, `X_new` es una matriz y no se usa `group_ids`. `predict()` integra la
incertidumbre posterior/q de la media; sus límites son intervalos creíbles
puntuales. `posterior_predict()` añade el ruido residual; sus límites son
intervalos predictivos puntuales y `draws` guarda una matriz observaciones ×
draws por grupo. No son bandas simultáneas.

Para grupos HLRM nuevos se genera un efecto y una varianza locales por grupo y
draw desde la jerarquía global. Para grupos CHLRM nuevos se genera una etiqueta
por grupo y draw desde `omega`; para grupos observados se usa la incertidumbre
posterior de su asignación. No se usa una clasificación rígida `z_pred`.

MCMC aproxima la posterior predictiva con draws de la cadena; VI/SVI integran
bajo q. Esto incluye ruido residual y variación de parámetros, pero no garantiza
calibración de q. El PPC usa la misma primitiva generativa en modo streaming;
`pp_check_from_prediction()` permite reutilizar draws ya almacenados sin volver
a simular. WAIC y los scores predictivos usan log-verosimilitudes y no se calculan desde
`y_rep`.

Los PPC procesan un dataset a la vez. LRM/HLRM procesan la likelihood por bloques
de grupo, sin retener una matriz observaciones × draws para todos los grupos a
la vez. Los presupuestos controlan arrays principales, no toda la RAM de R.

### Ejecución reproducible y parada

En las clases LRM, HLRM y CHLRM, `run_settings$seed` controla
**inicialización y ajuste**, conservando el RNG del llamador. Sin semilla se
mantiene el RNG global de R. El argumento opcional `fit(seed = ...)` cambia sólo
la secuencia del ajuste: no reconstruye una inicialización ya realizada.

| Control | VI/CAVI | SVI |
|---|---|---|
| `epsilon` por defecto | `1e-8` | `1e-4` |
| `elbo_every` | 1, obligatorio | 50 |
| `elbo_window` | 1, obligatorio | 5 |
| `patience` | 5 cambios consecutivos | 5 checkpoints consecutivos elegibles |
| `min_iter` | 1 | 1000 |

El cambio relativo es `abs(actual - anterior) / (1 + abs(anterior))`; en SVI
se comparan medias móviles. `epsilon = NULL` o `0` desactiva la parada por
tolerancia. `iterations` sigue siendo un máximo explícito; cada análisis debe
fijarlo externamente según su protocolo. Una evaluación final fuera del
calendario se registra pero no suma paciencia.
El ELBO se calcula sobre todos los grupos con los factores variacionales
actuales; no es el ELBO de un minibatch ni una certificación de fidelidad posterior.

`fit_model()` devuelve `diagnostics`: razón de
parada, iteraciones, convergencia por tolerancia y controles efectivos. Alcanzar
el máximo no se marca como convergencia. Un ELBO no finito provoca un error.

Cada objeto representa un solo ajuste. CAVI/SVI conserva todos los resultados
multistart, pero todavía no reanuda un loop interrumpido a mitad de iteración;
MCMC sí admite continuación y checkpoint tipados. El resumen posterior, el
ajuste descriptivo, el PPC y los
criterios predictivos exigen un ajuste exitoso, pero son operaciones
independientes en los tres modelos.

### Criterios generales

- Los modelos trabajan en la escala suministrada o usan explícitamente
  `fit_model(..., scaling = data_scaling(...))`; el mapa aprendido en
  entrenamiento es auditable y los resultados regresan a la escala original.
- CHLRM es invariante ante permutaciones de etiquetas. MCMC ofrece relabeling
  conjunto y reporta partición de Dahl y matriz de similitud posterior.
- Los tres modelos permiten WAIC con unidad finca, especificando grupos en LRM
  e integración marginal en HLRM. No es intercambiable con criterios
  condicionales a efectos locales. La evaluación por validación cruzada
  dejando fincas completas fuera sigue siendo una comprobación independiente.
- DIC fue retirado; la comparación usa WAIC con unidad explícita, LOO/CV/LOGO
  y scoring predictivo fuera de muestra según corresponda.
