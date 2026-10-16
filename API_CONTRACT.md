# Contrato de la API pública

## Propósito

Este documento define la interfaz pública reusable de LRM, HLRM y CHLRM. Es
normativo para el código, las pruebas, la tesis y la evaluación del paper: un
análisis externo debe poder usar estas operaciones sin conocer kernels Gibbs,
actualizaciones variacionales, nombres internos ni la estructura mutable de
los objetos R6.

El repositorio ya construye un paquete R instalable con `DESCRIPTION`,
`NAMESPACE`, documentación inicial y `R CMD check --no-manual` en estado OK.
La superficie pública es el namespace exportado; los kernels y clases R6 son
detalles internos.

## Principios

1. La API se organiza por preguntas estadísticas, no por archivos ni por
   detalles de MCMC/CAVI/SVI.
2. LRM, HLRM y CHLRM usan el mismo nombre cuando una operación tiene el mismo
   significado.
3. Una cantidad que no existe para un modelo se omite; no se rellena con
   `NULL`, cero ni un estimador ficticio.
4. El objetivo predictivo y la unidad de evaluación quedan registrados en el
   resultado.
5. Evaluar no modifica el ajuste: no existe un campo mutable
   `fit_evaluation`.
6. Las matrices potencialmente grandes sólo se conservan cuando el usuario lo
   solicita y deben respetar un presupuesto explícito.
7. Una semilla suministrada hace reproducible la operación y restaura el RNG
   del llamador.
8. Ninguna función pública conoce experimentos, escenarios, tablas, rutas ni
   tamaños fijados por una tesis o artículo.

## Capas de la API

| Capa | Responsabilidad | Pública |
|---|---|---:|
| Construcción | validar datos, prior y configuración; crear el modelo | sí |
| Inferencia | ajustar mediante MCMC, CAVI o SVI | sí |
| Resumen | resumir la distribución posterior o aproximación `q` | sí |
| Predicción | distribución de la media o de una respuesta futura | sí |
| Evaluación | ajuste descriptivo, PPC y criterios predictivos | sí |
| Likelihood pointwise | permitir WAIC, CV, PSIS-LOO y auditoría avanzada | sí, avanzada |
| Orquestación | multichain, multistart, checkpoint y reanudación | sí |
| Kernels | pasos Gibbs, actualizaciones CAVI/SVI, álgebra auxiliar | no |
| Protocolo del paper | escenarios, réplicas, tablas y figuras E1--E5 | no; consume la API |

Que el paper necesite una capacidad no obliga a exportar su receta exacta. La
API exporta el mecanismo general; `paper/evaluation/` fija sus argumentos y
combina sus resultados.

## Modelos y motores disponibles

| Modelo | MCMC actual/final | CAVI actual/final | SVI actual/final |
|---|---:|---:|---:|
| LRM | sí/sí | sí/sí | sí/sí |
| HLRM | sí/sí | sí/sí | sí/sí |
| CHLRM | sí/sí | sí/sí | sí/sí |

Los tres motores están implementados en los tres modelos. SVI para LRM usa
minibatches de observaciones, actualizaciones naturales y reproduce CAVI bajo
batch completo con tasa uno.
`CAVI` es el nombre canónico del algoritmo determinista; `VI` describe la
familia general que contiene tanto CAVI como SVI y no será un valor ambiguo del
campo `method`.

## Contrato canónico de construcción y ajuste

La API final no expondrá `hyperparams` ni `run_settings` como listas abiertas.
Separará objetos con responsabilidades y validadores propios:

| Objeto | Pregunta que responde | No contiene |
|---|---|---|
| `model_spec` | ¿Cuál es el modelo probabilístico y su likelihood? | prior, inicialización o controles del algoritmo |
| `prior_spec` | ¿Qué distribución previa tiene cada parámetro? | estimadores OLS o puntos iniciales |
| `initialization_spec` | ¿Cómo se construye el estado inicial de cada cadena/inicio? | información que altere la prior |
| `inference_control` | ¿Cómo se ejecuta y cuándo se detiene el motor? | parámetros del modelo o del DGP |
| `data_scaling_spec` | ¿Qué columnas/respuesta se centran o escalan? | una copia modificada silenciosamente de los datos |
| `simulation_spec` | ¿Cómo se generan diseño, tamaños y particiones? | decisiones del algoritmo de inferencia |

El término *hiperparámetro* se reservará para parámetros de una distribución
previa o de la jerarquía probabilística. `learning_rate`, `minibatch_size`,
tolerancias, iteraciones y número de cadenas son controles del algoritmo, no
hiperparámetros del modelo.

La interfaz objetivo será:

```r
model <- chlrm_model(component_count = 3)
prior <- chlrm_prior(
  coefficient_mean = rep(0, 3),
  coefficient_covariance = diag(3),
  coefficient_dispersion_df = 7,
  coefficient_dispersion_scale = 3 * diag(3),
  residual_degrees_of_freedom = 6,
  residual_variance_scale_shape = 3,
  residual_variance_scale_rate = 3,
  mixture_concentration = rep(1, 3)
)

fit <- fit_model(
  model = model,
  data = data,
  prior = prior,
  method = "svi",
  initialization = initialization_kmeans(),
  control = svi_control(
    max_iterations = 20000,
    minibatch_size = 20,
    learning_rate = power_schedule(offset = 100, exponent = 0.7),
    stopping = svi_stopping(check_every = 50, patience = 5)
  ),
  scaling = data_scaling(predictors = "continuous", response = TRUE),
  seed = 1
)

posterior <- summarize_posterior(fit)
in_sample <- evaluate_in_sample_fit(fit)
ppc <- posterior_predictive_check(fit, n_draws = 500, seed = 2)
criteria <- predictive_criteria(
  fit,
  n_draws = 500,
  seed = 3,
  predictive_target = "new_groups",
  return_log_likelihood = TRUE
)

mean_prediction <- predict(fit, newdata, target = "mean", seed = 4)
response_prediction <- posterior_predict(fit, newdata, seed = 5)
```

`method` acepta `"mcmc"`, `"cavi"` o `"svi"`. El ajuste conserva una copia
validada de las especificaciones, la transformación aprendida exclusivamente
en entrenamiento, las semillas y la procedencia. El backend R6 podrá seguir
siendo mutable, pero el objeto público será un ajuste S3 estable.

Cada objeto backend representa un solo ajuste. Una segunda llamada interna a
`fit()` se rechaza. MCMC se extiende explícitamente mediante `resume_mcmc()`;
no se vuelve a ejecutar el ajuste desde cero sobre el mismo estado.

### Backends variacionales compartidos

| Modelo | MCMC | CAVI/SVI |
|---|---|---|
| LRM | `mcmc_lrm$new(...)` | `vi_lrm$new(...)`; CAVI y SVI |
| HLRM | `mcmc_hlrm$new(...)` | `vi_hlrm$new(...)` |
| CHLRM | `mcmc_chlrm$new(...)` | `vi_chlrm$new(...)` |

En los tres modelos, `run_settings$stochastic = FALSE` selecciona internamente
CAVI y `TRUE` selecciona SVI. Los contratos nativos y públicos usan `CAVI` o
`SVI`; `VI` sólo nombra la familia estadística. El booleano y las clases R6 no
forman parte de la API pública.

### Inicialización uniforme

`start_point` no será un nombre público: un único punto no describe las cadenas
MCMC, las covarianzas de los factores variacionales ni las responsabilidades
CHLRM. El argumento común será `initialization`, construido mediante:

- `initialization_from_prior()`;
- `initialization_ols()`;
- `initialization_kmeans()` para CHLRM; o
- `initialization_values()` para valores completos suministrados.

Cada estrategia materializa un estado válido para el modelo y motor elegido,
puede generar estados independientes por cadena/inicio y registra semilla y
procedencia. Una estrategia incompatible falla explícitamente; nunca cambia la
prior. CAVI y SVI podrán recibir el mismo conjunto de estados iniciales para
comparaciones pareadas.

### Control y stopping uniformes

El argumento siempre se llama `control`, pero su clase depende del motor:

- `mcmc_control()` contiene warmup, draws, thinning, cadenas, bloques y
  `mcmc_stopping()` basado en R-hat, bulk/tail ESS y MCSE relativo a la
  desviación posterior;
- `cavi_control()` contiene máximo de iteraciones y `cavi_stopping()` basado en
  ELBO determinista, cambio de parámetros, paciencia y mínimo de iteraciones;
- `svi_control()` contiene minibatch, schedule de aprendizaje y
  `svi_stopping()` basado en evaluaciones periódicas comparables del ELBO
  completo, suavizado, cambio de parámetros, paciencia y mínimo de iteraciones.

SVI tiene stopper. No usa como criterio la ELBO de minibatches diferentes,
porque esas cantidades no son directamente comparables. La tasa de aprendizaje
y el stopper son responsabilidades separadas: el schedule decide el tamaño del
paso; el stopper decide si existe evidencia suficiente de estabilización. Todo
motor conserva además un límite duro de trabajo.

Los tres controles devuelven un mismo contrato `convergence-1`:

```r
list(
  status = "converged",       # converged, max_work, failed, interrupted
  reason = "objective_and_parameters_stable",
  work_completed = 4200L,
  checks_completed = 64L,
  metrics = list(...),
  thresholds = list(...),
  history = data.frame(...)
)
```

En la implementación vigente, `mcmc_stopping()` se exige sólo cuando hay dos o
más cadenas. Una cadena alcanza su presupuesto pero no puede producir R-hat y
por ello no se etiqueta como convergida. `starts > 1` activa la orquestación
multistart de CAVI/SVI; cada inicio conserva semillas y diagnóstico, y la
selección usa el mayor ELBO terminal finito entre los inicios que satisfacen la
regla declarada.

`metrics` sólo contiene diagnósticos aplicables. R-hat/ESS no se inventan para
CAVI/SVI y la tolerancia ELBO no se presenta como diagnóstico MCMC.

La equivalencia está en la interfaz y el significado del resultado, no en
forzar la misma regla matemática a algoritmos distintos:

| Motor | Cantidad monitoreada | Regla de elegibilidad para detenerse |
|---|---|---|
| MCMC | R-hat rank-normalized, bulk/tail ESS y MCSE en bloques | todos los parámetros monitoreados cumplen umbrales y el trabajo mínimo |
| CAVI | ELBO completo determinista y cambio de parámetros variacionales | ambos son estables durante `patience` comprobaciones después de `min_iterations` |
| SVI | ELBO completo evaluado periódicamente, tendencia suavizada y cambio de parámetros naturales | objetivo y parámetros son estables durante `patience`; nunca se comparan minibatches distintos |

Los nombres transversales no cambian entre modelos:

| Responsabilidad | Nombre público |
|---|---|
| límite duro | `max_iterations` en CAVI/SVI; `max_post_warmup_draws` en MCMC |
| trabajo mínimo | `min_iterations` o `min_post_warmup_draws` |
| frecuencia diagnóstica | `check_every` |
| estabilidad consecutiva | `patience` |
| objetivo variacional | `objective_tolerance` |
| estado variacional | `parameter_tolerance` |
| minibatch | `minibatch_size`; `sampling_unit` lo determina el modelo y se registra |
| tasa SVI | objeto `learning_rate`, por ejemplo `power_schedule()` |

`iterations`, `epsilon`, `elbo_every`, `minibatch_m`, `tau` y `kappa` podrán
seguir apareciendo transitoriamente dentro de los kernels, pero no serán nombres
públicos alternativos. `offset` y `exponent` son los nombres de usuario para
los parámetros de `power_schedule()`; el objeto registra también su notación
teórica $(\tau,\kappa)$. En SVI la unidad es observación para LRM y grupo para
HLRM/CHLRM, porque los parámetros locales de un grupo deben actualizarse en
conjunto; el usuario controla el tamaño, no redefine esa unidad.

`max_iterations`/`max_post_warmup_draws` es siempre un límite de seguridad, no evidencia de
convergencia. El valor final de `status` distingue `converged` de `max_work`.
Para SVI, `learning_rate` y `stopping` son independientes: una tasa pequeña no
se interpreta por sí sola como convergencia. El ELBO estocástico por minibatch
puede guardarse como diagnóstico del optimizador, pero no activa la parada.

## Operaciones públicas sobre ajustes

### `fit_model()`

```r
fit_model(model, data, prior, method, initialization, control,
          scaling = NULL, seed = NULL, verbose = TRUE)
```

Pregunta: **¿cuál es el resultado de ejecutar el algoritmo de inferencia?**

Devuelve `fit-result-2`, con:

- `parameters`: draws MCMC o parámetros de los factores variacionales;
- `objective`: log-likelihood/log-joint para MCMC o ELBO para CAVI/SVI;
- `diagnostics`;
- `timing`;
- `metadata`;
- `latent_variables` únicamente para CHLRM.

La log-likelihood total almacenada por MCMC y el ELBO almacenado por CAVI/SVI son
objetivos algorítmicos diferentes. No deben compararse como una misma métrica.

### `summarize_posterior()`

```r
summarize_posterior(fit)
```

Pregunta: **¿qué aprendimos sobre los parámetros después de observar los
datos?**

Devuelve `posterior-summary-1`. Cada parámetro aplicable usa `mean`, `sd`,
`lower`, `upper`, `level`, `distribution`, `mean_exists` y `contains_zero`.
CHLRM añade asignaciones, entropía, coclustering y metadatos de relabeling.

### `posterior_draws()`

```r
posterior_draws(
  fit,
  n_draws = NULL,
  seed = NULL,
  include_latent = FALSE,
  max_output_mb = 256
)
```

Pregunta: **¿qué realizaciones conjuntas de la distribución aproximada puedo
usar para tablas, figuras o cálculos posteriores?**

Devuelve `posterior-draws-1`. Los nombres de parámetros son los mismos de
`fit-result-2` y la primera dimensión siempre identifica el draw. Para MCMC,
`draw_index` conserva cadena e iteración; para CAVI/SVI, `distribution =
"variational_approximation"` deja explícito que son muestras de (q), no
cadenas MCMC y por tanto no reciben R-hat ni ESS ficticios. Si
`include_latent = TRUE`, sólo CHLRM añade `group_allocations` en el bloque
separado `latent_variables`. El presupuesto se valida antes de materializar
los arrays.

Con `n_draws = NULL`, MCMC devuelve todos los draws almacenados y CAVI/SVI usa
`evaluation_draws` del control. `posterior::as_draws()` acepta tanto el ajuste
como `posterior-draws-1` y devuelve `draws_df`. Los tensores usan nombres
indexados (`component_coefficients[k,p]`). En MCMC, `.chain` se conserva;
`posterior` reindexa `.iteration` para formar una estructura rectangular, por
lo que la iteración y posición almacenada originales también se preservan en
`source_iteration` y `source_stored_draw`. Esta conversión no modifica el
contrato semántico ni convierte draws de `q` en cadenas MCMC.

### `evaluate_in_sample_fit()`

```r
evaluate_in_sample_fit(fit, return_details = FALSE)
```

Pregunta: **¿qué tan cerca está la media ajustada de los datos usados para
entrenar?**

Devuelve `in-sample-fit-1`:

- `statistics$mse`;
- `statistics$r2`;
- definición del estimador en `metadata`;
- `fitted_values` y `residuals` sólo si `return_details = TRUE`.

Es ajuste descriptivo clásico dentro de muestra. No mide calibración posterior
ni generalización.

### `posterior_predictive_check()`

```r
posterior_predictive_check(
  fit,
  n_draws = NULL,
  seed = NULL,
  return_yrep = FALSE,
  max_output_mb = 256,
  discrepancies = ppc_discrepancies()
)
```

Pregunta: **si el modelo ajustado generara nuevos datos, produciría
características parecidas a las observadas?**

Devuelve `posterior-predictive-check-2`:

- discrepancias replicadas;
- discrepancias observadas;
- probabilidades predictivas posteriores;
- `mcse_iid` nominal y `mcse`/`effective_sample_size`: para MCMC se calculan
  desde la secuencia indicadora ordenada por cadena; para CAVI/SVI corresponden
  a draws independientes de la aproximación `q`;
- `y_rep` sólo si `return_yrep = TRUE`.

El preset estándar calcula media global (T1), desviación estándar global (T2)
y MSE frente a un OLS observado fijo (T5). En datos agrupados añade varianza de
medias grupales (T3) y promedio de varianzas intragrupo (T4) cuando sus tamaños
las hacen definibles. `ppc_discrepancies()` también permite seleccionar un
subconjunto; una discrepancia no aplicable produce un error explícito, no un
campo ficticio.

Es una comprobación de adecuación generativa, no una prueba de hipótesis
calibrada ni una evaluación fuera de muestra.

### `predictive_criteria()`

Pregunta: **¿qué capacidad predictiva estimada tiene el modelo para el objetivo
y la unidad declarados?**

LRM:

```r
predictive_criteria(
  fit,
  n_draws = NULL,
  seed = NULL,
  predictive_target = "observed_data",
  predictive_unit = c("observation", "group"),
  return_log_likelihood = FALSE,
  max_output_mb = 256,
  groups = NULL
)
```

HLRM:

```r
predictive_criteria(
  fit,
  n_draws = NULL,
  seed = NULL,
  predictive_target = c("observed_groups", "new_groups"),
  return_log_likelihood = FALSE,
  max_output_mb = 256,
  marginal_control = list()
)
```

CHLRM:

```r
predictive_criteria(
  fit,
  n_draws = NULL,
  seed = NULL,
  predictive_target = "new_groups",
  return_log_likelihood = FALSE,
  max_output_mb = 256
)
```

Devuelve `predictive-criteria-2`, con WAIC y sus diagnósticos aplicables dentro
de `criteria`. DIC se retira de la API final: depende de un punto posterior y de
la parametrización, puede ser inestable en jerarquías/mezclas y no añade una
base defendible para seleccionar componentes.

Con `return_log_likelihood = TRUE` añade una lista nombrada de objetos
`log-likelihood-1`. Cada matriz tiene unidades predictivas en filas y draws en
columnas:

$$
L_{us}=\log p(y_u\mid\theta^{(s)}).
$$

Esta matriz permite auditar WAIC y alimentar `psis_loo()` u otros scores sin
volver a ajustar el modelo. CV sí reajusta por definición en cada fold.

### `predict()`

```r
predict(
  fit,
  newdata,
  group_ids = NULL,
  n_draws = NULL,
  level = 0.95,
  seed = NULL,
  return_draws = FALSE,
  max_output_mb = 256
)
```

Pregunta: **¿cuál es la distribución posterior de la media de respuesta?**

Devuelve `posterior-prediction-2` con `target = "mean"`, resumen e intervalo
creíble puntual. LRM no usa `group_ids`; HLRM y CHLRM distinguen grupos
observados de nuevos.

### `posterior_predict()`

Usa los mismos argumentos de `predict()`, pero por defecto conserva draws:

```r
posterior_predict(
  fit,
  newdata,
  group_ids = NULL,
  n_draws = NULL,
  level = 0.95,
  seed = NULL,
  return_draws = TRUE,
  max_output_mb = 256
)
```

Pregunta: **¿cuál es la distribución de una respuesta futura incluyendo ruido
observacional?**

Devuelve `posterior-prediction-2` con `target = "response"` e intervalo
predictivo puntual. Con `return_draws = TRUE` conserva `y_rep` agrupado.

## Propiedad y selección automática de la likelihood

La likelihood pertenece a `model_spec`, no a MCMC, CAVI, SVI ni a GOF. Cada
familia tiene un único núcleo probabilístico validado que alimenta:

- el log joint o ELBO del ajuste;
- `log_likelihood()` y los criterios predictivos;
- `posterior_predict()`; y
- la parte condicional de `simulate_data()`.

Las representaciones computacionales pueden diferir —estadísticas suficientes,
vectorización, cuadratura o `log-sum-exp`—, pero deben aprobar pruebas de
equivalencia contra la misma densidad. Cambiar de motor de inferencia nunca
cambia el modelo observacional.

El usuario de alto nivel declara qué quiere predecir. El backend decide qué
debe condicionarse o integrarse:

| Modelo | Objetivo | Unidad | Likelihood |
|---|---|---|---|
| LRM | `observed_data` | observación | condicional en coeficientes y varianza |
| LRM | `observed_data` | grupo suministrado | misma likelihood agregada; no es marginal |
| HLRM | `observed_groups` | observación y grupo | condicional en parámetros locales |
| HLRM | `new_groups` | grupo | marginaliza coeficientes y varianzas locales |
| CHLRM | `new_groups` | grupo | marginaliza la asignación de componente |

La API avanzada `log_likelihood()` permite declarar directamente `unit` y
`conditioning`. Es pública porque paquetes de modelado bayesiano normalmente
exponen la log-likelihood pointwise para herramientas de diagnóstico y
comparación. Los usuarios ordinarios deberían preferir
`predictive_criteria()`, que evita combinaciones incompatibles.

## Contratos de resultados

La tabla describe el contrato objetivo. `fit_model()` ya produce
`fit-result-2`; los backends directos todavía producen `fit-result-1` y no
constituyen la superficie que se congelará. El simulador ya produce únicamente
`simulation-result-2`; sus consumidores fueron migrados sin conservar el schema
anterior como alias.

| Esquema | Clase | Productor |
|---|---|---|
| `fit-result-2` | `vb_fit_result` | `fit_model()` |
| `posterior-summary-1` | `vb_posterior_summary` | `summarize_posterior()` |
| `posterior-draws-1` | `vb_posterior_draws` | `posterior_draws()` |
| `in-sample-fit-1` | `vb_in_sample_fit` | `evaluate_in_sample_fit()` |
| `posterior-predictive-check-2` | `vb_posterior_predictive_check` | `posterior_predictive_check()` |
| `predictive-criteria-2` | `vb_predictive_criteria` | `predictive_criteria()` |
| `posterior-prediction-2` | `vb_prediction_result` | `predict()`, `posterior_predict()` |
| `log-likelihood-1` | `vb_log_likelihood` | `log_likelihood()` o `predictive_criteria()` |
| `simulation-result-2` | `vb_simulation_result` | `simulate_data()` |
| `fit-resource-estimate-1` | `vb_fit_resource_estimate` | `estimate_fit_resources()` |

Cambiar la estructura o el significado de un esquema requiere incrementar su
versión y migrar pruebas y consumidores. No se cambia únicamente la etiqueta
de resultados antiguos.

## Funciones públicas adicionales

### Construcción, datos, simulación y priors

- `lrm_model()`, `hlrm_model()` y `chlrm_model()` describen el modelo y su
  likelihood, pero no la prior ni el algoritmo;
- `lrm_prior()`, `hlrm_prior()` y `chlrm_prior()` construyen priors explícitas;
- `prior_moments()` permite auditar escala y existencia de momentos antes de
  ajustar; la simulación prior-predictiva usa
  `prior_predictive_dgp()` + `simulate_data()` sin un wrapper redundante;
- `initialization_from_prior()`,
  `initialization_ols()`, `initialization_kmeans()` e
  `initialization_values()` construyen inicializaciones sin alterar la prior;
- `mcmc_control()`, `cavi_control()` y `svi_control()` contienen únicamente
  decisiones del algoritmo; los tres aceptan `max_memory_mb`, incluso con una
  sola cadena o inicio;
- `estimate_fit_resources()` valida modelo/datos/control y separa pico
  conservador de ajuste, estado retenido y almacenamiento potencial de una
  evaluación posterior. Es una guarda previa, no una promesa de RSS del
  sistema operativo;
- `grouped_regression_data()` crea desde un `data.frame` la representación
  agrupada que comparten HLRM y CHLRM, con mapeo y procedencia de columnas;
- `data_scaling()`, `fit_scaling()`, `transform_data()`, `transform_prior()` y
  `scaling_map()` implementan y permiten auditar el escalado explícito; la
  restauración de resultados ocurre uniformemente en la fachada del ajuste; y
- `conditional_dgp()`, `hierarchical_dgp()` y `prior_predictive_dgp()` describen
  mecanismos generadores con condicionamiento explícito;
- `simulation_design()`, `simulation_sample()` y `simulate_data()` forman la
  interfaz común de simulación.

No se exportarán listas abiertas llamadas `hyperparams` ni constructores de
“valores iniciales por modelo”. La clase de `model` permite que una misma
estrategia de inicialización materialice el estado que corresponde a
LRM/HLRM/CHLRM y MCMC/CAVI/SVI. Una prior calibrada mediante OLS, si se conserva
como conveniencia, tendrá un constructor separado y explícitamente marcado como
empírico-Bayes; nunca será el default ni una opción denominada “no informativa”.

Los simuladores históricos y los helpers nombrados por el protocolo ya fueron
retirados. El mecanismo general usa los tres DGP tipados y los valores exactos
del PDF se construyen en `paper/evaluation/R/scenarios.R`. `hyperparams_*()`
permanece sólo en los backends/consumidores aún no migrados y no formará parte
del namespace público final.

#### Contrato de simulación

`truth` no es una entrada ni una estrategia: es el conjunto completo de
parámetros y variables latentes realizados que devuelve cada simulación. La
entrada es un mecanismo generador de datos (DGP) que declara qué nivel se
condiciona:

| DGP | Se mantiene fijo entre réplicas | Se vuelve a generar | Pregunta principal |
|---|---|---|---|
| `conditional_dgp()` | parámetros que entran directamente en la likelihood | latentes observacionales no suministradas y respuesta | recuperación/estrés en una configuración concreta |
| `hierarchical_dgp()` | parámetros poblacionales de la jerarquía | parámetros locales/de componente, latentes y respuesta | desempeño integrado para nuevas unidades de la población |
| `prior_predictive_dgp()` | diseño y prior | todos los parámetros, latentes y respuesta desde $p(\theta,y)$ | prior predictiva y calibración computacional SBC |

`hierarchical_dgp()` sólo aplica a HLRM y CHLRM. En LRM no existe un nivel local
intermedio: se usa `conditional_dgp()` para recuperación con
$(\beta,\sigma^2)$ fijos o `prior_predictive_dgp()` para muestrearlos desde la
prior. El modo condicional sigue siendo probabilístico: condiciona parámetros
conocidos y genera las variables aleatorias restantes.

En CHLRM los parámetros directos son coeficientes/varianzas de componente y
pesos. Las asignaciones de grupo se muestrean desde los pesos por defecto, pues
son latentes por dataset, pero pueden suministrarse para diseños pareados o
train/test que condicionen en las mismas asignaciones. Pesos poblacionales y
proporciones empíricas realizadas siempre se conservan por separado.

Todos los generadores:

- producen el diseño `X` y la respuesta `y` desde una distribución de
  covariables explícita;
- aíslan/restauran RNG y validan dimensiones, soporte y matrices SPD;
- devuelven `simulation-result-2` con `train`, `truth`, `configuration`,
  `conditioning`, `parameter_source` y `design_source`, y
  `test` sólo cuando fue generado;
- no estiman parámetros, no ejecutan réplicas y no conocen escenarios del
  paper.

La interfaz uniforme es:

```r
sim <- simulate_data(
  generator = conditional_dgp(
    model = chlrm_model(component_count = 3),
    parameters = list(
      component_coefficients = beta,
      component_variances = sigma2,
      mixture_weights = weights
    )
  ),
  design = simulation_design(
    distribution = "multivariate_normal",
    covariance = Sigma_x,
    include_intercept = TRUE
  ),
  sample = simulation_sample(group_count = 20, observations_per_group = 30),
  seed = 1
)
```

Por defecto `truth` conserva parámetros y latentes necesarios para evaluar
recuperación, pero no duplica medias condicionales ni innovaciones para cada
observación. Éstas pueden solicitarse con
`retain_observation_truth = TRUE`. `max_output_mb` fija un presupuesto local y
rechaza la solicitud antes de materializar arreglos excesivos.

Los tres DGP usan el mismo `simulate_data()` y devuelven la misma envoltura, sin
campos ficticios. La geometría que construye $\beta_k$ a distancia mínima `d`,
los valores S1--S7c y sus repeticiones no forman parte de la API: la suite
construye esos coeficientes y los entrega a `conditional_dgp()`. La generación
predictiva posterior corresponde a `posterior_predict()` sobre un ajuste y no
es un cuarto modo del generador Monte Carlo.

La correspondencia por modelo es:

| Modelo | `conditional_dgp()` | `hierarchical_dgp()` | `prior_predictive_dgp()` |
|---|---|---|---|
| LRM | fija $\beta,\sigma^2$; genera $y$ | no aplica | genera $\beta,\sigma^2,y$ desde prior y likelihood |
| HLRM | fija $\beta_j,\sigma_j^2$; genera $y_j$ | fija $\beta,\Sigma,\xi^2$; genera efectos/varianzas de grupo y $y_j$ | genera jerarquía completa y datos |
| CHLRM | fija $\beta_k,\sigma_k^2,\omega$; genera $z_j,y_j$ | fija $\beta,\Sigma,\xi^2,\omega$; genera componentes, asignaciones y datos | genera jerarquía completa, pesos, asignaciones y datos |

En cada modo `simulation_design()` declara si $X$ es fijo o aleatorio. Esa
decisión no se infiere del modelo y queda en metadatos, porque condicionarse en
un diseño observado y promediar sobre una población de covariables responden
preguntas distintas.

Los modos no son sustitutos y sus métricas no se interpretan igual:

- repetir $y$ con los mismos parámetros estima desempeño frecuentista
  **condicional** en ese escenario; no es calibración bayesiana marginal;
- volver a muestrear parámetros locales desde una población permite evaluar
  recuperación de parámetros poblacionales y predicción de grupos nuevos;
- SBC exige que el DGP conjunto use exactamente la misma prior y likelihood que
  el ajuste. Si se altera una de ellas, el resultado pasa a ser un estudio de
  robustez/mala especificación y no una prueba SBC.

Para HLRM, un DGP condicional con $\beta_j,\sigma_j^2$ fijos no basta para
afirmar recuperación repetida de $\beta,\Sigma,\xi^2$. Para CHLRM, recuperación
de componentes se evalúa después de relabeling o mediante cantidades
invariantes a etiquetas; SBC de etiquetas crudas no es interpretable. El
resultado registra `conditioning`, `parameter_source`, `design_source` y los
estimandos que el DGP permite evaluar.

#### Contrato de priors e hiperparámetros

Los parámetros generadores y los hiperparámetros de inferencia deben compartir
familia, soporte y parametrización, pero no valores automáticamente:

| Modelo | Verdad/DGP | Prior de inferencia |
|---|---|---|
| LRM | `coefficients`, `residual_variance` | media/covarianza de coeficientes y parámetros inverse-gamma de la varianza |
| HLRM | `population_coefficients`, `population_covariance`, `residual_variance_scale`, `residual_degrees_of_freedom` | `Mu0`, `Lambda0`, `n0`, `S0`, `nu0`, `a0`, `b0` internamente; nombres descriptivos en el constructor público |
| CHLRM | lo anterior más `mixture_weights`, parámetros de componente y asignaciones | lo anterior más concentración Dirichlet para pesos |

El simulador devuelve la verdad para evaluación, pero el ajuste recibe una
prior construida por separado. Copiar la verdad a la prior sería fuga de
información salvo en una prueba algebraica explícitamente marcada. MCMC, CAVI
y SVI de un mismo modelo reciben exactamente el mismo objeto prior validado.
En SBC, el DGP y el ajuste sí comparten el mismo **objeto prior**, pero el ajuste
nunca recibe el valor $\theta$ realizado que generó ese dataset.
`shape`/`rate`, covariance/precision y grados de libertad se nombran sin
ambigüedad; inversas derivadas no son entradas públicas independientes.

Los nombres públicos serán descriptivos y el traductor interno será único. Un
nombre compartido representa la misma distribución y parametrización:

| Bloque de prior | LRM | HLRM | CHLRM |
|---|---:|---:|---:|
| `coefficient_mean`, `coefficient_covariance` para el coeficiente global | sí | sí | sí |
| `coefficient_dispersion_df`, `coefficient_dispersion_scale` para la inverse-Wishart de la covarianza local | — | sí | sí |
| `residual_degrees_of_freedom` | sí | sí | sí |
| `residual_variance_scale` de la inverse-gamma residual | sí | — | — |
| `residual_variance_scale_shape`, `residual_variance_scale_rate` para $\xi^2$ | — | sí | sí |
| `mixture_concentration` para los pesos Dirichlet | — | — | sí |

En forma matemática:

$$
\begin{aligned}
\boldsymbol\beta &\sim N(\texttt{coefficient\_mean},
  \texttt{coefficient\_covariance}),\\
\boldsymbol\Sigma &\sim IW(\texttt{coefficient\_dispersion\_df},
  \texttt{coefficient\_dispersion\_scale}),\\
\sigma_j^2\mid\xi^2 &\sim IG(\nu/2,\nu\xi^2/2),\qquad
\xi^2\sim Gamma(a,b),\\
\boldsymbol\omega &\sim Dirichlet(\texttt{mixture\_concentration}).
\end{aligned}
$$

La inverse-Wishart sólo pertenece a HLRM/CHLRM y la Dirichlet sólo a CHLRM. La
inverse-gamma residual pertenece a los tres; en LRM se reemplaza $\xi^2$ por
`residual_variance_scale` fija. Los guiones de la tabla significan “no
pertenece al modelo”, no un campo enviado como `NULL`.

`Beta0`, `Sigma0`, `Mu0`, `Lambda0`, `n0`, `S0`, `nu0`, `a0`, `b0` y
`alpha0` pueden conservarse dentro de los kernels porque coinciden con las
ecuaciones, pero no serán la única documentación de la escala. Cada prior
registrará familia y parametrización (`shape--scale` para inverse-gamma,
`shape--rate` para gamma, escala —no precisión— para inverse-Wishart).

No existirá una opción denominada `NON_INFORMATIVE`: el código actual calcula
su escala con residuos OLS y por tanto es empírico-Bayes. Se separarán:

1. `lrm_prior()`/`hlrm_prior()`/`chlrm_prior()` para una prior explícita e
   independiente de los datos;
2. `initialization_ols()` para inicialización, sin modificar la prior.

La primera API estable no expondrá `empirical_prior()`: evita convertir una
receta OLS dependiente de datos en una aparente recomendación general. Si una
aplicación necesita empírico-Bayes, deberá construir y registrar explícitamente
la prior resultante fuera del ajuste; esa política puede revisarse después de
validación prior-predictiva específica del dominio.

Los objetos de prior reportarán la existencia de sus momentos. Por ejemplo,
con $\sigma^2\mid\xi^2\sim IG(\nu_0/2,\nu_0\xi^2/2)$, la media condicional no
existe cuando $\nu_0\leq2$; en ese caso $\xi^2$ no se describirá erróneamente
como la media de $\sigma^2$. `prior_moments()` y simulación prior-predictiva
permitirán auditar la escala antes de ajustar.

No se promete una “mejor prior” universal en una escala arbitraria. Para la
primera API estable se exige una prior explícita. Un futuro
`regularizing_prior()` sólo podrá ofrecerse en una escala declarada, con
familias conjugadas compartidas por los tres motores, momentos documentados y
evidencia prior-predictiva. La v1 no expondrá un constructor empírico-Bayes.

#### Contrato de escalado de datos

Se elimina el nombre `vb_standardizer`. El constructor público será
`data_scaling()` y producirá un `data_scaling_spec`; `fit_scaling()` producirá
un `data_scaling_map` aprendido sólo con entrenamiento. Los nombres describen
la responsabilidad sin ligar el preprocesamiento a VI ni presentarlo como un
modelo.

La selección de columnas continuas es explícita:

```r
scaling <- data_scaling(
  predictors = c("diameter", "height"),
  response = TRUE,
  intercept_column = "intercept"
)
map <- fit_scaling(scaling, training_data)
scaled_training <- transform_data(training_data, map)
scaled_prior <- transform_prior(prior, map)
```

`predictors = "all"` selecciona todas las columnas excepto el intercepto.
Columnas no seleccionadas —por ejemplo indicadores categóricos— permanecen
sin modificar. El paquete no intenta adivinar si una columna numérica es
continua o categórica. Nombres, orden e intercepto se validan también al
transformar datos nuevos.

Cuando `fit_model(..., scaling = data_scaling(...))` recibe la especificación,
el paso es explícito aunque su aplicación sea automática: el ajuste aprende el
mapa sobre entrenamiento, transforma datos, prior e inicialización, y conserva
el mapa. Las priors públicas se declaran siempre en la escala de entrada. Los
parámetros, draws, resúmenes, predicciones, PPC y criterios se devuelven
siempre en esa escala. `scaling_map(fit)` permite auditar el mapa; la escala
numérica del backend no se ofrece como una segunda modalidad pública. El
usuario avanzado puede ejecutar por separado `fit_scaling()`,
`transform_data()` y `transform_prior()`.

No existe un booleano ambiguo `standardize = TRUE` ni transformación implícita
por motor. MCMC, CAVI y SVI reciben exactamente los mismos datos y prior ya
transformados. En CV, cada fold aprende su propio mapa sólo con su conjunto de
entrenamiento. Los simuladores siempre devuelven la escala del DGP.

La transformación cubre toda la jerarquía: medias de coeficientes mediante el
mapa afín; covarianzas `Sigma0`, `Lambda0`, `S0` y covarianzas posteriores
mediante el mapa lineal; varianzas residuales y $\xi^2$ mediante $s_y^2$; y el
rate gamma de $\xi^2$ mediante la transformación de escala correspondiente.
Pesos Dirichlet, responsabilidades, asignaciones y $K$ no cambian. Aplicar el
mismo objeto y la prior transformada hace que los tres motores ajusten el mismo
modelo en la misma escala numérica.

La recuperación sí depende del **tipo de bloque paramétrico**, y por eso se
despacha por modelo, pero no usa tres álgebra distintas:

| Modelo | Bloques restaurados a escala de entrada | Bloques invariantes |
|---|---|---|
| LRM | coeficientes y varianza residual | ninguno |
| HLRM | coeficientes poblacionales/locales, covarianza entre coeficientes, varianzas residuales y $\xi^2$ | identificadores de grupo |
| CHLRM | coeficientes poblacionales/de componente, covarianza entre coeficientes, varianzas residuales y $\xi^2$ | pesos, responsabilidades, asignaciones y $K$ |

Un motor no implementa su propia inversión. Métodos internos por modelo
enumeran los bloques aplicables y llaman un único motor afín. Esto evita tanto
campos `NULL` como fórmulas duplicadas entre MCMC/CAVI/SVI.

Si $\boldsymbol\beta=A+T\boldsymbol\beta^*$ transforma coeficientes de escala
estandarizada a original, entonces
$m^*=T^{-1}(m-A)$ y $V^*=T^{-1}VT^{-\mathsf T}$ transforman una prior original
a la escala de ajuste. Para la jerarquía,
$\Sigma^*=T^{-1}\Sigma T^{-\mathsf T}$. Si
$\xi^2\sim\mathrm{Gamma}(a,b)$ usa `rate = b` y
$\xi^{2*}=\xi^2/s_y^2$, entonces $b^*=b s_y^2$. Estas identidades se probarán
en ambos sentidos y no se rederivarán dentro de cada motor. Para
$\sigma^2\sim IG(a,b)$ con parametrización shape--scale,
$\sigma^{2*}=\sigma^2/s_y^2$ implica $b^*=b/s_y^2$.

La log-likelihood pública siempre se reporta en la escala de entrada. Si se
evalúa en escala transformada, se aplica $-\log s_y$ por observación —o
$-n_j\log s_y$ por grupo— al regresar a la densidad original.

### Ejecución robusta

- `fit_model()` orquesta internamente cadenas MCMC o inicios CAVI/SVI sin
  exponer listas `hyperparams`/`run_settings`;
- `fit_orchestration()` devuelve `fit-orchestration-1`, una auditoría de sólo
  lectura: metadatos/diagnósticos por cadena o resumen por inicio; CHLRM
  multistart añade clasificaciones por inicio para medir estabilidad. No expone
  backends R6, parámetros nativos ni estados parciales;
- `resume_mcmc()` continúa un `vb_model_fit` sin mutarlo y recalcula
  relabeling/diagnósticos multichain;
- `save_mcmc_checkpoint()` y `read_mcmc_checkpoint()` persisten/restauran el
  contrato tipado mediante estado plano, sin serializar entornos R6.

### Evaluación avanzada

- `log_likelihood()`;
- `predictive_log_likelihood()` y `predictive_score()` para respuestas
  retenidas, con ELPD/MLPD y objetivo observado/nuevo explícito;
- `predictive_performance()` para RMSE, MAE, R² de test, cobertura/ancho de
  intervalos predictivos y scores logarítmicos en un contrato separado del
  ajuste dentro de muestra y del PPC;
- `cross_validation_plan()`, `cross_validation_split()` y
  `cross_validate_model()` para CV por observación, dentro de grupo, por grupos
  o LOGO, con semillas y presupuesto de ajustes explícitos; el resultado
  conserva log scores pointwise y, cuando `compute_performance = TRUE`, también
  predicciones retenidas, RMSE, MAE, cobertura y ancho de intervalos;
- `psis_loo()` para PSIS-LOO sobre MCMC convergente, con eficiencia relativa,
  Pareto-k y unidad/condicionamiento preservados;
- `ppc_discrepancies()` para declarar estadísticas PPC aplicables;
- `pp_check_from_prediction()`;
- `waic_from_loglik()` como primitiva transparente mientras
  `predictive_criteria()` es la interfaz recomendada.

## Qué no pertenece a la API pública

- `sample_*`, `update_q_*` e `initialize_*`;
- `loglik_total_*`, `logprior_*` y funciones especializadas de integración;
- validadores y operaciones Cholesky internas;
- preparadores de draws y constructores `new_*_result`;
- helpers de relabeling usados por las clases;
- constructores de escenarios, presets de priors y geometrías exclusivas del
  protocolo;
- rutas, números de réplicas, semillas maestras, tablas y figuras particulares
  de la tesis o del paper.

Estas piezas pueden cambiar para optimizar rendimiento sin romper análisis de
usuarios, siempre que los contratos públicos y su significado permanezcan.

## Comparabilidad

WAIC, PSIS-LOO y log scores sólo se comparan cuando coinciden:

- dataset;
- unidad predictiva;
- objetivo predictivo;
- condicionamiento o marginalización;
- convención de likelihood.

No son comparables directamente WAIC por observación y por grupo, ni HLRM
condicional y marginal. WAIC calculado con draws de `q` evalúa la aproximación
variacional ajustada; no demuestra que CAVI/SVI reproduzca la posterior MCMC.

## Relación con el flujo bayesiano estándar

La separación sigue las responsabilidades habituales de un paquete bayesiano:

1. ajustar y diagnosticar el algoritmo;
2. resumir la posterior;
3. generar cantidades predictivas;
4. comprobar adecuación mediante PPC;
5. exponer log-likelihood pointwise;
6. evaluar predicción y comparar modelos bajo un estimando declarado.

`evaluate_in_sample_fit()` añade el puente con estadística clásica y aprendizaje
automático mediante MSE, R² y residuos. Es útil, pero se mantiene separado del
PPC bayesiano y de los criterios predictivos para no presentar ajuste dentro de
muestra como generalización.

### Correspondencia con el ecosistema R

La interfaz canónica actual es la fachada S3 `vb_model_fit`; los motores R6 son
implementación privada. P4 no copiará las estructuras internas de `lme4`,
`brms`, `mclust` o `mixtools`, pero sí adoptará sus convenciones útiles.
La fachada no vuelve a calcular inferencia: conserva datos, fórmula/diseño,
metadatos y referencia al resultado necesario para despachar métodos estándar.

| Convención externa | Interfaz objetivo | Semántica en este paquete |
|---|---|---|
| `lme4`/base R | `print()`, `summary()`, `coef()`, `fitted()`, `residuals()` y `nobs()` implementados; accesores jerárquicos por decidir | sólo donde la cantidad exista; no se fuerza una matriz de modelo única sobre entradas agrupadas |
| `brms`/`posterior` | `as_draws_*()`, `posterior_epred()`, `posterior_predict()` | draws MCMC conservan cadena; draws CAVI/SVI se identifican explícitamente como muestras de `q` y no reciben R-hat/ESS ficticios |
| `bayesplot` | extractor `as_yrep()` | matriz externa draws × observaciones, sin obligar a guardar `y_rep` en cada ajuste |
| `loo` | extractor `as_log_lik_matrix()` y, después, método `loo()` | matriz externa draws × unidad con objetivo, unidad y condicionamiento preservados |
| `mclust`/`mixtools` | `components()`, `responsibilities()`, `classify_groups()`, `classification_uncertainty()` y `co_clustering()` implementados | sólo CHLRM; no se fuerzan variables latentes sobre LRM/HLRM |

Los contratos actuales de predicción y log-likelihood usan internamente
observación/unidad × draw. `bayesplot`, `brms` y `loo` usan normalmente draw ×
observación/unidad. La compatibilidad se resuelve mediante `as_yrep()` y
`as_log_lik_matrix()`, probados como draws × unidad; no con transposiciones
implícitas sin metadatos.

Se mantendrá la API matricial de bajo nivel. Una interfaz opcional
`formula + data` deberá almacenar `terms`, contrastes y niveles para validar
`newdata`; no cambiará la parametrización matemática del modelo.

No se debe usar `stats::logLik()` para esconder la matriz pointwise: ese
genérico representa una cantidad escalar y perdería la unidad predictiva y el
condicionamiento. La función explícita `log_likelihood()` conserva esa
información.

`posterior` es la dependencia canónica para draws y `loo` implementa
PSIS-LOO. `bayesplot` podrá ser una dependencia sugerida de interoperabilidad.
`lme4`, `brms`,
`mclust` y `mixtools` son referencias de diseño, no dependencias necesarias.

### Puente predictivo

`predictive_score()` ya evalúa respuestas fuera de muestra mediante ELPD y
MLPD, conservando valores pointwise. `cross_validate_model()` usa esa operación
o `predictive_performance()` en folds tipados y aprende cualquier
`data_scaling()` exclusivamente dentro de cada entrenamiento. El paquete
distingue cuatro esquemas:

| Esquema | Pregunta |
|---|---|
| `observation_kfold` | generalización de una regresión LRM a observaciones nuevas |
| `within_group_kfold` | predicción de observaciones nuevas en grupos HLRM/CHLRM ya conocidos |
| `group_kfold` | predicción marginal de grupos completos nuevos |
| `leave_one_group_out` | LOGO marginal, un grupo nuevo por refit |

El runner es secuencial y exige `max_fits`; la paralelización, checkpoints y
grillas masivas pertenecen a la suite del paper. `psis_loo()` delega el
suavizado y Pareto-k al paquete estándar `loo`: no acepta CAVI/SVI porque la
identidad de importancia leave-one-out requiere draws de la posterior, no de
una aproximación `q`. Para esos métodos se usa CV con refit. Las métricas
comunes de test están en `predictive_performance()` y el PPC configurable
reporta MCSE corregido por ESS en MCMC. La agregación entre réplicas pertenece
a la evaluación del paper.

## Serialización y frontera de I/O

El paquete devuelve objetos serializables y sólo escribe archivos cuando el
usuario llama una operación explícita como `save_mcmc_checkpoint()`. No decide
rutas, nombres, compresión, política de retención ni qué subconjunto debe
persistirse.

La suite del paper es responsable del manifiesto, la procedencia, los
checkpoints por tarea y la selección de objetos que necesita para recalcular
sus métricas. Esas obligaciones están en
[`paper/EVALUATION_PLAN.md`](paper/EVALUATION_PLAN.md), no forman parte del
contrato público del paquete.
