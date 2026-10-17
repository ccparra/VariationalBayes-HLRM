# Roadmap del paquete R

## Objetivo

Convertir `R/` en un paquete instalable, probado y documentado que permita
simular, ajustar, diagnosticar, comparar y predecir con LRM, HLRM y CHLRM. El
paquete no conoce E1--E5, tablas del paper, números de réplicas ni rutas de la
tesis.

La frontera se decide por generalidad, no por quién solicitó originalmente una
función. Un simulador parametrizable, una likelihood pointwise o un PPC
configurable pertenecen al paquete. Una receta con S1--S7c, hiperparámetros
fijos, 200 réplicas o nombres de tablas pertenece a `paper/evaluation/`.

La lista normativa de operaciones que deberán exportarse está en
[`API_CONTRACT.md`](API_CONTRACT.md). P4 debe convertir esa lista —y no todos
los símbolos actualmente cargados mediante `source()`— en el `NAMESPACE`.

## Estado actual

Ya están conectados y probados:

- kernels MCMC, CAVI y SVI existentes, con identificadores de método canónicos;
- simuladores LRM/HLRM/CHLRM con contrato versionado y DGP condicional,
  jerárquico y prior-predictivo separados del protocolo del paper;
- likelihood marginal CHLRM y HLRM;
- contratos públicos versionados para ajuste, resumen, draws posteriores,
  ajuste descriptivo, PPC, criterios predictivos y predicción;
- likelihood con unidad/condicionamiento explícitos;
- `predict()` y `posterior_predict()` uniformes;
- PPC streaming y reutilización de `y_rep`;
- RNG/estado/parada robustos en LRM, HLRM y CHLRM;
- MCMC multichain reanudable y CAVI/SVI multistart comparable;
- tests algebraicos, dimensiones unitarias y presupuestos locales.

LRM, HLRM y CHLRM implementan MCMC/CAVI/SVI. En LRM, SVI muestrea
observaciones; en HLRM/CHLRM muestrea grupos completos.

La siguiente tabla distingue código terminado de decisiones que sólo están
especificadas. Una tarea no se considera implementada por aparecer en
`API_CONTRACT.md`.

| Bloque | Estado verificable | Brecha principal |
|---|---|---|
| P1, inferencia para motores existentes | Cerrado | Los nueve motores MCMC/CAVI/SVI están conectados; su profiling final pertenece a P5 |
| P2, resultados y estado de ejecución | Cerrado | fachada tipada, convergencia uniforme y consumidores activos migrados |
| P2.1, interoperabilidad R | Parcial avanzado | fachada S3, métodos base, draws `posterior`, extractores orientados, accesores jerárquicos/CHLRM y salidas tidy existen; falta decidir fórmula/datos |
| P2.2, especificaciones estadísticas | Cerrado | DGP tipados, prior-predictivo, consumidores, procedencia, primitiva observacional y smoke SBC están implementados |
| P2.3, escalado reversible | Cerrado | `fit_model()` comparte el mapa entre motores y restaura parámetros, draws, resúmenes, predicción, PPC, ajuste y densidades; CV aprende el mapa sólo en cada fold de entrenamiento |
| P3, evaluación predictiva | Cerrado | WAIC, scoring retenido, ELPD/MLPD, métricas e intervalos de test, LOGO/CV, PSIS-LOO MCMC y PPC configurable con MCSE dependiente están implementados |
| P4, paquete instalable | Parcial avanzado | paquete plano, imports/exports, documentación inicial y `R CMD check` OK; faltan vignettes y convertir la documentación manual a roxygen |
| P5, rendimiento | Parcial avanzado | existe estimador público y presupuesto uniforme por ajuste; falta profiling empírico en escenarios piloto |

La frontera de simulación ya fue separada: el paquete conserva los DGP y
constructores genéricos; la geometría de coeficientes y la prior exacta del PDF
viven en `paper/evaluation/R/scenarios.R`. Los nombres históricos fueron
retirados después de migrar tests, tesis y experimentos, sin aliases.

## P1. Inferencia comparable

Entregables del paquete:

1. **Implementado:** relabeling CHLRM iterativo con referencia, escalas y
   criterio de convergencia explícitos;
2. **Implementado:** objeto de referencia de etiquetas reutilizable para
   alinear MCMC, CAVI y SVI;
3. **Implementado:** ejecución MCMC multichain general, sin imponer cuatro
   cadenas, con semillas independientes, ejecución secuencial segura y
   combinación bloqueada hasta aprobar diagnósticos;
4. **Implementado:** diagnósticos rank-normalized R-hat, bulk ESS, tail ESS y
   ESS/segundo mediante `posterior`;
5. **Implementado:** ejecución CAVI/SVI multistart general con estados iniciales
   serializables y reutilizables entre métodos, resultados de todos los intentos
   y selección explícita por ELBO entre ajustes elegibles;
6. **Implementado:** resultados alineados con metadatos de cadena/inicio,
   semilla y permutación.

Criterio de aceptación: un ejemplo pequeño con etiquetas permutadas y varios
inicios/cadenas produce alineamiento idempotente, diagnósticos reproducibles y
no combina cadenas antes de validarlas.

### P1.1. SVI para LRM

Es un requisito de uniformidad de la versión final del paquete, aunque no es un
requisito para ejecutar E1--E5. Está implementado con minibatches de
observaciones y actualización de parámetros naturales; no introduce grupos
ficticios para reutilizar el SVI jerárquico.

- **Implementado:** usar un control público `minibatch_size` y registrar
  `sampling_unit = "observation"`;
- **Implementado:** escalar `X'X`, `X'y`, `y'y` y el residuo esperado por
  `n / |B|`;
- **Implementado y probado:** hacer que batch completo y tasa uno reproduzcan
  la actualización CAVI;
- **Implementado:** conectar multistart, contratos, predicción y evaluación sin ramas
  especiales visibles al usuario;
- **Implementado:** demostrar reproducibilidad y estadísticas suficientes
  insesgadas. La utilidad computacional por tamaño se medirá en P5/piloto.

Este bloque no impide construir el `smoke` del paper, que evalúa SVI sobre
CHLRM, pero debe cerrarse antes de publicar la primera versión completa.

### P1.2. Nombre canónico CAVI

**Implementado.** `MCMC`, `CAVI` y `SVI` son los únicos identificadores de
método en contratos, resultados, diagnósticos, orquestación y consumidores.
“VI” se usa únicamente para la familia variacional o en nombres privados de
maquinaria compartida por CAVI y SVI; no es un alias de método aceptado. El
booleano `stochastic` permanece encapsulado en el adaptador hacia esos kernels
compartidos y no forma parte de la API pública.

## P2. Contratos y estado uniformes

- **Implementado:** migrar LRM/HLRM al contrato
  `initialized/running/fitted/failed` de CHLRM.
- **Implementado:** aislar RNG de inicialización, ajuste, evaluación y
  predicción en todas las clases.
- **Implementado en constructores:** validar datos, hiperparámetros, settings y
  matrices SPD antes de iniciar loops. Las operaciones intermedias de los
  kernels también fallan ante pérdida material de simetría o definitud.
- **Implementado para MCMC:** checkpoint tipado plano, guardado atómico y
  reanudación exacta sin mutar el ajuste de entrada; multichain conserva el
  estado individual de cada cadena y puede extenderse por bloques hasta
  convergencia o máximo.
  CAVI/SVI conserva estados iniciales y resultados, pero no reanuda loops
  interrumpidos.
- **Implementado:** sustituir `logdet()` e inversas explícitas por operaciones
  Cholesky que exigen SPD; la simetrización sólo corrige error de redondeo.
- **Implementado:** uniformar `fit()`, `summarize_posterior()`, `predict()`,
  `posterior_predict()`, `evaluate_in_sample_fit()`,
  `posterior_predictive_check()` y `predictive_criteria()` mediante contratos
  versionados, un núcleo común y nombres canónicos sólo para cantidades
  aplicables. No se exponen parámetros ficticios ni campos opcionales rellenos
  con `NULL`.
- **Implementado:** separar parámetros, variables latentes y objetivos de
  optimización/monitoreo; retirar nombres públicos sustituidos sin wrappers
  permanentes.
- **Implementado:** contrato `simulation-result-2` uniforme con nombres
  semánticos de verdad/configuración; se retiraron `sim_*`, `real_parameters`
  y `cluster_info` sin wrappers permanentes.
- **Implementado:** los presets del protocolo salieron de la API; sus valores
  concretos viven en `paper/evaluation/`.
- **Implementado:** política de nombres privados: los símbolos que coinciden
  con las ecuaciones (`Beta`, `Betaj`, `Betak`, `sigma2*`, `m`, `p`, `K`,
  `nj`) permanecen internos; `rho` designa sólo responsabilidades CHLRM y
  `learning_rate` el paso SVI. Se eliminaron aliases deterministas mutables
  `Betaj`, `sigma2j` y `nk` del estado Gibbs CHLRM.

Criterio de aceptación: reanudar produce el mismo resultado que una ejecución
continua bajo la misma secuencia RNG y los errores no dejan estados ambiguos.

### P2.0. Especificaciones canónicas de ajuste y convergencia

La fachada final reemplazará `hyperparams`, `run_settings`, `start_point` y el
booleano `stochastic` por cinco responsabilidades tipadas:

1. `model_spec`: familia LRM/HLRM/CHLRM, dimensión y likelihood;
2. `prior_spec`: distribuciones previas y sus hiperparámetros;
3. `initialization_spec`: estrategia o valores iniciales, independientes de la
   prior;
4. `inference_control`: presupuesto, ejecución y criterio de parada del motor;
5. `data_scaling_spec`: preprocesamiento explícito aprendido en entrenamiento.

La llamada canónica será `fit_model(model, data, prior, method,
initialization, control, scaling, seed)`. Los nombres internos de las
ecuaciones pueden conservarse en kernels, pero existe un solo traductor entre
nombres públicos semánticos y notación interna. No se aceptarán listas abiertas
con campos ignorados o valores que cambien de significado entre motores.

Controles finales:

- `mcmc_control(..., stopping = mcmc_stopping(...))`: ejecución/extensión por
  bloques, con R-hat rank-normalized, bulk/tail ESS y MCSE relativo después de
  un mínimo;
- `cavi_control(..., stopping = cavi_stopping(...))`: ELBO completo
  determinista más cambio de parámetros, paciencia y mínimo de iteraciones;
- `svi_control(..., learning_rate = power_schedule(...), stopping =
  svi_stopping(...))`: ELBO completo periódico, tendencia suavizada y cambio de
  parámetros naturales. El valor ruidoso de minibatch no detiene el algoritmo.

La fachada tipada produce `convergence-1` con `status`, `reason`, trabajo realizado,
métricas, umbrales e historial; sólo incluyen diagnósticos aplicables. El
stopper CAVI/SVI ya está conectado mediante controles tipados y exige, cuando
se solicita, estabilidad conjunta del objetivo y de coordenadas naturales del
estado variacional. MCMC multichain aplica trabajo mínimo, R-hat, bulk/tail ESS
y límite máximo; una ejecución de una sola cadena se declara honestamente como
presupuesto fijo, no como convergencia diagnosticada. SVI declara
`sampling_unit = "observation"` en LRM y `"group"` en HLRM/CHLRM.

Aceptación: las tres familias reciben el mismo tipo de prior por todos sus
motores; CAVI/SVI pueden comenzar desde estados pareados; ninguna
inicialización altera la prior; y cada método informa convergencia o agotamiento
del presupuesto sin traducir un diagnóstico ajeno.

### P2.1. Interoperabilidad y salidas idiomáticas de R

Los contratos versionados son el núcleo semántico, pero todavía falta una
superficie cómoda para usuarios de paquetes estadísticos. Antes de congelar el
namespace se debe:

1. **Implementado en su núcleo:** crear una fachada de ajuste S3 común para
   LRM/HLRM/CHLRM y mantener los motores mutables R6 como implementación;
2. **Implementado:** están registrados `print()`, `summary()`, `predict()`,
   `fitted()`, `residuals()`, `coef()` y `nobs()`; los accesores semánticos
   `population_parameters()`, `group_parameters()` y `group_deviations()`
   evitan imponer terminología de efectos fijos/aleatorios no equivalente;
3. **Implementado:** `posterior-draws-1` conserva procedencia MCMC, nombres
   semánticos y distingue muestras de `q`; `posterior::as_draws()` lo convierte
   a `draws_df` con nombres indexados y conserva la procedencia MCMC original
   en `source_iteration` y `source_stored_draw`;
4. **Implementado:** `as_yrep()` y `as_log_lik_matrix()` entregan draws ×
   observación/unidad, validan el target y preservan el estimando; PSIS-LOO
   reutiliza el segundo extractor;
5. **Implementado:** `components()`, `responsibilities()`,
   `classify_groups()`, `classification_uncertainty()` y `co_clustering()`
   exponen sólo en CHLRM las cantidades que realmente aplican;
6. **Implementado:** `tidy_parameters()` y `tidy_truth()` ofrecen
   representaciones largas para simulación, tablas y gráficos sin volver
   `dplyr` una dependencia del núcleo;
7. conservar la interfaz matricial y evaluar una interfaz `formula + data` que
   almacene términos, contrastes y niveles para predicción segura.

Se usará `posterior` en `Imports`; `loo` y `bayesplot` pueden vivir en
`Suggests`. No se añadirán `lme4`, `brms`, `mclust` o `mixtools` sólo para
imitar sus nombres. El criterio de aceptación es que un mismo ajuste alimente
resumen, predicción, PPC y LOO sin reconstruir parámetros ni perder la unidad
predictiva.

### P2.2. Modelo, likelihood, simulación y priors coherentes

El modelo, no el motor, es dueño de la likelihood. MCMC, CAVI y SVI consumen el
mismo kernel probabilístico y la misma prior; sólo cambia el algoritmo de
inferencia. `log_likelihood()` reutiliza ese modelo con el objetivo predictivo
declarado y no una fórmula duplicada dentro de GOF.

El paquete expone una operación `simulate_data()` con tres DGP tipados y
mantendrá la predictiva posterior fuera de ellos:

1. `conditional_dgp()` fija los parámetros directos de la likelihood y genera
   latentes no condicionadas y respuestas;
2. `hierarchical_dgp()` fija parámetros poblacionales y muestrea la jerarquía
   local completa; sólo aplica a HLRM/CHLRM;
3. `prior_predictive_dgp()` genera parámetros y datos desde la distribución
   conjunta; sirve para prior predictiva y como primitiva de SBC;
4. `posterior_predict()` genera respuestas futuras después del ajuste y no se
   reutiliza como generador Monte Carlo preajuste.

`truth` será siempre una salida realizada y no el nombre de una especificación
de entrada. Esto evita confundir parámetros fijos con una distribución que los
genera.

Entregables de simulación:

- **Implementado:** hacer obligatorios coeficientes/varianzas en `conditional_dgp()`; no
  activar una geometría editorial cuando falten;
- **Implementado:** permitir en CHLRM pesos poblacionales y asignaciones opcionales, conservando
  ambas cantidades por separado;
- **Implementado:** generar internamente `X` y `y` desde `simulation_design()` y tamaños/partición
  desde `simulation_sample()`;
- **Implementado:** uniformar soporte opcional de train/test sin inventar campos nulos;
- **Implementado:** devolver verdad/configuración suficientes y `tidy_truth()`;
- **Implementado:** registrar para cada bloque si fue condicionado, muestreado desde la jerarquía
  o muestreado desde la prior;
- **Implementado:** reutilizar la misma primitiva observacional que `posterior_predict()` y las mismas
  parametrizaciones que prior/likelihood, sin duplicar fórmulas;
- **Implementado en paquete/tests:** mover la geometría al paper, migrar tests y
  retirar los nombres provisionales.

La nueva procedencia cambió el significado del resultado: `simulation-result-2`
reemplazó al schema v1 en el paquete, tests y consumidores activos de tesis,
sin aliases ni wrappers históricos. La verdad observación por observación es
opt-in mediante `retain_observation_truth`; `max_output_mb` impide materializar
simulaciones que excedan el presupuesto explícito.

Pruebas obligatorias de los DGP:

- momentos condicionales de $y$ y momentos de parámetros locales contra sus
  fórmulas teóricas;
- identidad entre la densidad usada para generar y la evaluada por la
  likelihood bajo la verdad realizada;
- diseños fijos y aleatorios, RNG aislado, dimensiones degeneradas válidas y
  rechazo de soportes/covarianzas inválidos;
- train/test de las mismas unidades conserva sus latentes; test de unidades
  nuevas vuelve a muestrearlas desde el nivel poblacional;
- el modo prior-predictivo genera pares $(\theta,y)$ desde la misma prior y
  likelihood que recibe el ajuste y soporta diagnósticos SBC, con cantidades
  invariantes a etiquetas para CHLRM.

El paquete genera un dataset por llamada. Réplicas, grillas factoriales, MCSE y
criterios ADEMP pertenecen a la suite del paper, que debe justificar cada valor
del DGP en vez de usar defaults editoriales ocultos.

Entregables de priors:

- **Implementado:** crear `lrm_prior()`, `hlrm_prior()` y `chlrm_prior()` con argumentos públicos
  descriptivos y traducción única a la notación interna;
- **Implementado/decidido:** separar prior explícita e `initialization_ols()`;
  no se expondrá `empirical_prior()` en v1 y ninguna opción usa un rótulo
  ambiguo como `NON_INFORMATIVE`;
- **Implementado:** usar el mismo objeto prior validado para MCMC, CAVI y SVI;
- **Implementado:** comprobar que soporte, unidades y convenciones `shape/rate` coinciden entre
  generador, prior y kernels;
- **Implementado:** impedir que un simulador copie automáticamente `truth` dentro de la prior.

Política de defaults: los constructores explícitos no estimarán escalas desde
los datos. Antes de la versión 1.0 se exigirá prior explícita. Sólo se añadirá
una receta regularizadora general si puede describirse en una escala conocida,
mostrar sus momentos y superar simulación prior-predictiva; nunca se llamará
“no informativa”. Las priors OLS serán siempre opt-in y se identificarán como
empírico-Bayes. La conjugación requerida por Gibbs/CAVI/SVI se conservará salvo
una decisión teórica explícita de cambiar los tres motores.

`prior_moments()` implementa la auditoría de escala sin simular: reporta la
existencia de medias/varianzas inverse-gamma, la media inverse-Wishart cuando
existe, momentos gamma/Dirichlet y los momentos marginales locales de la
jerarquía. La prior predictiva se construye directamente con el DGP tipado, sin
un wrapper `prior_predict()` adicional.

#### Migración de la auditoría

| Implementación actual | Clasificación correcta | Destino |
|---|---|---|
| `simulate_lrm()` | DGP condicional LRM | retirado; reemplazado por `conditional_dgp()` |
| `simulate_hlrm()` | DGP condicional HLRM | retirado; reemplazado por `conditional_dgp()` |
| `simulate_hlrm_hierarchy()` | DGP jerárquico HLRM | retirado; reemplazado por `hierarchical_dgp()` |
| `simulate_chlrm_protocol()` | DGP condicional CHLRM más receta editorial | retirado; núcleo en paquete y receta en `paper/evaluation/` |
| `simulate_chlrm_hierarchy()` | DGP jerárquico CHLRM | retirado; reemplazado por `hierarchical_dgp()` |

El DGP prior-predictivo ya genera la jerarquía desde el mismo `prior_spec` que
recibe el ajuste. La orquestación y los diagnósticos de rango SBC siguen siendo
tarea de la suite científica.

- Los generadores históricos usaban diseños inconsistentes. Fueron reemplazados
  por `simulation_design()`, capaz de
  representar grilla con jitter, normal independiente, normal multivariada,
  binaria y diseño suministrado, con nombres de columnas e intercepto
  explícitos.
- El contrato común distingue test de los mismos grupos —misma verdad local/asignación,
  nuevas observaciones— de test con grupos nuevos —nuevos efectos y, en CHLRM,
  nuevas asignaciones— para HLRM y CHLRM, y observaciones nuevas para LRM.
- La jerarquía generadora HLRM/CHLRM sí coincide en forma con los kernels:
  coeficientes normales, covarianza inverse-Wishart, varianzas inverse-gamma y
  escala gamma. Se conservará esa correspondencia y se probará por momentos.
- `hyperparams_*(type = "NON_INFORMATIVE")` no es independiente de los datos:
  usa residuos OLS para fijar escalas. Se eliminará ese rótulo y se clasificará
  como empírico-Bayes.
- `Betaj_ols` viaja hoy dentro del objeto de hiperparámetros aunque sólo debe
  controlar inicialización; se moverá a un objeto de valores iniciales.
- `sigma20` aparece en el constructor HLRM como escala auxiliar de OLS y se lee
  durante inicialización, pero no define la prior jerárquica. Se retirará del
  contrato de prior HLRM/CHLRM.
- Con `nu0 = 2`, usado por constructores actuales, la media prior condicional de
  la varianza residual no existe. Esto puede ser una prior propia de cola
  pesada, pero no se describirá como “varianza media”; los presets deberán
  declarar qué cuantil o momento calibran.
- `n0 = p + 2` hace que la media inverse-Wishart sea `S0` bajo la
  parametrización usada. Esta convención es coherente y se documentará junto a
  las condiciones de existencia de media y varianza.
- Constantes como `jitter = 1e-6`, `nu0 = 2`, `a0 = 2`, `n0 = p + 2` y
  `alpha0 = 1` están hoy incrustadas en `hyperparams_*()`. Sólo el jitter
  numérico podrá tener un default técnico documentado; las demás serán
  argumentos o presets generales con significado estadístico explícito.

#### Auditoría y consolidación de likelihoods

- LRM conserva likelihood normal por observación; agrupar filas sólo cambia la
  unidad predictiva, no introduce latentes.
- HLRM conserva dos estimandos válidos: condicional a parámetros locales para
  grupos observados y marginal sobre coeficientes/varianzas locales para grupos
  nuevos. La API elegirá el estimando a partir del objetivo predictivo y lo
  registrará, sin pedir al usuario ordinario que seleccione detalles internos.
- CHLRM usa una única asignación por grupo. El scoring principal integra esa
  asignación mediante `log-sum-exp`; la likelihood condicional usada por Gibbs
  seguirá siendo privada y no se confundirá con WAIC marginal.
- `loglik_pointwise_chlrm()` reemplaza silenciosamente una varianza inválida por
  `1e-12`. Esa corrección se retirará: el contrato debe rechazar varianzas no
  positivas. En las SSE sólo se truncarán negativos atribuibles a redondeo
  después de comprobar una tolerancia relativa; un negativo material será un
  error.
- Se retirarán helpers públicos/duplicados `loglik_pointwise_*` y
  `loglik_total_*` cuando los kernels consuman nombres privados inequívocos.
- Se conservarán las pruebas ya existentes de agregación pointwise, cuadratura
  HLRM e invariancia de etiquetas. Se añadirán igualdad de likelihood entre
  motores dados los mismos parámetros, límite CHLRM `K=1`, límite HLRM sin
  heterogeneidad y acuerdo simulador--likelihood bajo parámetros verdaderos.

Los simuladores devuelven datos y verdad en la escala del DGP. El escalado
opcional se aplica después mediante el contrato explícito de P2.3. Aceptación:
pruebas de RNG, dimensiones, momentos del diseño, ruido condicional,
serialización de la verdad, igualdad del prior entre motores y rechazo de
parametrizaciones ambiguas.

### P2.3. Escalado explícito y reversible

El paquete ofrecerá `data_scaling()` reusable por LRM, HLRM y CHLRM y aplicado
de forma idéntica antes de MCMC, CAVI o SVI. No se usará el nombre
`vb_standardizer` ni se agregará un booleano suelto a cada motor.

- **Implementado:** `fit_scaling()` aprende centros/escalas sólo de entrenamiento,
  excluye el intercepto y exige declarar qué columnas se escalan;
- **Implementado:** `transform_data()` transforma matrices, datos agrupados y
  respuesta con orden/nombres validados;
- **Implementado:** `transform_prior()` y la frontera de `fit_model()` llevan
  priors e inicializaciones declaradas en escala de entrada a la escala del
  backend;
- **Implementado:** parámetros, draws y resúmenes restauran interceptos,
  pendientes, varianzas y covarianzas poblacionales/locales/de componente;
- **Implementado:** predicción, predictiva posterior, PPC y ajuste dentro de
  muestra vuelven a la escala de entrada;
- **Implementado:** `Lambda0`, `S0`, la escala/rate asociado a $\xi^2$ y los
  valores iniciales usan las mismas identidades afines; y
- **Implementado:** la log-likelihood aplica el jacobiano por observación o
  grupo y WAIC se recalcula desde la matriz restaurada.

La fachada `fit_model()` puede orquestar esas operaciones cuando recibe
`scaling = data_scaling(...)`, pero conserva el mapa y reporta la escala; no es
un cambio silencioso. Priors públicas, resúmenes y predicciones se expresan por
defecto en la escala de entrada. La restauración comparte un motor afín y usa
métodos por modelo sólo para enumerar bloques: LRM restaura coeficientes y
varianza residual; HLRM añade parámetros locales, covarianza y $\xi^2$; CHLRM
añade componentes y deja invariantes pesos/asignaciones/responsabilidades.

El objeto debe ser serializable, guardar la procedencia de entrenamiento y
fallar ante desviación cero, columnas ausentes/reordenadas o transformación
incompleta. La prueba crítica de CV ajusta un objeto distinto por fold y
demuestra ausencia de fuga.

## P3. Evaluación predictiva general

- **Implementado:** contrato `log-likelihood-1` único para LRM, HLRM y CHLRM,
  con matriz pointwise, identificadores de unidad, condicionamiento y variables
  integradas. `predictive_criteria()`/WAIC consumen este contrato; los escalares
  internos del Gibbs permanecen separados.
- **Implementado:** `predictive_log_likelihood()` y `predictive_score()`
  calculan log density pointwise, ELPD y MLPD sobre respuestas retenidas.
- **Implementado:** LOGO/CV por grupos nuevos, CV configurable dentro de grupo
  y K-fold por observación mediante planes reproducibles; el runner secuencial
  limita `max_fits`, aprende escalado dentro de cada fold y puede devolver
  log scores pointwise junto con RMSE, MAE, cobertura y ancho predictivo.
- **Implementado:** `psis_loo()` usa el paquete de referencia `loo`, exige
  draws MCMC y convergencia multichain por defecto, estima eficiencia relativa
  con cadenas balanceadas y conserva Pareto-k. Para CAVI/SVI se exige CV
  refitted en vez de presentar draws de `q` como una posterior exacta.
- **Implementado:** PPC streaming con discrepancias configurables; el preset
  estándar calcula T1--T5 cuando son aplicables. Reporta MCSE nominal y, para
  MCMC, MCSE basado en el ESS de cada secuencia indicadora por cadena.
- **Implementado:** `predictive_performance()` reúne RMSE, MAE, R² de test,
  cobertura y ancho de intervalos predictivos, ELPD y MLPD sin confundirlos
  con ajuste dentro de muestra o PPC.
- **Implementado:** se retiraron DIC, `include_dic`, sus puntos posteriores,
  campos de resultado y rutas de reporting; todos los consumidores usan
  `predictive-criteria-2` sin wrappers permanentes.
- **Implementado:** resultados pointwise preservados para agregación externa.
  `predictive_criteria(return_log_likelihood = TRUE)` los devuelve bajo el
  contrato `log-likelihood-1`; los runners deciden si los persisten para no
  inflar todos los objetos ajustados.

### P3.1. Retiro completo de DIC

Esta migración es parte obligatoria de P3 y precede cualquier congelamiento de
la API. No consiste únicamente en ocultar el valor o devolver `NULL`:

1. retirar `include_dic` de las firmas públicas de LRM, HLRM y CHLRM para todos
   sus motores;
2. eliminar los cálculos y helpers DIC de `gaussian-evaluation.R`,
   `hlrm-marginal.R` y `chlrm-evaluation.R`;
3. retirar `dic_conditional`, `dic_group_marginal`, puntos posteriores,
   estados y metadatos DIC de los objetos de resultado;
4. publicar únicamente el contrato `predictive-criteria-2` y migrar sus
   validadores, impresores y consumidores;
5. eliminar DIC de reporting, runners activos y tablas del paper; los archivos
   históricos de `thesis/v0` no se convierten en dependencias del paquete;
6. sustituir pruebas DIC por pruebas de WAIC, log-likelihood pointwise y
   equivalencia del estimando marginal; y
7. verificar con una búsqueda final que `R/`, `tests/` y los runners activos no
   contengan rutas ejecutables DIC ni aliases de compatibilidad.

Estado: **implementado**. Las seis clases, helpers, contrato, reporting y
pruebas fueron migrados conjuntamente. La suite actual ejecuta 150 bloques sin
fallos, errores, warnings ni skips.

Criterio de aceptación: las pruebas completas pasan; ninguna API acepta
`include_dic`; ningún resultado contiene campos DIC; WAIC conserva unidad,
objetivo, condicionamiento y procedencia explícitos; y las comparaciones
principales usan WAIC junto con PSIS-LOO o CV/LOGO y MLPD cuando corresponda.

Pertenecen al paquete WAIC/PPC/LOGO genéricos, T1--T5 configurables y el scoring
sobre respuestas observadas. No pertenecen al paquete el sesgo o cobertura
agregados entre réplicas, 5000 draws, folds concretos ni comparaciones contra
el MCMC de referencia: ésos son decisiones y métricas de
`paper/evaluation/`.

## P4. Arquitectura instalable

- **Implementado:** crear `DESCRIPTION`, `NAMESPACE`, licencia y documentación inicial.
- **Implementado:** eliminar `source()` y `library()` de archivos internos y
  aplanar `R/` con nombres por responsabilidad.
- **Implementado:** exponer únicamente constructores genéricos; mover presets y geometrías
  específicas del paper a `paper/evaluation/`.
- **Implementado:** exportar el preprocesamiento reversible de P2.3, pero no argumentos booleanos
  de estandarización dentro de los modelos.
- **Implementado:** mantener kernels puros internos y una API exportada pequeña.
- **Implementado:** declarar dependencias mínimas en `Imports`; herramientas de evaluación en
  `Suggests` cuando corresponda.
- Añadir vignettes pequeñas de simulación, ajuste, diagnóstico y predicción.
- **Implementado para el paquete:** hacer que tests y ejemplos sean
  independientes del directorio de trabajo; la suite del paper tiene su propio
  runner desde la raíz y consume sólo exports.

Criterio de aceptación alcanzado para la arquitectura: `R CMD check` pasa sin
errores, warnings ni notes y el smoke del paper instala el paquete en una
biblioteca temporal antes de probarlo. Las vignettes siguen siendo una mejora
documental previa a publicación, no una dependencia oculta del paper.

## P5. Rendimiento seguro

- **Implementado:** `estimate_fit_resources()` estima por dimensiones la
  memoria pico conservadora, el estado retenido y el almacenamiento potencial
  de evaluación; `fit_model()` rechaza solicitudes que exceden
  `max_memory_mb` antes de construir backends.
- **Implementado:** el mismo límite se acepta en MCMC/CAVI/SVI con una o varias
  cadenas/inicios y el estimado queda en `fit$metadata$resource_estimate`.
- **Implementado en smoke E1:** el manifiesto registra tiempo de ajuste,
  estimación de memoria, presupuesto y tamaño real del RDS.
- **Implementado en preflight E3:** los controles confirmatorios se
  materializan en un plan; `screen` acota el perfilado a 10 tareas y
  `factorial` declara 109, ambos secuenciales y bloqueados hasta aprobación.
- Perfilar kernels y RSS con el `screen` antes de optimizarlos.
- Procesar draws/log-likelihood por bloques cuando sea necesario.
- Evitar tensores observación × grupo × componente × draw.
- Documentar qué límites controlan arrays locales y cuáles no acotan RAM total.

El scheduler de miles de tareas no pertenece al paquete; sí pertenecen los
ajustes individuales eficientes, serializables y predecibles en recursos.

## Orden vigente

P1 está cerrado para los motores existentes. P2 está cerrado para estado de
ejecución y resultados aplicables, pero no para construcción pública.
P2.0--P2.3 ya cerraron la superficie tipada, las especificaciones de
modelo/prior/inicialización/control, los generadores y el escalado reversible.
P3 separa ajuste descriptivo, PPC, criterios y predicción, y conecta el
contrato de log-likelihood con scoring retenido, LOGO/CV y PSIS-LOO. I5/P1.1
cerró SVI para LRM. Después se retiran los loaders `source()` y se completa la
estructura instalable de P4; así cada cambio físico queda protegido por pruebas
semánticas ya estables.

El trabajo se ejecutará en incrementos atómicos. Un incremento no se cierra si
deja dos APIs públicas, dos versiones activas del mismo schema o consumidores
que dependan del nombre retirado. Los adaptadores privados hacia los kernels
actuales pueden permanecer cuando sean la frontera deliberada de la
arquitectura, pero no se publican como compatibilidad legacy.

1. **I0 — P3.1, retirar DIC (cerrado).** Se eliminaron firmas, cálculos, campos,
   tests y consumidores DIC; `predictive-criteria-2` conserva WAIC y la
   log-likelihood pointwise con su unidad/objetivo explícitos.
2. **I1 — P2.0, construcción tipada (cerrado).** Se implementaron `model_spec`,
   `prior_spec`, `initialization_spec` e `inference_control`, el traductor
   privado a los kernels y `fit-result-2`. `CAVI` es canónico hasta los
   contratos internos; los nombres privados `vi_*` identifican código
   variacional compartido por CAVI y SVI, no una API anterior.
3. **I2 — P2.2, modelo y simulación (cerrado).** Se consolidaron las
   primitivas de prior, likelihood y RNG; los tres DGP, `simulate_data()`,
   `simulation-result-2`, pruebas de momentos/densidad/unidades y el smoke SBC
   ya usan el contrato nuevo sin simuladores públicos sustituidos.
4. **I3 — P2.3, escalado reversible (cerrado).** `data_scaling()` y las
   transformaciones de datos, priors, inicializaciones, parámetros,
   predicciones, PPC y densidades se integraron a `fit_model()` con mapa
   aprendido sólo sobre entrenamiento.
5. **I4 — resto de P3 (cerrado).** PSIS-LOO MCMC, LOGO/CV, MLPD, métricas de
   error/cobertura e intervalos y PPC configurable usan contratos explícitos.
   El PPC diferencia MCSE IID de MCSE corregido por ESS para MCMC.
6. **I5 — P1.1, SVI para LRM (cerrado).** El motor usa minibatches de
   observaciones, parámetros naturales y equivalencia CAVI con batch completo.
7. **I6 — P1.2/P2.1/P4, superficie idiomática e instalación (en curso).** La
   estructura plana, instalación, `R CMD check`, fachada S3, métodos base,
   draws `posterior`, extractores orientados y accesores CHLRM ya están
   implementados. Los consumidores activos de tesis/paper usan sólo exports y
   una guarda estática verifica esa frontera. `fit_orchestration()` ya devuelve
   una auditoría sanitizada y no el backend; CHLRM conserva clasificaciones por
   inicio para estabilidad. Faltan vignettes y decidir si una interfaz adicional
   `formula + data` justifica su costo de mantenimiento.
8. **I7 — P5, presupuesto seguro (en curso).** La guarda pública y uniforme de
   memoria y el registro por tarea smoke están implementados. E3 ya midió un
   máximo de 496 MiB RSS en `screen`; su factorial piloto completó 109/109
   tareas en 17.1 minutos, con 488.7 MiB RSS máximo y 143.9 MiB de resultados.
   Falta extender el profiling a E1/E2/E4/E5 y optimizar sólo los kernels que
   el perfil identifique.

P5 es transversal: cada incremento registra pruebas y costo local; antes del
`pilot` se perfila la API congelada y se fijan límites seguros. La suite del
paper puede empezar su especificación E-A durante I1--I4, pero no ejecutará
ajustes confirmatorios hasta completar I6, instalar el paquete y aprobar la
Puerta P de `paper/EVALUATION_PLAN.md`.

I1/P2.0 está cerrado. Existe la construcción tipada completa:
especificaciones tipadas de modelo/prior/inicialización/control, un único
adaptador privado, `fit_model()` y `fit-result-2` para las nueve combinaciones
implementadas. CAVI/SVI puede detenerse por ELBO y cambio de parámetros
naturales; CHLRM-MCMC materializa correctamente la inicialización k-means; y
los valores explícitos usan nombres públicos semánticos y validación de soporte.
La fachada integra MCMC multichain, CAVI/SVI multistart, `convergence-1` y las
operaciones públicas de resumen, ajuste dentro de muestra, PPC, criterios y
predicción. `posterior_draws()` ya entrega todos los parámetros aplicables con
nombres semánticos, procedencia cadena/iteración para MCMC y muestras de (q)
identificadas para CAVI/SVI; las asignaciones CHLRM son latentes opcionales.

La migración de consumidores se cerró por dueño estadístico, no mediante una
sustitución textual que conservara decisiones retiradas:

| Consumidor | Momento de migración | Razón |
|---|---|---|
| tests unitarios de ajustes y operaciones | I1, implementado | fijan el contrato sin depender del paper |
| `runs/01-farms-hlrm` y `runs/02-farms-chlrm` | implementado | priors explícitas, escalado y fachada tipada |
| `runs/03-simulation` | implementado | DGP tipados y verdad versionada |
| tuning SVI y selección de K | implementado | controles tipados, WAIC/LOGO y perfil smoke |
| recetas del artículo | implementado en `smoke` y piloto E3 | E1--E5 consumen el paquete instalado; screen y factorial E3 midieron recursos y produjeron agregados/figuras; `full` y productos confirmatorios siguen siendo responsabilidad de la suite |

`dev/check-consumer-boundary.R` rechaza carga directa de kernels, `:::`,
constructores históricos de hiperparámetros y orquestadores privados desde
`paper/evaluation` o `thesis/final/code`. `thesis/v0` permanece fuera de esta
regla por ser un artefacto histórico congelado.
