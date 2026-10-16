# Trazabilidad del protocolo computacional frente al repositorio

## Propósito

Este documento contrasta los requisitos de
[`Computacion_revision_paper_Camilo.pdf`](references/Computacion_revision_paper_Camilo.pdf) con:

- el código principal de `R/`;
- los análisis exploratorios históricos de `thesis/final/code/`;
- la suite confirmatoria que se construirá en `paper/evaluation/`;
- las pruebas en `tests/`;
- el contrato público `API_CONTRACT.md`, el roadmap maestro `ROADMAP.md` y el
  plan del paquete `PACKAGE_ROADMAP.md`.

La evaluación corresponde al estado local revisado el 29 de septiembre de
2026. No certifica resultados numéricos todavía: identifica qué existe, qué es
parcial, qué es incorrecto para el protocolo y qué falta implementar.

El PDF usa “VI” para el algoritmo variacional determinista. En este documento
se conserva esa etiqueta al describir literalmente sus requisitos o el código
actual; la API final lo registra como **CAVI**, mientras **VI** queda como la
familia que contiene CAVI y SVI.

## Marco académico de simulación

La suite se documentará mediante ADEMP: *Aims, Data-generating mechanisms,
Estimands, Methods* y *Performance measures*. Esta estructura evita confundir
el simulador reusable con el experimento que lo llama:

| ADEMP | Propiedad en este repositorio |
|---|---|
| Aims | preguntas E1--E5 y comparación de exactitud, incertidumbre, clustering, predicción y costo |
| Data-generating mechanisms | constructores genéricos del paquete más configuraciones S1--S7c de esta suite |
| Estimands | parámetros poblacionales/de componente, partición, cantidades predictivas y criterios con unidad declarada |
| Methods | MCMC, VI y SVI con configuración, inicialización y reglas de parada registradas |
| Performance measures | sesgo, RMSE, cobertura, ancho, fidelidad frente a MCMC, ARI/CPE, MLPD, runtime, memoria y sus MCSE |

Los resultados por réplica se guardarán en formato largo para poder resumir
métricas escalares y sus errores Monte Carlo con herramientas como `rsimsum`.
Las métricas de clustering, fidelidad distribucional y predicción jerárquica
seguirán siendo específicas de esta suite. No es necesario que el paquete
dependa de un framework de simulación: la separación
generar--ajustar--resumir y el scheduler reanudable pertenecen a
`paper/evaluation/`.

El artículo de Mehrotra y Maity sobre el paquete `vir` es un precedente
complementario para la comparación Gibbs--CAVI--SVI. Evalúa MSE de
coeficientes, cobertura, error predictivo y tiempo al variar el tamaño de los
datos, y demuestra por qué la calidad puntual no sustituye la calibración de
incertidumbre. Se adopta esa lógica, no su API de funciones por combinación ni
su diseño completo: usa 50 datasets para sus simulaciones, sólo cinco para los
tiempos y en algunos benchmarks fija iteraciones SVI. Nuestro protocolo debe
mantener escenarios pareados, MCSE, reglas de convergencia y presupuesto de
memoria explícitos. Véase [Mehrotra y Maity (2021)](https://arxiv.org/abs/2102.08877).

En particular, el artículo presenta las condiciones de Robbins--Monro
$\sum_t\rho_t=\infty$ y $\sum_t\rho_t^2<\infty$, pero también ofrece una tasa
constante elegida empíricamente. Aquí la familia decreciente
$\rho_t=(t+\tau)^{-\kappa}$, $0.5<\kappa\leq1$, será la configuración
confirmatoria; las tasas constantes podrán aparecer sólo como comparador de
sensibilidad finita y nunca como garantía asintótica o default universal.

## Alcance y propiedad

Este documento especifica **qué debe ejecutar y demostrar la evaluación**, no
la superficie pública del paquete. La separación vinculante es:

- el paquete recibe datos y parámetros arbitrarios y devuelve simulaciones,
  ajustes, draws, predicciones, likelihoods y diagnósticos estructurados;
- esta suite fija S1--S7c, priors principal/sensibilidad, réplicas, semillas,
  cadenas, inicios, folds, tamaños de Monte Carlo y reglas de aceptación;
- esta suite calcula métricas contra verdad simulada o contra un método de
  referencia, administra archivos y produce tablas/figuras;
- el manuscrito importa artefactos agregados y nunca vuelve a ajustar modelos.

Una brecha del paquete se registra aquí sólo como **dependencia externa**. Su
implementación y pruebas se planifican en `PACKAGE_ROADMAP.md`; no se copian
kernels dentro de `paper/evaluation/`.

## Relación con el plan de refactor

El núcleo separa inferencia, ajuste descriptivo, chequeo predictivo posterior,
criterios predictivos y predicción. `fit()`, `summarize_posterior()`,
`evaluate_in_sample_fit()`, `posterior_predictive_check()`,
`predictive_criteria()`, `predict()` y `posterior_predict()` usan contratos
versionados y sólo incluyen cantidades aplicables. HLRM dispone de likelihood
marginal numérica por finca, con integración de efectos locales, draws globales
MCMC/q y controles de cuadratura. En LRM agrupar observaciones cambia la unidad
de scoring, pero no crea una marginalización inexistente. El núcleo PPC streaming
es común, con replicados observacionales para las tres familias y métodos
implementados. El preset estándar ya calcula T1--T5 cuando
son aplicables y corrige el MCSE mediante ESS de la secuencia indicadora MCMC.
Esto no sustituye los estudios E1--E5 ni fija todavía sus 5000 replicados.
La API predictiva ya separa media y respuesta, grupos observados y nuevos, y
MCMC de q mediante metadatos. `posterior_predict()` genera `y_rep` reutilizable
por `pp_check_from_prediction()`; `posterior_predictive_check()` conserva una
ruta streaming basada en la misma primitiva para limitar memoria. La API
también calcula ELPD/MLPD retenido, error/cobertura de test y CV/LOGO genérico;
cobertura Monte Carlo entre réplicas, error agregado y regret siguen siendo
métricas de esta suite.
`predictive_criteria(return_log_likelihood = TRUE)` conserva la matriz
pointwise exacta usada por WAIC y validación predictiva. Los runners activos la escriben comprimida
para permitir auditoría y métricas posteriores sin repetir el ajuste.
Los seis constructores de inferencia validan ahora los hiperparámetros y las
matrices SPD antes de los loops; los kernels usan Cholesky para sistemas y
log-determinantes y rechazan asimetrías que no sean mero redondeo. Esto reduce
fallos silenciosos en las réplicas, pero no sustituye el catálogo de escenarios,
las métricas ni la orquestación E1--E5.

Este documento es la **especificación de aceptación** del trabajo experimental.
`PACKAGE_ROADMAP.md` es el **plan de implementación reusable** y `ROADMAP.md`
fija sus dependencias con este protocolo. Deben evolucionar de
forma coordinada:

1. cada requisito de esta revisión debe apuntar a un módulo, una prueba y una
   salida reproducible;
2. cada fase del plan debe cerrar una o más filas de la matriz de trazabilidad;
3. una fila sólo cambia a `Implementado` cuando pasa su criterio de aceptación;
4. los experimentos definitivos no comienzan mientras existan brechas
   `Incorrecto` en el estimando que utilizarán;
5. los resultados `smoke` y `pilot` validan el pipeline, pero no se presentan
   como resultados confirmatorios.

La meta no es escribir un script monolítico para replicar una tabla. La suite
debe declarar y ejecutar el protocolo mediante configuración, consumiendo
componentes generales del paquete.

El protocolo se interpreta además como validación científica del método, no
sólo como una lista de exigencias editoriales. Recuperación paramétrica,
cobertura, fidelidad posterior, clustering, predicción, convergencia y
escalabilidad son dimensiones que una implementación publicable debe poder
evaluar. El paquete R contendrá las capacidades reutilizables; E1--E5 serán una
suite externa de validación que use únicamente su API pública. Así se evita
mezclar lógica estadística con código específico del manuscrito sin perder la
posibilidad de reproducirlo.

## Leyenda

| Estado | Significado |
|---|---|
| Implementado | Existe y su contrato coincide sustancialmente con el PDF |
| Parcial | Hay una pieza reutilizable, pero faltan casos, orquestación o validación |
| Incorrecto | Existe una salida con ese nombre, pero calcula otro estimando |
| Ausente | No se encontró implementación reproducible |
| No aplica | El requisito era opcional y se decidió no conservarlo |

## Veredicto ejecutivo

El repositorio ya contiene un paquete instalable para LRM, HLRM y CHLRM con
MCMC, CAVI y SVI, más una suite separada que ejecuta E1--E5 en perfil `smoke`.
Eso valida contratos, persistencia, métricas y conexión entre piezas, pero aún
no constituye el estudio confirmatorio solicitado por el protocolo.

Las brechas bloqueantes son:

1. la prior exacta, el generador y el catálogo S1--S7c ya son usados por los
   runners smoke; cada manifiesto registra huellas de configuración y código,
   versión del paquete, commit y estado del árbol;
2. E1--E5 ya poseen expansión de tareas y runners reanudables en `smoke`;
   `pilot/full` permanecerán bloqueados hasta congelar controles y presupuestos;
3. MCMC administra cadenas independientes, extensión exacta por bloques,
   checkpoints, diagnósticos modernos, parada automática y referencia común
   CHLRM; falta fijar los tamaños exactos en configuración del paper;
4. VI y SVI ya tienen multistart con estados iniciales compartibles, registro de
   todos los intentos, selección por ELBO y referencia CHLRM reutilizable; falta
   fijar cinco inicios y configuraciones exactas en el runner del paper;
5. la parada VI/SVI es uniforme en los tres modelos disponibles y multistart
   conserva sus diagnósticos; falta materializar la configuración del paper;
6. las métricas puras de recuperación, cobertura, fidelidad posterior,
   clustering, selección y predicción existen; falta completar su agregación
   Monte Carlo entre réplicas y convertirla en tablas finales;
7. `predictive_criteria()` CHLRM ya usa WAIC marginal y DIC fue retirado; E5
   smoke integra WAIC y CV, pero faltan los presupuestos confirmatorios y sus
   diagnósticos de precisión Monte Carlo;
8. los perfiles confirmatorios todavía no han sido ejecutados con las réplicas,
   iteraciones y draws exigidos;
9. la reanudación por tarea funciona, pero faltan figuras/tablas automáticas,
   manifiesto de entorno completo y auditoría final de procedencia.

Por tanto, las salidas actuales deben considerarse exploratorias. Son útiles
para verificar que los algoritmos corren y para orientar el refactor, pero no
responden todavía las quince preguntas empíricas de la sección 30.

No se mantiene un porcentaje a partir del número de filas: una capacidad de
orquestación puede bloquear varias salidas mientras una métrica aislada puede
ser pequeña. El estado y la acción pendiente de cada requisito se registran en
la matriz siguiente y deben actualizarse junto con el código.

## Matriz de trazabilidad global

| PDF | Requisito | Código actual | Plan actual | Estado | Acción necesaria |
|---|---|---|---|---|---|
| 1 | Un pipeline común para MCMC, VI y SVI | Las clases comparten contratos versionados y E1--E5 tienen runners smoke sobre el paquete instalado y con huella verificable | Paquete común más runner declarativo | Parcial | Añadir un scheduler único para `pilot/full` y un manifiesto global |
| 2 | Generador CHLRM común | `simulate_data(generator = conditional_dgp(...))` produce train/test, verdad, condicionamiento y procedencia; E1--E4 lo consumen sin defaults del review en el paquete | DGP genérico del paquete más escenarios en esta suite | Implementado | Mantener los valores del protocolo exclusivamente en esta suite |
| 2.1 | Coeficientes con separación mínima controlada por `d` | `make_review_component_coefficients()` vive en `paper/evaluation/R/scenarios.R`, entrega parámetros al DGP y está probado para K=2,3,6,10 | Receta exclusiva de la suite | Implementado | K=1, cuando aplique fuera de E1, usa coeficientes explícitos |
| 3 | Prior independiente: `mu0=0`, `Lambda0=100 I`, `n0=p+2`, `S0=I`, `nu0=a0=b0=1`, `alpha0=1/K` | `review_chlrm_prior()` configura `chlrm_prior()` desde la suite; la API genérica no conoce el preset y la configuración queda incluida en la huella del run | Constructor genérico en paquete y receta aquí | Implementado | Conservar el preset versionado dentro de la suite |
| 4 | E1 con nueve escenarios (S1--S6 y S7a--S7c) y 200 réplicas | El catálogo fija el DGP; `expand_e1_tasks()` declara 5400 tareas full y `run_e1_smoke()` prueba persistencia, reanudación, recuperación alineada, clustering y predicción | Configuración declarativa más runner por perfiles | Parcial | Congelar controles en piloto y habilitar después las 200 réplicas |
| 4.1 | Test independiente con mismos grupos/asignaciones | `simulation_sample(..., test_group_mode = "same")` conserva asignaciones y genera nuevas covariables/respuestas; los runners simulados lo puntúan fuera de muestra | DGP tipado | Implementado | No confundir este test con grupos nuevos |
| 5 | Cuatro cadenas MCMC independientes | `fit_model()` con `mcmc_control(chains=...)` ejecuta cadenas configurables y `fit_orchestration()` conserva semillas, tiempos y diagnóstico | P1 del paquete | Implementado | El runner del paper debe fijar cuatro cadenas |
| 5.1 | Inicialización desde priors, sin verdad | MCMC inicializa estocásticamente desde hiperparámetros | Contrato RNG general | Parcial | Validar independencia y registrar semillas/estados iniciales |
| 5.2 | Warm-up 10000; bloques 5000; sin thinning; máximo 100000 | El runner exige `thin=1`, extiende exactamente por `block_size`, diagnostica cada bloque y corta en `max_samples`; checkpoint/reanudación está probado | P2 del paquete | Parcial | Fijar 10000/5000/100000 en la configuración del paper y persistir checkpoints por tarea |
| 5.2 | R-hat < 1.01, bulk/tail ESS >= 1000 y ESS/sec | La orquestación interna de `fit_model()` usa `posterior`, conserva resultados por variable y aplica esos umbrales por defecto | P1 del paquete | Implementado | Consumir el diagnóstico versionado desde el runner del paper |
| 5.3 | Monitorear componentes relabelados e invariantes | Se guardan logjoint, ocupación y PSM | Contempla invariantes | Parcial | Construir lista canónica de variables y diagnósticos por cadena |
| 6 | Relabeling iterativo con draw de mayor log joint y escala combinada | El orquestador aprende una sola referencia sobre todas las cadenas y las alinea antes del diagnóstico | P1 del paquete | Implementado | Usar el objeto de referencia producido por el ajuste del dataset |
| 6 | Aplicar referencia común a MCMC, VI y SVI | Un objeto `chlrm_label_reference` alinea draws MCMC y resultados VI/SVI; ambas clases lo exponen | P1 del paquete | Implementado | El runner multistart/multichain debe propagar una única referencia por dataset |
| 7 | Cinco inicializaciones VI con responsabilidades Dirichlet | `fit_model()` con `starts=5` materializa estados independientes, conserva cada intento y registra semillas mediante la orquestación interna | P1 del paquete | Implementado | El runner del paper debe fijar `starts=5` |
| 7 | Parada relativa `1e-8` durante cinco iteraciones; máximo 10000 | Control CAVI y diagnósticos implementados; el preflight E3 fija 10000, tolerancia y paciencia exactas | Segundo incremento | Parcial | Propagar el preset ya materializado a los runners confirmatorios restantes |
| 7 | Elegir mayor ELBO y guardar ELBO, iteraciones y runtime de cada inicio | El runner conserva todos los intentos y selecciona el mayor ELBO entre ajustes convergentes/finitos | P1 del paquete | Implementado | Persistir su resumen desde el runner del paper |
| 8 | Cinco inicializaciones SVI con estados comparables a VI | Un mismo `vi_initial_states` alimenta VI y SVI; semillas de inicialización y ajuste están separadas | P1 del paquete | Implementado | El runner del paper debe reutilizar el mismo objeto y fijar cinco inicios |
| 8 | Minibatch uniforme sin reemplazo y escalamiento por `m / tamaño_batch` | Está presente en los updates principales | Incluido | Implementado | Añadir pruebas de insesgamiento de estadísticas suficientes |
| 8.1 | Full ELBO cada 50; media móvil de 5; mínimo 1000; paciencia 5; máximo 50000 | Control SVI independiente y diagnósticos implementados; el preflight E3 materializa exactamente 50/5/1000/5/50000 | Segundo incremento | Parcial | Medir el `screen`, congelar la configuración y propagarla después de E3 |
| 9--15 | Sistema completo de métricas | La suite implementa recuperación alineada, cobertura, fidelidad respecto a MCMC, ARI, CPE, entropía, MLPD, regret y MCSE; E1 ya persiste métricas por tarea y agregados básicos | Métricas experimentales en esta suite | Parcial | Completar agregación por parámetro/réplica y tablas confirmatorias con escalas explícitas |
| 16--17 | Tablas 2--4 y figuras 1--5 | No existen como productos Monte Carlo | Producción automática genérica | Ausente | Generarlas desde resultados por réplica, con MCSE |
| 18 | E2 de recuperación de K | `run-e2-smoke.R` prueba CAVI/SVI sobre K candidatos y MCMC sobredimensionado, con selección por ELBO/MLPD y ocupación | Protocolo E2 versionado | Parcial | Ejecutar piloto y luego 200 datasets con la grilla confirmatoria |
| 19 | E3 factorial 4x3x3 con 20 semillas | El protocolo declara el factorial exacto; el preflight materializa 10 tareas `screen` y 109 `factorial`, con controles/presupuesto y ejecución todavía bloqueada | Protocolo E3 versionado | Parcial | Ejecutar `screen`, medir RSS/tiempo, congelar y luego ejecutar las 720 corridas SVI |
| 20 | E4 con datasets anidados y 20 réplicas | `run-e4-smoke.R` usa prefijos anidados de un dataset máximo y compara CAVI/SVI, con MCMC en el tamaño permitido | Protocolo E4 versionado | Parcial | Ejecutar tamaños y 20 réplicas confirmatorias después del piloto de recursos |
| 21 | E5 Farms completo | `run-e5-smoke.R` ajusta MCMC/CAVI/SVI, K candidatos, CV dentro de finca, WAIC marginal, estabilidad multistart y MCMC sobredimensionado | Protocolo E5 versionado | Parcial | Congelar preprocessing/controles y ejecutar el presupuesto confirmatorio |
| 22 | WAIC marginal por finca con 5000 draws | Criterio marginal integrado, metadatos y S configurable | Tercer incremento | Parcial | Configurar 5000 draws y registrar diagnósticos en runner final |
| 23 | DIC marginal sólo si se conserva | Retirado de API, helpers, resultados, tests y reporting | P3.1 cerrado; `predictive-criteria-2` contiene sólo criterios conservados | No aplica | Ninguna; el PDF lo hacía opcional |
| 24 | PPC con cinco discrepancias, 5000 datasets y MCSE | PPC streaming T1--T5 común; `posterior_predict()` y `pp_check_from_prediction()` permiten conservar/reutilizar `y_rep`; reporta MCSE IID y corregido por ESS MCMC | Tercer y quinto incrementos | Parcial | Fijar 5000 replicados en el runner y construir tabla 8 |
| 25 | Sensibilidad prior P1/P2 | Los runners usan la prior fija del review, independiente de los datos; las variantes P1/P2 aún no están declaradas como tareas | El plan solicita sensibilidad | Ausente | Declarar y ejecutar baseline, separación débil y Farms bajo ambas variantes sin empirical Bayes |
| 26 | Cuatro figuras Farms | Existen predictivas y PSM separadas | Producción automática prevista | Parcial | Construir paneles comparables y figura de intervalos post-relabeling |
| 27 | Benchmark controlado y metadatos de hardware/software | Se mide `proc.time()`; no hay manifiesto de entorno completo ni 20 repeticiones | Benchmark reproducible genérico | Parcial | Capturar entorno, controlar threads/compilación y separar tiempos |
| 28--30 | Pipeline final, orden y quince preguntas | Hay runners smoke reanudables E1--E5 y manifiestos locales; todavía no existe DAG confirmatorio ni reporte editorial | Plan general por fases | Parcial | Implementar scheduler global, auditoría de procedencia y reporte final |

## Auditoría detallada por bloque

### Modelo generador y prior

#### Lo reutilizable

- `simulate_data()` acepta DGP condicional, jerárquico o prior-predictivo y
  devuelve nombres semánticos para parámetros, componentes y asignaciones.
- Todos los modelos devuelven `simulation-result-2` con `model`, `train`,
  `truth`, `configuration`, `conditioning`, `parameter_source` y
  `design_source`; `test` sólo aparece cuando se solicita.
- `validate_grouped_data()` comprueba dimensiones y finitud.

El paper usa principalmente el DGP condicional parametrizable. La construcción
exacta de coeficientes, la prior y los valores S1--S7c son presets de esta suite,
no funciones que la API deba conocer por nombre.

#### Brechas vigentes

El DGP genérico admite covarianza arbitraria, coeficientes explícitos,
asignaciones fijas o muestreadas y train/test condicionalmente independientes.
Los simuladores históricos y los presets del protocolo fueron retirados del
paquete. La geometría y la prior exactas del review viven en
`paper/evaluation/R/scenarios.R`. El catálogo S1--S7c, semillas deterministas,
identificadores de tarea y persistencia atómica ya están materializados. Falta
una huella de configuración más fuerte que el identificador legible y la
orquestación `pilot/full` con presupuestos congelados.

- El DGP jerárquico no reemplaza el DGP condicional del protocolo:
  responde una pregunta distinta sobre variabilidad entre datasets.
- El prior de inferencia se construye por separado de `truth`. Aunque simulador
  y ajuste usan la misma jerarquía y parametrización, copiar los valores
  verdaderos a la prior contaminaría la evaluación.
- La separación euclídea de coeficientes no garantiza separación predictiva;
  debe evaluarse junto con la distribución de `X` y el ruido.
- `hyperparams_chlrm(..., "NON_INFORMATIVE")` no es independiente de los datos:
  calcula `s2_prior` desde residuos OLS y lo utiliza en `Lambda0`, `S0` y `b0`.
- Las constantes tampoco coinciden con el PDF (`nu0`, `Lambda0`, `S0` y
  `alpha0`).

#### Papel de cada DGP en la revisión

El diseño del review es un DGP condicional legítimo: fija componentes para
controlar separación, tamaño, balance y ruido, y vuelve a generar datos. No es
SBC ni una evaluación integrada sobre la jerarquía. Se implementará como un
adaptador pequeño que construye `conditional_dgp()`; no tendrá un simulador
alternativo.

En esos escenarios, $\beta_k$, $\sigma_k^2$ y $\omega$ sí tienen verdad
generadora. $\beta$, $\Sigma$ y $\xi^2$ no adquieren una verdad poblacional sólo
por resumir los componentes fijados. Si el paper reporta recuperación de esos
hiperparámetros deberá usar `hierarchical_dgp()` o declarar explícitamente un
estimando de población finita, distinto del parámetro generador jerárquico.

Para no sobregeneralizar sus conclusiones, la evaluación completa añadirá dos
validaciones pequeñas y separadas:

1. `hierarchical_dgp()` para HLRM/CHLRM, comprobando recuperación poblacional y
   predicción de grupos nuevos al volver a muestrear efectos/componentes;
2. `prior_predictive_dgp()` para LRM/HLRM/CHLRM como base de SBC de MCMC/CAVI/SVI,
   usando estadísticas invariantes a etiquetas en CHLRM.

Estas validaciones no sustituyen S1--S7c ni cambian los escenarios solicitados;
demuestran aspectos diferentes del método usando la misma API y kernels.

#### Composición objetivo entre capas

```r
scenario <- scenario_e1(id = "S1", replicate = 1L)
data <- do.call(package_simulator, scenario$generator_arguments)
prior <- do.call(package_prior_constructor, scenario$prior_arguments)
```

Los nombres definitivos de los dos constructores se fijarán en P4 y se
documentarán en `API_CONTRACT.md`. La salida debe contener `train`, `test`,
`truth` y `configuration`; el runner añade `scenario_id`, semilla de réplica y
una huella estable de configuración.

### MCMC y label switching

#### Lo reutilizable

- El Gibbs sampler guarda parámetros de componentes, pesos, asignaciones,
  parámetros globales, log-likelihood y log joint.
- `dahl_partition_chlrm()` produce una partición invariante a etiquetas.
- `fit_label_reference_chlrm()` aprende la referencia iterativa y
  `align_to_reference_chlrm()` permuta conjuntamente `Betak`, `sigma2k`,
  `omega` y `z`; `Betaj`/`sigma2j` son cantidades inducidas invariantes.
- Hay pruebas de que la permutación conjunta conserva el coclustering.

#### Brechas

- Las clases base representan una sola cadena; `fit_model()` las orquesta
  secuencialmente mediante una primitiva interna y registra formalmente
  `diagnosed/not_converged`.
- El runner exige `thin=1`, pero los scripts exploratorios antiguos aún aplican
  thinning y no deben alimentar el paper.
- La extensión por bloques y el checkpoint tipado del ajuste completo ya
  existen; el artefacto conserva internamente el estado de cada cadena. Falta
  que el scheduler del paper persista un checkpoint por tarea.
- La referencia común se aprende sobre todas las cadenas y puede pasarse al
  multistart VI/SVI; falta persistir esa relación en el manifiesto del paper.

#### Acción

Configurar el runner del paper con warm-up 10000, bloques 5000, máximo 100000 y
una ruta versionada de checkpoint por tarea. El kernel Gibbs no debe
duplicarse.

### VI y SVI

#### Lo reutilizable

- Los updates CAVI/SVI comparten estadísticas suficientes precalculadas.
- SVI muestrea grupos sin reemplazo y escala contribuciones por `m/batch`.
- Las normales se mezclan en parámetros naturales en varios updates globales.
- Se registra ELBO y tiempo acumulado.

#### Brechas

- La orquestación interna llamada por `fit_model()` permite materializar
  estados, ejecutar cualquier número de inicios, conservar todos los intentos
  y seleccionar por ELBO. Falta fijar cinco inicios y persistir semillas y
  diagnósticos en el manifiesto del paper.
- CAVI y SVI ya tienen políticas separadas con paciencia; falta materializar
  en la configuración confirmatoria sus máximos y calendarios exactos.
- VI/SVI no reanuda un loop interrumpido a mitad de iteración. Debe decidirse,
  a partir del piloto, si hace falta checkpoint interno o si la unidad de
  reanudación segura será la tarea/inicio completo.
- El análisis de tesis aún varía `rho_initial` y `kappa` en Farms, pero la suite
  del paper ya declara por separado la grilla requerida de `|S|`, `chi` y
  `tau`, y su smoke verifica una celda contra la referencia CAVI.

#### Acción

Conectar al runner las cuatro responsabilidades ya separadas en el paquete:
estado inicial serializable, ejecución de un inicio, política de convergencia y
selección multistart. El runner debe persistir su configuración y resultados,
no duplicar estas implementaciones.

### Métricas de recuperación e incertidumbre

#### Disponible

- Medias, desviaciones estándar e intervalos marginales mediante
  `tidy_parameters()`;
- alineamiento húngaro de componentes con la verdad antes de recuperación;
- bias/RMSE/cobertura/ancho por bloque, fidelidad respecto a draws MCMC y
  errores de correlación;
- ARI, CPE, PSM, probabilidades de asignación y entropía normalizada;
- MLPD, RMSE, MAE, R², cobertura y ancho predictivo fuera de muestra;
- predictive regret, recuperación de K y MCSE básicos;
- persistencia pointwise y agregación básica de E1.

#### Pendiente

- integrar razones de desviación estándar CAVI/MCMC y SVI/MCMC en los runners;
- añadir Wasserstein marginal suplementaria si el costo piloto lo permite;
- agregar todos los parámetros alineados sobre réplicas con niveles
  `per_parameter`, `per_replication` y `aggregate`;
- elevar draws y réplicas a los tamaños confirmatorios sólo después del piloto.

Las cantidades que comparan componentes se calculan desde draws o resúmenes
alineados y nunca desde etiquetas crudas. Los agregados no reemplazarán la
evidencia pointwise necesaria para recalcular MCSE y coberturas.

### WAIC, validación predictiva y PPC

El bloqueo detectado era semántico, no cosmético. La función retirada
`waic_chlrm(unit="group")` seleccionaba `Betaj/sigma2j` según una asignación.
La ruta `predictive_criteria()` aplica la fórmula marginal requerida por el PDF:

$$
\ell_j^{(s)}=
\operatorname{logsumexp}_k\left[
\log\omega_k^{(s)}+
\sum_i\log N(y_{ij}\mid x_{ij}^{\top}\beta_k^{(s)},\sigma_k^{2(s)})
\right].
$$

Por tanto:

- se eliminó el WAIC condicional CHLRM; no se conserva alias `*_legacy`;
- `loglik_group_marginal_chlrm()` ya implementa esta fórmula y tiene pruebas;
- `vi_sample_draws_chlrm()` ya incluye draws Dirichlet de `omega`; para esta
  likelihood se debe usar `include_allocations = FALSE`;
- MCMC, VI y SVI ya pasan la misma definición de matriz `m x S` a
  `waic_from_loglik()`;
- DIC fue retirado: su punto en la media posterior depende del relabeling y no
  representa de forma estable una posterior multimodal/singular;
- WAIC se comparará sólo con igual dataset, escala, unidad predictiva y
  marginalización; PSIS-LOO con diagnóstico o CV/LOGO y MLPD serán la evidencia
  predictiva principal cuando sea viable.

El PPC estándar reproduce T1--T5 cuando son aplicables: media global,
desviación estándar, varianza de medias grupales, promedio de varianzas
intragrupo y MSE respecto al OLS fijado con los datos observados. Reporta tanto
MCSE IID nominal como MCSE corregido con ESS de la secuencia indicadora por
cadena MCMC. Los scripts exploratorios siguen usando menos de 5000 draws. El
PPC guarda sólo resúmenes por defecto. La API pública
`posterior_predict()` puede conservar replicados agrupados y
`pp_check_from_prediction()` los reutiliza sin ajustar ni simular de nuevo. No
se fuerza al PPC a materializarlos: ambas rutas comparten la misma primitiva de
respuesta normal y tienen responsabilidades separadas.

### Estado de los cinco experimentos

#### E1: estudio Monte Carlo principal

Estado: infraestructura `smoke` implementada. El catálogo declara S1--S7c y el runner E1
`smoke` instala el paquete, ejecuta MCMC/CAVI/SVI, persiste cada tarea de forma
atómica y reanuda sin repetirla. El schema vigente incluye ARI, CPE, entropía y
desempeño predictivo realmente retenido, además de recuperación paramétrica
alineada con verdad. Faltan todos los escenarios en `pilot/full`, la agregación
paramétrica completa y la corrida confirmatoria de 200 réplicas.

#### E2: identificación y selección de K

Estado: infraestructura `smoke` implementada. La suite genera las dos
separaciones, ajusta una grilla reducida CAVI/SVI, ejecuta MCMC
sobredimensionado, calcula MLPD retenido, ARI/CPE/entropía y produce selección
por ELBO, MLPD o moda de componentes ocupados. El smoke usa datasets simulados
y un test independiente; no debe confundirse con la selección exploratoria de
K en Farms. Falta expandirlo a las 200 réplicas y la grilla confirmatoria.

#### E3: sensibilidad de SVI

Estado: infraestructura `smoke` y preflight implementadas. El protocolo declara
el factorial exacto; `run-e3-smoke.R` ajusta una referencia CAVI y una celda
SVI, mientras `run-e3-preflight.R` materializa un `screen` de 10 tareas y el
piloto factorial de 109 sin hacer ajustes. La referencia CAVI usa cinco inicios;
cada semilla SVI es deliberadamente una corrida independiente con un inicio,
para conservar las 720 corridas especificadas y medir fallos de convergencia.
El análisis de tuning de la tesis permanece separado y no alimenta E3. Falta
ejecutar/medir el `screen`, congelar y luego correr las 720 SVI confirmatorias.

#### E4: escalabilidad

Estado: infraestructura `smoke` implementada. `run-e4-smoke.R` simula una sola
población máxima y usa prefijos anidados para comparar CAVI/SVI; MCMC se limita
al tamaño pequeño previsto. Persiste tiempos, memoria estimada, MLPD, ARI y
pendientes log--log. Faltan los tamaños exactos, 20 réplicas pareadas y
mediciones empíricas de RSS.

#### E5: Farms

Estado: infraestructura `smoke` implementada. `run-e5-smoke.R` usa Farms con
escalado reversible, ajusta MCMC/CAVI/SVI para una grilla mínima de K, calcula
WAIC marginal por finca y CV de cinco folds dentro de finca con MLPD/MSE, mide
estabilidad multistart y ejecuta MCMC sobredimensionado para ocupación. Falta
la grilla K confirmatoria con su regla de extensión, cuatro cadenas/veinte
inicios, la distribución `p(K*|y)` si el protocolo la mantiene y los paneles
Farms finales.

## Brechas entre paquete y protocolo

`PACKAGE_ROADMAP.md` cubre la coherencia bayesiana, contratos, relabeling,
WAIC/LOO/CV marginales, intervalos, organización y performance. Para
responder completamente al PDF, la suite del paper debe declarar o
materializar:

1. extender el runner reanudable E1, ya probado en `smoke`, a `pilot/full`
   después de congelar la configuración SVI y los límites de recursos;
2. el uso explícito de la prior del protocolo ya implementada;
3. las configuraciones exactas de parada MCMC, CAVI y SVI sobre los controles ya
   disponibles;
4. la persistencia del multistart y de la referencia común ya implementados;
5. definiciones de datos crudos y métricas para cada réplica;
6. los tamaños R=200, RK=100, RS=20 y las grillas exactas;
7. configurar para E5 los planes ya disponibles de CV dentro de finca y LOGO
   de finca nueva, sin reimplementar sus particiones ni su scoring;
8. las tablas y figuras enumeradas como artefactos con fuentes de datos;
9. captura del entorno computacional y control de threads;
10. rutas por tarea para checkpoints MCMC, reanudación del pipeline, ejecución
    paralela controlada y presupuestos globales de memoria;
11. perfiles de ejecución `smoke`, `pilot` y `full`;
12. criterios de aceptación numéricos para cada experimento;
13. una decisión explícita sobre escalado; la suite entrega
    `data_scaling()` y los kernels nunca transforman silenciosamente.

### Política de estandarización

No existe ya un conflicto de arquitectura: los kernels no escalarán
silenciosamente, pero el paquete ofrecerá `data_scaling()` y un mapa reversible.
El paper lo configurará para Farms y guardará uno distinto por conjunto de
entrenamiento/fold. E1--E4 no se reescalan de forma
dependiente de los datos: sus escalas, $\Sigma_x$, $d$, coeficientes y
varianzas forman parte del DGP y transformarlas alteraría el estimando
declarado.

Para E5, cada conjunto de entrenamiento —y cada fold— aprende
$\mu_x,s_x,\mu_y,s_y$ y aplica

$$
x_r^*=(x_r-\mu_{x,r})/s_{x,r},\qquad
y^*=(y-\mu_y)/s_y.
$$

La API del paquete ya aprende y aplica el mapa dentro de cada fold mediante
`cross_validate_model(..., scaling = data_scaling(...))`. La suite sólo declara
las columnas y el esquema; no crea `preprocessing.R` ni duplica las fórmulas.
Debe verificar en sus artefactos que:

- los predictores seleccionados con varianza nula sean rechazados antes del
  ajuste;
- cada mapa use únicamente estadísticas de entrenamiento y quede disponible
  mediante `scaling_map()` cuando se retengan ajustes piloto;
- los resultados públicos satisfagan $y=\mu_y+s_y y^*$;
- las pendientes satisfagan
  $\beta_r=(s_y/s_{x,r})\beta_r^*$ e interceptos por
  $\beta_0=\mu_y+s_y\beta_0^*-\sum_r\beta_r\mu_{x,r}$;
- las varianzas residuales satisfagan $\sigma^2=s_y^2\sigma^{2*}$ y las
  covarianzas usen la transformación lineal completa; y
- coeficientes poblacionales, de grupo y de componente vuelvan a escala de
  entrada sin cambiar pesos ni asignaciones.

Los criterios predictivos se calculan en una sola escala declarada por tabla.
No se compararán log densities obtenidas bajo transformaciones distintas sin
incluir el jacobiano correspondiente.

## Escala computacional y ejecución segura

El protocolo completo no debe lanzarse como un único script en un portátil.
Sólo E1 requiere aproximadamente

$$
9\times200\times(4+5+5)=25200
$$

ajustes si se cuenta cada cadena e inicialización. E2 añade aproximadamente
32000 ajustes variacionales, además de MCMC. E3 añade 720 corridas SVI. E5
añade al menos 528 ajustes completos antes de contar la validación cruzada.

Se requieren tres perfiles:

| Perfil | Propósito | Ejemplo |
|---|---|---|
| `smoke` | validar contratos y archivos | 1 escenario, 2 réplicas, cadenas cortas, 2 inicios |
| `pilot` | estimar tiempo/memoria y depurar métricas | todos los escenarios, 5--10 réplicas |
| `full` | resultados publicables | tamaños exactos del PDF |

Cada tarea mínima debe ser una combinación independiente
`experimento/escenario/replica/metodo/inicio-o-cadena`, escribir primero a un
archivo temporal y promoverse atómicamente al finalizar. El agregador nunca
debe volver a ajustar modelos.

## Arquitectura vigente de la suite del paper

```text
paper/evaluation/
├── R/
│   ├── scenarios.R
│   ├── tasks.R
│   ├── e2.R ... e5.R
│   ├── metrics.R
│   ├── sbc.R
│   ├── provenance.R
│   ├── reporting.R
│   ├── aggregation.R
├── run-smoke.R
├── run-e2-smoke.R ... run-e5-smoke.R
├── run-smoke-figures.R
├── run-validation.R
├── tests/testthat/
└── output/smoke/
    ├── tasks/                 # RDS E1
    └── E2/ ... E5/           # RDS y CSV por experimento
```

La lógica estadística reusable —priors, estados iniciales, multichain,
multistart, relabeling, draws, scores y diagnósticos— debe permanecer en el
paquete `R/`. La suite sólo define escenarios/configuraciones, expande y ejecuta
tareas, calcula métricas contra verdad o referencia y produce artefactos. Antes
del piloto se añadirá una configuración serializable compartida y un scheduler
global; no se creará otro árbol `experiments/` que separe el código de sus
salidas sin necesidad.

## Esquema de resultados requerido

Como mínimo:

- `task_manifest.csv`: configuración, semilla, commit, estado y rutas;
- `provenance.csv`/`provenance.rds`: huellas de configuración, código y run,
  versión del paquete, commit, R, plataforma y threads;
- `fit_diagnostics.csv`: convergencia, iteraciones, runtime, R-hat y ESS;
- `posterior_summary.csv`: método, réplica, parámetro, verdad, media, sd e
  intervalo;
- `clustering_metrics.csv`: ARI, CPE, entropía y componentes ocupados;
- `predictive_pointwise.csv`: observación/grupo, log score y fold;
- `model_criteria.csv`: ELBO, WAIC marginal, PSIS-LOO/diagnósticos cuando sea
  válido y resultados de CV/LOGO;
- `environment.txt`: R, paquetes, compilador, BLAS, OpenMP, CPU, RAM y threads.
- `figures/figure-manifest.csv`: ruta, perfil, huella de todos los CSV fuente y
  marca explícita de que las figuras smoke no son resultados científicos.

Los objetos de draws deben guardarse de forma selectiva o comprimida. Las
tablas agregadas no pueden ser la única evidencia porque impedirían recalcular
MCSE, coberturas, correlaciones o criterios pointwise.

## Pruebas de aceptación antes de experimentar

Las pruebas actuales cubren los 17 puntos en ejemplos pequeños mediante tests
del paquete, 32 bloques propios de la suite y ejecuciones smoke reanudables:
generador/prior, equivalencia CAVI--SVI, reglas de parada, relabeling, WAIC
marginal, predicción, reutilización de `y_rep`, métricas manuales, escalado sin
fuga, persistencia y los tres DGP. También cubren contratos SPD, RNG,
multistart, multichain y reanudación MCMC exacta. Esto valida implementación e
integración, pero no demuestra calibración inferencial ni propiedades
asintóticas; esas requieren `pilot/full` con muchas réplicas.

Los puntos 2--9, 13, 15--17 son pruebas de dependencia que deben vivir con el
paquete. Los puntos 1, 10--12 y 14 son aceptación de la suite y deben vivir
bajo `paper/evaluation/`; aquí sólo se registra conjuntamente la puerta previa
al perfil confirmatorio.

1. El generador alcanza la separación mínima solicitada para distintos K y p.
2. `Sigma_x` empírica converge a la configurada.
3. Train/test comparten asignaciones y son condicionalmente independientes.
4. La prior principal no accede a datos.
5. Con minibatch completo y tasa uno, un paso SVI coincide con CAVI.
6. Las reglas de parada requieren exactamente la paciencia especificada.
7. El relabeling iterativo es idempotente y aplica una permutación conjunta.
8. WAIC marginal es invariante a permutaciones y coincide con K=1.
9. Los adaptadores WAIC/LOO consumen la misma likelihood marginal y conservan
   unidad, escala y condicionamiento.
10. Cobertura, CPE, ARI, MLPD y MCSE coinciden con ejemplos manuales pequeños.
11. La estandarización de CV no usa información del fold de prueba.
12. Reanudar una tarea produce el mismo resultado que ejecutarla de una vez.
13. Predicción de grupos observados y nuevos comparte estados latentes dentro
    de grupo, propaga asignaciones CHLRM y reutiliza `y_rep` en PPC.
14. Estandarizar y desestandarizar parámetros/predicciones reproduce la escala
    original; las log densities incluyen el jacobiano si cambia su escala.
15. MCMC, CAVI y SVI reciben el mismo objeto prior y éste no contiene valores de
    `truth` salvo en pruebas algebraicas marcadas.
16. El DGP jerárquico reproduce los momentos poblacionales teóricos; train/test
    de las mismas unidades conserva latentes y el test de unidades nuevas los
    vuelve a muestrear.
17. El DGP prior-predictivo genera $(\theta,y)$ desde exactamente la prior y la
    likelihood usadas por el ajuste; SBC usa cantidades identificables o
    invariantes a etiquetas en CHLRM.

El smoke ejecutable de 17 vive en `paper/evaluation/R/sbc.R` y sus pruebas.
Comprueba conexión prior-predictiva--ajuste--draws y la invariancia a
permutaciones, pero no evalúa uniformidad de rangos: eso requiere muchas
réplicas y pertenece a los perfiles `pilot`/`full`.

## Plan de ejecución de la evaluación

### Puerta P: dependencias del paquete

Esta puerta se consulta aquí, pero se implementa en `PACKAGE_ROADMAP.md`:

1. `conditional_dgp()`, `hierarchical_dgp()`, `prior_predictive_dgp()` y
   `simulate_data()` exportados con `simulation-result-2`; el review consume el
   DGP condicional y las validaciones complementarias usan los otros dos;
2. draws, likelihood marginal, predicción, PPC y criterios con contratos
   estables;
3. multichain, multistart, relabeling y checkpoint explícitos;
4. scoring fuera de muestra y CV con objetivo/unidad declarados;
5. paquete instalable que aprueba `R CMD check`, estima/rechaza presupuestos
   por ajuste y contrasta esos estimados con mediciones piloto.

Los puntos 1--4 y la instalación/check/guarda preventiva del punto 5 están
probados. El smoke actual instala el paquete en una biblioteca temporal y la
suite ejecuta 24 pruebas de frontera sólo con exports. Cada tarea E1 smoke
registra tiempo de ajuste, memoria estimada y tamaño real del artefacto. La
medición empírica de RSS y costo en todos los tamaños continúa pendiente de
P5/piloto, por lo que la Puerta P completa aún no autoriza el perfil
confirmatorio.

SVI para LRM ya pertenece al paquete, pero no a esta puerta porque el PDF
evalúa SVI en CHLRM; la suite no debe reimplementarlo.

### Fase E-A: especificación ejecutable

Responsabilidad exclusiva de `paper/evaluation/`:

1. **Implementado:** declarar S1--S7c y las configuraciones E2--E5;
2. **Implementado:** `make_review_component_coefficients()` verifica la
   distancia antes de entregar parámetros a `conditional_dgp()`;
3. **Implementado:** la prior principal se construye con argumentos explícitos
   de `chlrm_prior()`; las sensibilidades quedan declaradas para el perfil final;
4. **Implementado en smoke:** Farms usa `data_scaling()` dentro de cada ajuste
   y CV aprende el mapa sólo con el fold de entrenamiento;
5. **Parcial avanzado:** CPE, ARI, entropía, MLPD, recuperación escalar,
   cobertura, fidelidad respecto a MCMC, regret y MCSE tienen funciones puras y
   ejemplos manuales; el adaptador CHLRM alinea componentes con la verdad y E1
   lo persiste y agrega por coordenada entre réplicas. Falta integrar las
   comparaciones de fidelidad CAVI/SVI contra MCMC en el runner confirmatorio;
6. **Implementado para smoke:** existen schemas crudos/agregados versionados y
   `paper-run-provenance-1` registra huellas de configuración/código, paquete,
   commit, R, plataforma y threads. El inventario detallado de hardware y RSS
   medido se completa en el piloto.

Criterio de salida: cada receta produce exactamente los argumentos, semillas y
estimandos descritos en el PDF sin ejecutar todavía miles de ajustes.

### Fase E-B: runner reanudable y perfil `smoke`

1. Expandir cada tarea como experimento/escenario/réplica/método/cadena-o-inicio.
2. Instalar y consumir el paquete; bloquear cualquier `source()` de kernels.
3. Persistir resultados atómicamente, checkpoints, errores y procedencia.
4. Ejecutar E1--E5 con tamaños mínimos y regenerar todas las tablas/figuras de
   prueba desde artefactos crudos.
5. Verificar que reanudar no repite tareas completas ni cambia resultados.

Criterio de salida: el pipeline puede detenerse y reanudarse sin intervención,
y cada celda de prueba tiene una fuente rastreable.

Estado: **parcial avanzado**. E1--E5 tienen runners smoke reanudables que
instalan el paquete y escriben RDS por tarea más CSV inspeccionables. E1 agrega
recuperación/clustering/predicción; E2 selecciona K; E3 compara SVI con CAVI;
E4 usa tamaños anidados; y E5 prueba WAIC, CV dentro de finca, estabilidad y
MCMC sobredimensionado. La agregación E1 ya conserva coordenadas alineadas,
sesgo, RMSE, cobertura, ancho y MCSE; seis figuras smoke se regeneran sólo desde
CSV. El preflight E3 ya materializa los controles, 10 tareas `screen` y 109
tareas `factorial`, y escribe presupuesto/puerta sin ejecutar ajustes. Falta
medir tiempo, RSS y almacenamiento reales con `screen`.

### Fase E-C: perfil `pilot`

1. Ejecutar todos los escenarios con 5--10 réplicas según costo.
2. Medir tiempo, RAM y almacenamiento por unidad de trabajo.
3. Ajustar paralelismo y límites globales sin cambiar estimandos.
4. Ejecutar primero el piloto E3 y congelar la configuración SVI.
5. Decidir con evidencia si VI/SVI necesita checkpoint intra-loop o si basta
   repetir únicamente el inicio incompleto.

### Fase E-D: experimentos confirmatorios

Congelar commit del paquete, configuración y semillas. Ejecutar, en orden, E3,
E1, E2, E4 y E5. E5 se ejecuta al final para evitar decisiones guiadas por los
datos reales.

### Fase E-E: auditoría y manuscrito

1. Generar tablas 1--8 y S1--S3.
2. Generar figuras 1--10 y las cuatro figuras Farms.
3. Comprobar procedencia de cada celda y panel.
4. Responder explícitamente las quince preguntas de la sección 30.
5. Actualizar el LaTeX únicamente desde artefactos del pipeline congelado.

## Prioridad inmediata de esta suite

Mientras el paquete cierra P3/P4, la evaluación puede avanzar sin duplicarlo:

1. ejecutar las 10 tareas del `screen` E3 de manera secuencial y medir
   RSS/tiempo/almacenamiento;
2. revisar fallos y congelar la huella de configuración antes de habilitar las
   109 tareas `factorial` del piloto;
3. fijar límites globales antes de E1/E2/E4/E5 y mantener `full` bloqueado.

No corresponde a esta suite implementar SVI para LRM —ya disponible en el
paquete—, likelihoods, kernels,
relabeling o clases R6. Tampoco se ejecutará `full` antes de completar la puerta
P y el piloto. Los tests unitarios del paquete y las corridas `smoke` no son
resultados del artículo.
