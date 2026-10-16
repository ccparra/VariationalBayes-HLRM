# Registro histórico de decisiones del desarrollo estadístico

> **Estado: histórico y no normativo.** Este archivo conserva decisiones,
> alternativas descartadas y estados observados durante distintos incrementos.
> Puede mencionar APIs, pendientes o estructuras que ya fueron sustituidos. No
> debe usarse para decidir qué sigue ni para inferir la API vigente.

El control actual está separado en [`PACKAGE_ROADMAP.md`](PACKAGE_ROADMAP.md) y
[`paper/EVALUATION_PLAN.md`](paper/EVALUATION_PLAN.md), coordinados por
[`ROADMAP.md`](ROADMAP.md).

## Propósito histórico

Este documento consolida la auditoría de `R/` y define el orden seguro para
corregir, validar y reorganizar LRM, HLRM y CHLRM. El objetivo no es cambiar
nombres o mover archivos inmediatamente, sino preservar primero el
comportamiento correcto, resolver ambigüedades estadísticas y verificar la
correspondencia entre teoría, código y experimentos.

El criterio de terminación es operacional: el repositorio debe poder ejecutar
el protocolo de `paper/references/Computacion_revision_paper_Camilo.pdf`
mediante una interfaz estable, reproducible y segura en recursos. La
trazabilidad requisito por requisito se mantiene en
`paper/EVALUATION_PLAN.md`. Ese documento define
**qué debe poder evaluarse**; este plan define **cómo se modificará el código
para hacerlo posible**. Una fase no se considera completa sólo porque exista
una función: debe satisfacer sus pruebas y artefactos de aceptación.

El núcleo conjugado MCMC y CAVI es razonablemente coherente. Los riesgos más
importantes están en los contratos públicos, la predicción jerárquica, las
priors por defecto, la parada de SVI, los diagnósticos, la reproducibilidad y la
estructura de dependencias mediante `source()`.

## Avance: primer incremento (24 de septiembre de 2026)

Se conserva el plan de ejecución completo; este incremento no da por cerrada
ninguna fase transversal. Se adelantaron contratos mínimos del protocolo para
tener ejemplos verificables antes de mover físicamente las clases.

| Entregable | Implementación | Evidencia |
|---|---|---|
| Prior independiente y separación geométrica | `R/chlrm-protocol.R` | valores exactos, dimensiones y distancia mínima |
| Simulación fija con covariables correlacionadas y train/test | `simulate_chlrm_protocol()` | reproducibilidad, RNG preservado, tamaños, covarianza y ruido |
| Inicialización sin dependencia OLS de la prior | `Initialize_chlrm_vi()` y opción de ejecución `initialization` | ajustes cortos VI/SVI con la nueva prior |
| Likelihood marginal por finca | `R/chlrm-likelihood.R` | ejemplo manual, K=1, ceros en pesos, estabilidad e invariancia a etiquetas |
| Draws compactos de componentes y pesos | `vi_sample_draws_chlrm(include_allocations = FALSE)` | momentos Dirichlet y compatibilidad con likelihood |
| Conexión con WAIC genérico | `waic_from_loglik(ll)` | cálculo manual y permutación conjunta por draw |

Validación: `Rscript tests/testthat.R` pasa con los tests existentes y
`tests/testthat/test-chlrm-protocol.R`. Los ajustes son de integración, no una
evaluación de convergencia, cobertura, rapidez comparativa o fidelidad posterior.
No se ejecutó E1--E5 ni se certifica todavía un paquete mediante `R CMD check`.

Limitaciones de aquel incremento: las clases mantuvieron sus loaders `source()`;
GOF/DIC aún no cambiaban (sustituidos en el tercer incremento); no había catálogo de escenarios, multistart,
checkpoint ni cota global de RAM. La nueva likelihood evita tensores
observación × componente × draw, pero sigue recorriendo datos: optimizarla
requiere primero medir tiempo/memoria y conservar la estabilidad numérica.

Los incrementos segundo y tercero siguientes detallan el control de ejecución
y la integración marginal en GOF. Las fases transversales siguen pendientes.

## Avance: segundo incremento, control de ejecución CHLRM

Implementación en `R/inference-control.R`, conectada a las clases CHLRM:

- `run_settings$seed` controla constructor y ajuste, restaura el RNG del
  llamador y conserva el estado del generador entre ambas llamadas. Sin semilla
  se conserva el comportamiento con RNG global; `fit(seed=...)` sólo reemplaza
  la secuencia del ajuste, no la inicialización.
- Validación temprana de ajustes VI/SVI, tamaños de minibatch y tasa.
- CAVI: cambio relativo con denominador `1 + abs(ELBO_anterior)`, tolerancia
  predeterminada `1e-8` y cinco iteraciones consecutivas.
- SVI: ELBO completo cada 50 iteraciones, media móvil de cinco checkpoints,
  mínimo 1000 iteraciones, tolerancia `1e-4` y cinco checkpoints elegibles.
  Los controles son configurables para tests; no se ejecutó el perfil completo.
- Diagnósticos explícitos de tolerancia alcanzada o máximo agotado; error ante
  ELBO no finito. Estado `initialized/running/fitted/failed`, comprobaciones
  antes de inferencias/GOF y rechazo de un segundo `fit()` sobre el mismo objeto.
- Prueba de serialización de un objeto inicializado y de equivalencia del
  paso SVI con CAVI para batch completo y tasa uno.

Evidencia: `tests/testthat/test-inference-control.R`; pasa la suite completa
con `Rscript tests/testthat.R`. La convergencia aquí sólo significa estabilidad
del ELBO bajo esa regla: no garantiza óptimo global ni calibración posterior.

**No cerrado:** el contrato de estado completo y reanudable, bloques MCMC,
multistart, streams paralelos, migración RNG a LRM/HLRM y utilidades de draws,
presupuesto global de recursos, reorganización física y empaquetado. Las
funciones de control son reutilizables, pero sólo CHLRM las usa por ahora.
En ese incremento GOF aún era condicional y dependía de construir primero el
resumen posterior; ambos puntos quedaron sustituidos en el incremento siguiente.

## Avance: tercer incremento, evaluación marginal y retirada de rutas obsoletas

- `R/chlrm-evaluation.R` centraliza GOF para CHLRM MCMC/VI/SVI. Las clases no
  duplican WAIC, DIC, PPC ni infieren ejes de draws por coincidencias de tamaños.
- La evaluación funciona directamente después de `fit()`, sin depender de
  `summarize_posterior()`. La API monolítica usada en este incremento fue
  sustituida posteriormente por tres responsabilidades independientes.
- DIC reutiliza la likelihood marginal: media MCMC tras alineamiento húngaro
  (iterativo desde el sexto incremento); media analítica de q para VI/SVI. Si `ak <= 1`, la media de
  varianza no existe y DIC se declara no disponible, sin inventar un promedio.
- Se eliminaron `waic_chlrm`, `dic_chlrm`, `ppval_chlrm`, `yrep_chlrm` y
  `logmeanexp_safe`. No se conservan wrappers ni nombres `*_legacy` en R/.
  La likelihood condicional del Gibbs sigue siendo necesaria para su log joint:
  no es código obsoleto ni se usa como WAIC marginal.
- PPC procesa un dataset a la vez; sólo conserva todos los replicados si se
  solicita expresamente y cabe su matriz en el presupuesto. Draws VI/SVI para
  GOF son compactos, sin duplicación por finca ni muestreo de z para WAIC.
- Los consumidores de `thesis/final/code` usan el nuevo contrato y columnas
  `WAIC_group_marginal`. LOGO reutiliza la misma likelihood en vez de duplicarla.
  El reporte compartido separa los campos condicionales de HLRM con NA; eso no
  es una ruta de compatibilidad con el GOF CHLRM retirado.

Evidencia: `tests/testthat/test-chlrm-evaluation.R`, además de las suites previas.
Pruebas de fórmulas manuales, permutaciones, K=1/p=1, presupuesto, reproducibilidad,
streaming PPC y evaluaciones sin inferencias previas. Se corrigieron además
colapsos de dimensiones para `p=1`, `K=1` y minibatch de un grupo en las
estadísticas VI, responsabilidades y parámetros inducidos de MCMC.
No se regeneraron outputs
experimentales ni se modificó `thesis/v0`; las salidas históricas no cambian de
significado retroactivamente.

Pendiente: alineamiento iterativo/común entre métodos, PPC T1--T5 y MCSE con
dependencia MCMC, CV/PSIS, múltiples cadenas/inicios, migración de LRM/HLRM y
reorganización completa. La existencia de esta ruta no certifica resultados
publicables ni completa por sí sola la fase de evaluación.

## Avance: cuarto incremento, evaluación LRM/HLRM conectada y validada

Se cierra el bloque de evaluación ejecutable entre familias, no el plan entero:

- LRM MCMC/VI y HLRM MCMC/VI/SVI evalúan después de `fit()` sin inferencias
  previas; RNG local durante GOF, presupuestos explícitos y medias analíticas
  variacionales para DIC/MSE. LRM rechaza SVI porque no está implementado.
- `R/gaussian-evaluation.R` separa criterios observacionales y grupales
  condicionales. `R/evaluation-utils.R` concentra WAIC y PPC streaming para
  todos los modelos, incluido CHLRM; se eliminan las ocho funciones GOF/PPC
  sustituidas de LRM/HLRM, sin mantener wrappers muertos.
- `R/hlrm-marginal.R` integra beta_j analíticamente y sigma2_j mediante
  cuadratura Gauss-Laguerre adaptada a la precisión Gamma. No forma matrices
  n_j × n_j; usa una descomposición de rango a lo sumo p. Tiene límite de
  nodos/trabajo, diagnóstico de acuerdo entre reglas y error si no converge.
- HLRM permite el objetivo `new_groups` con draws globales de MCMC o q e
  incluye WAIC marginal. DIC fue retirado posteriormente de toda la API.
- LRM admite `groups` para sumar likelihood por las mismas fincas. Los tres
  exponen **el mismo nombre `waic_group_marginal`**, sin alias `waic_marginal`.
  La comparación requiere igual conjunto de datos, escala y agrupación; no
  implica que q represente exactamente la posterior ni que MCMC haya convergido.
- Reportes y figuras activos migrados a los nuevos contratos. Un reporte
  rechaza esquemas antiguos y conserva los campos condicional/marginal separados.
- La posterior predictiva de observaciones de grupos conocidos está disponible
  en las tres familias: MCMC usa draws posteriores y VI/SVI draws de q. Se pueden
  pedir replicados completos para densidades/intervalos dentro del presupuesto;
  los PPC guardan sólo resúmenes por defecto. `predict()` no sustituye esa distribución.

Aceptación: pruebas manuales de WAIC, acuerdo de cuadratura con el límite
analítico Sigma=0 y con integración adaptativa independiente para p=1 y p=2,
momentos de q global, reproducibilidad, ausencia de mutación parcial al fallar,
dimensiones unitarias, límites de memoria/trabajo y reportes. Suites nuevas:
`test-gaussian-evaluation.R` y `test-hlrm-marginal.R`; también pasan las anteriores.
La cuadratura comprueba acuerdo numérico, no ofrece una cota matemática rigurosa
del error. No se ejecutan estudios grandes ni se reescriben outputs históricos.

La marginalización HLRM dejó de ser un hueco de este bloque. El resto del plan
sigue requiriendo módulos completos: PPC T1--T5/MCSE dependiente, validación
cruzada/PSIS, alineamiento iterativo, cadenas/inicios múltiples, reanudación,
uniformidad de los contratos de ajuste LRM/HLRM y empaquetado. Cada bloque debe
cerrarse con consumidores y pruebas, no dejar implementaciones provisionales
sin conectar ni considerarse terminado por tener sólo una función aislada.

## Avance: quinto incremento, interfaz predictiva posterior uniforme

Se cerró la separación entre media latente, respuesta futura y PPC para las
nueve combinaciones implementadas: LRM, HLRM y CHLRM con MCMC/CAVI/SVI. SVI
usa observaciones como unidad de minibatch en LRM y grupos en HLRM/CHLRM.

- `R/posterior-prediction.R` define un único esquema
  `posterior-prediction-2`, con un núcleo común y campos específicos sólo
  cuando corresponden.
  `predict()` resume draws de
  $E(\widetilde y\mid\theta,X)$; `posterior_predict()` genera además el ruido
  observacional y resume/puede conservar $\widetilde y$.
- Ninguna de las dos llamadas depende de `summarize_posterior()`. Ambas usan draws
  conjuntos de la cadena o de q, semilla aislada y presupuesto explícito.
- En HLRM, `group_ids` diferencia grupos observados de nuevos. Para cada grupo
  nuevo y draw se genera una sola pareja $(\beta_j,\sigma_j^2)$ desde la
  jerarquía global y se comparte entre todas sus filas.
- En CHLRM, un grupo observado usa su draw posterior de etiqueta MCMC o su
  responsabilidad variacional; un grupo nuevo genera una etiqueta desde
  $\omega^{(s)}$. La etiqueta se comparte por todo el grupo. La predicción ya
  no clasifica rígidamente con `z_pred`, por lo que propaga incertidumbre de
  asignación y permanece invariante a permutaciones conjuntas por draw.
- Los intervalos devueltos por `predict()` son creíbles puntuales de la media;
  los de `posterior_predict()` son predictivos puntuales de la respuesta. No
  se presentan como bandas simultáneas ni como intervalos frecuentistas.
- `pp_check_from_prediction()` reutiliza draws almacenados para PPC y gráficos.
  El GOF automático comparte la misma primitiva de generación gaussiana, pero
  conserva su ruta streaming para no materializar todos los replicados.
- Se retiraron los seis helpers `yhat_*` sustituidos, sin wrappers legacy, y se
  migraron todos los consumidores activos de `thesis/final/code` al esquema
  nuevo. `thesis/v0` permanece histórico e independiente.

Evidencia: `tests/testthat/test-posterior-prediction.R` y ampliaciones de las
suites gaussianas/CHLRM prueban reproducibilidad, dimensiones, presupuesto,
grupos observados/nuevos, reutilización PPC e invariancia a etiquetas. La suite
completa pasa. No se ejecutaron experimentos grandes ni se regeneraron outputs.

Este incremento cierra la API predictiva de la fase 1 y los intervalos básicos
de la fase 3. No cierra la fase 3 completa: faltan MLPD/test unificado,
LOGO/PSIS-LOO, discrepancias T1--T5, MCSE dependiente para MCMC y bandas
simultáneas si el manuscrito finalmente las exige.

## Avance: sexto incremento, separación de productos y relabeling común

El trabajo se dividió en un roadmap del paquete y un protocolo de evaluación
del paper, coordinados por `ROADMAP.md`. Las referencias editoriales activas se
movieron a `paper/references/`; la tesis permanece como producto separado.

En el paquete, `R/chlrm-relabeling.R` implementa el algoritmo del PDF: ubicación
$\eta_k=(\beta_k^\top,\log\sigma_k^2)^\top$, referencia inicial del draw con
mayor log joint, desviaciones combinadas de todos los draws, asignación húngara,
actualización de la referencia y repetición hasta que las permutaciones no
cambian. `align_to_reference_chlrm()` aplica el mismo objeto a draws MCMC o a
resultados VI/SVI; para estos últimos usa $E_q[\log\sigma_k^2]$, que existe aun
cuando la media aritmética de la varianza no existe.

Se retiró la función de una pasada. Las clases MCMC y VI/SVI exponen y guardan
la referencia común; DIC MCMC usa el nuevo alineamiento. Las pruebas verifican
escala combinada, máximo log joint, idempotencia, permutación conjunta,
preservación del coclustering y alineación MCMC--VI. Falta el orquestador
multichain que proporcione draws combinados de cadenas ya diagnosticadas.

### Estado consolidado del plan tras el sexto incremento

| Fase | Ya disponible y probado | Falta para cerrarla | Estado |
|---|---|---|---|
| 0. Caracterización y notación | tests algebraicos, dimensiones degeneradas, simuladores y esquemas GOF/predicción | mapa único de notación; colinealidad y catálogo completo de bordes | Parcial |
| 1. Seguridad y contratos | validación agrupada, presupuesto local, RNG/persistencia y parada CHLRM, API predictiva uniforme | estado/RNG/parada uniforme de LRM/HLRM, validadores completos, reanudación y `logdet` SPD | Parcial |
| 2. Coherencia bayesiana | prior del protocolo, separación prior/inicialización, equivalencia SVI-CAVI y relabeling iterativo/común | auditoría ELBO/condicionales completa y sensibilidad de priors | Parcial |
| 3. Evaluación | WAIC con unidad explícita, HLRM marginal, PPC streaming, `predict()`/`posterior_predict()` e intervalos puntuales | métricas Monte Carlo, T1--T5, MCSE MCMC, MLPD, LOGO/CV/PSIS y estudios de cobertura | Parcial avanzado en infraestructura |
| 4. Arquitectura de paquete | separación inicial de núcleo, tests, tesis v0/final y consumidores activos migrados | `DESCRIPTION`/`NAMESPACE`, retirar `source()` internos, nombres/exportaciones y configuración declarativa | Parcial inicial |
| 5. Experimentos | scripts exploratorios y diseños parciales | runner seguro, smoke/pilot, E1--E5 completos, agregación y artefactos finales | No iniciada como fase confirmatoria |

La infraestructura estadística básica ya permite ajustar, puntuar y predecir,
pero todavía no permite afirmar que el método haya sido evaluado según el PDF.
La distancia restante está principalmente en comparación rigurosa de métodos,
orquestación reproducible y ejecución experimental, no en añadir otra salida
puntual a las clases.

### Puerta de entrada a los experimentos finales

Los experimentos definitivos se ejecutarán **después** de completar las
correcciones, reorganización, API y capacidades experimentales del plan, no
durante el refactor. Antes sólo se permiten tests y ejecuciones pequeñas de
verificación. Un piloto acotado estimará tiempo y memoria; si obliga a corregir
el código, se repetirá la validación antes de fijar versión, configuración y
semillas para el estudio final. Estos incrementos no dan por completada una
fase transversal ni certifican un paquete publicable.

## Decisiones estadísticas previas al refactor

### Priors e hiperparámetros

Los modos actuales `OLS` y `NON_INFORMATIVE` utilizan estadísticas calculadas
con los mismos datos. Ambos son procedimientos empirical Bayes; el segundo no
es una prior no informativa.

La interfaz objetivo será:

```r
prior = "ols_empirical_bayes"
prior = "weakly_informative"
prior = list(...)  # especificación explícita
```

`ols_empirical_bayes` se conservará para inicialización, depuración y
comparaciones computacionales. Debe documentarse que reutiliza los datos y no
representa una prior completamente especificada antes de observarlos.

`weakly_informative` deberá ser independiente de OLS. Su función es regularizar
valores extremos sin concentrar excesivamente la masa dentro del rango
plausible. No existe un preset universal si `X` y `y` permanecen en escalas
arbitrarias: las escalas deben justificarse usando conocimiento sustantivo o
una convención de unidades documentada. Si se usan desviaciones observadas para
calibrarlas, el nombre correcto sería `data_scaled_empirical_bayes`.

La opción recomendada para resultados definitivos es una prior explícita. OLS
puede utilizarse separadamente como inicialización:

```r
model <- chlrm(
  prior = list(...),
  initial_values = "ols"
)
```

Para la parametrización conjugada actual se debe garantizar:

- matrices de covarianza simétricas definidas positivas;
- grados inverse-Wishart que produzcan los momentos que se desean reportar;
- shapes y rates Gamma/inverse-gamma positivos;
- `nu0 > 2` si se pretende una media finita para la prior de la varianza;
- escalas de coeficientes con unidades compatibles con cada columna de `X`;
- análisis de sensibilidad para las conclusiones del caso Farms.

### Intervalos de VI y SVI

Los intervalos producidos por VI y SVI son **intervalos creíbles
variacionales**, no intervalos de confianza frecuentistas ni intervalos de la
posterior exacta.

Si

$$
q(\beta)=N(\mu_\beta,\Sigma_\beta),
$$

el intervalo marginal del 95 % para el coeficiente $r$ es

$$
\mu_{\beta_r}\pm z_{0.975}\sqrt{\Sigma_{\beta,rr}}.
$$

Este cálculo es exacto dentro del factor normal $q(\beta)$. La limitación es la
aproximación mean-field: si omite dependencias posteriores importantes puede
subestimar la incertidumbre, aun cuando los cuantiles estén bien calculados.

Para

$$
q(\sigma^2)=IG(a,b),
$$

con $1/\sigma^2\sim Gamma(a,\text{rate}=b)$, los límites correctos son

$$
L=\frac{1}{Q_{\Gamma}(0.975)},\qquad
U=\frac{1}{Q_{\Gamma}(0.025)}.
$$

En CHLRM, si $q(\omega)=Dirichlet(\alpha)$, el intervalo marginal para
$\omega_k$ se obtiene de

$$
\omega_k\sim Beta(\alpha_k,\alpha_0-\alpha_k).
$$

Para la matriz `Sigma` inverse-Wishart se usarán draws del factor variacional y
cuantiles empíricos de varianzas, desviaciones, correlaciones o elementos de
interés. No se aplicará una aproximación normal elemento a elemento.

Los intervalos predictivos se construirán muestreando desde la aproximación
variacional completa y después desde la distribución observacional. En CHLRM
se debe propagar también la incertidumbre de asignación. Los intervalos por
componente son condicionales a la etiqueta elegida; deben acompañarse con
resúmenes invariantes a etiquetas, como coclustering y entropía de asignación.

### `sample_rhoj_chlrm`

La implementación calcula las probabilidades condicionales de asignación y
muestrea un valor discreto `z[j]`. Por tanto, el kernel no muestrea `rho`; el
problema es semántico, no algebraico.

Se renombrará como:

```r
sample_allocations_chlrm()
```

La convención será:

- `allocation`: asignación discreta $z_j$;
- `responsibility`: probabilidad variacional $q(z_j=k)$;
- `learning_rate`: paso $\rho_t$ de SVI.

### `nj`

$n_j$ es el número de observaciones del grupo $j$. Es necesario en el
normalizador de la likelihood, las formas posteriores de las varianzas, los
conteos efectivos por componente, el ELBO y la generación predictiva.

No debe suministrarse repetidamente como argumento público cuando puede
derivarse de los datos. Se calculará una vez:

```r
n_group <- lengths(y)
```

Se validará contra `nrow(X[[j]])` y se guardará en el objeto de datos. Los
kernels internos podrán reutilizarlo. El argumento `nj` se retirará de
funciones como `dic_chlrm()`, donde actualmente no se utiliza, y de funciones
que pueden derivarlo inequívocamente.

## DIC y WAIC

### Evidencia de las referencias primarias

La revisión de `spiegelhalter2002.pdf`, `watanabe10a.pdf`,
`vehtari2016.pdf` y `watanabe2021.pdf` fortalece las decisiones de esta
sección:

1. **DIC depende del foco jerárquico.** Spiegelhalter et al. muestran que en
   un modelo jerárquico la likelihood, la complejidad efectiva y DIC no quedan
   definidos de forma única hasta decidir qué nivel se mantiene fijo y qué
   parámetros se integran. Por ello, el DIC condicional en $z_j$ y el DIC que
   marginaliza $z_j$ responden preguntas distintas y no deben compartir
   nombre ni compararse directamente.
2. **La justificación original de DIC es aproximada.** El propio artículo lo
   presenta tentativamente y con una justificación de decisión heurística.
   Además, $p_D$ depende de la parametrización y puede ser negativo cuando la
   media posterior representa mal una posterior multimodal.
3. **CHLRM es singular.** Watanabe identifica explícitamente las mezclas de
   distribuciones y los modelos con variables ocultas como modelos singulares:
   distintos parámetros pueden inducir la misma distribución y la información
   de Fisher puede ser degenerada. El label switching de CHLRM es una
   manifestación directa de esa no identificabilidad.
4. **En modelos singulares DIC y WAIC dejan de ser equivalentes.** Watanabe
   demuestra que WAIC y validación cruzada bayesiana son asintóticamente
   equivalentes, incluso bajo singularidad, mientras los criterios DIC no
   tienen en general el mismo comportamiento que el error de generalización.
   Esto convierte a DIC en una referencia descriptiva, no en evidencia para
   seleccionar $K$ en CHLRM.
5. **La ventaja de WAIC es asintótica, no garantía finita.** Vehtari, Gelman y
   Gabry muestran que PSIS-LOO es más robusto que WAIC con priors débiles u
   observaciones influyentes. Recomiendan desconfiar de WAIC cuando alguna
   contribución satisface
   $\operatorname{Var}_s[\log p(y_i\mid\theta^{(s)})]>0.4$ y pasar a LOO o
   K-fold cuando los diagnósticos fallen.
6. **La unidad pointwise forma parte del estimando.** Las fórmulas de WAIC y
   LOO presuponen una factorización explícita de la likelihood. En CHLRM la
   finca completa debe ser la fila pointwise si el objetivo es predecir una
   unidad jerárquica completa; tratar observaciones internas como independientes
   responde a otra pregunta.
7. **En mezclas debe utilizarse la likelihood observada, no la completada.**
   Watanabe (2021) distingue explícitamente
   $\operatorname{Var}[\log p(X_i\mid w)]$ de
   $\operatorname{Var}[\log p(X_i\mid w,Y_i)]$ y advierte que confundirlas es
   un error común en experimentos numéricos. En CHLRM, $Y_i$ corresponde a la
   asignación latente. Esto confirmó que el antiguo WAIC basado en
   `Betaj` y `sigma2j` no era el WAIC de la mezcla marginal; esa ruta fue retirada.
8. **WAIC no garantiza recuperar el número mínimo de componentes.** En una
   mezcla singular, añadir componentes redundantes puede cambiar muy poco la
   pérdida de generalización. Watanabe (2021) muestra que WAIC puede seleccionar
   modelos redundantes con mayor frecuencia que en modelos regulares. Por ello,
   buen WAIC predictivo no implica recuperación correcta de $K$.
9. **La prior Dirichlet afecta la redundancia.** El comportamiento de pesos de
   componentes sobrantes cambia con `alpha0` y con el tamaño muestral; también
   existe una transición análoga bajo mean-field VI. La selección de $K$ debe
   incluir sensibilidad a `alpha0`, número de componentes ocupados y estabilidad
   del coclustering.
10. **WAIC y WBIC tienen objetivos distintos.** WAIC estima pérdida de
    generalización predictiva; WBIC aproxima la energía libre o evidencia
    marginal y requiere una posterior templada con
    $\beta=1/\log n$. WBIC puede considerarse en un experimento separado de
    selección de $K$, pero no debe añadirse reutilizando draws posteriores
    ordinarios ni presentarse como sustituto de WAIC.

Consecuencia: en los resultados principales se priorizarán MLPD de validación
cruzada y, cuando sea viable, PSIS-LOO por finca. WAIC marginal por finca será
un estimador complementario con diagnóstico. DIC marginal se conservará sólo
para continuidad con la literatura y siempre marcado como secundario. La
recuperación de $K$ se evaluará aparte de la capacidad predictiva; opcionalmente
se estudiará WBIC en MCMC como experimento metodológico independiente.

### Qué pregunta responde DIC

Sea la devianza

$$
D(\theta)=-2\log p(y\mid\theta).
$$

La implementación utiliza

$$
\bar D=E[D(\theta)\mid y],\qquad
p_D=\bar D-D(\bar\theta),
$$

y

$$
DIC=D(\bar\theta)+2p_D=\bar D+p_D.
$$

DIC fue propuesto como un compromiso entre ajuste y complejidad efectiva. Su
principal atractivo es que puede calcularse fácilmente a partir de draws
posteriores. Sin embargo, depende de representar el modelo mediante un único
valor posterior medio $\bar\theta$.

#### Problemas de DIC

1. **Depende de la parametrización.** Transformaciones equivalentes de los
   parámetros pueden cambiar $D(\bar\theta)$ y $p_D$.
2. **El promedio puede no representar la posterior.** En distribuciones
   asimétricas o multimodales, $\bar\theta$ puede caer en una región de baja
   densidad.
3. **La complejidad efectiva puede ser negativa o inestable.** Esto sucede en
   modelos débiles, jerárquicos y con parámetros latentes.
4. **Es problemático en mezclas.** El promedio de parámetros con etiquetas que
   cambian no representa un componente. La ruta CHLRM retirada usaba `Betaj`
   y `sigma2j`, invariantes ante permutación conjunta, pero condicionaba las
   asignaciones. La ruta nueva es marginal y promedia componentes alineados
   en MCMC; ese alineamiento no elimina los demás problemas de DIC ni lo
   vuelve evidencia principal para seleccionar $K$.
5. **No es pointwise.** No ofrece directamente contribuciones por unidad ni un
   error estándar de comparación.
6. **VI y MCMC no usan la misma distribución.** En MCMC la esperanza se estima
   bajo la posterior muestreada; en VI/SVI se estima bajo $q$. Los números son
   calculables con la misma fórmula, pero la diferencia puede reflejar tanto el
   modelo como el error de aproximación variacional.

Decisión: conservar DIC únicamente como diagnóstico histórico secundario,
etiquetado como `dic_conditional`. No usarlo para seleccionar $K$, afirmar
superioridad predictiva o comparar modelos con likelihood/unidad predictiva
distinta.

### Qué pregunta responde WAIC

Para cada unidad predictiva $i$ y draws $\theta^{(s)}$:

$$
\operatorname{lppd}_i=
\log\left(\frac{1}{S}\sum_{s=1}^{S}p(y_i\mid\theta^{(s)})\right),
$$

$$
p_{WAIC,i}=\operatorname{Var}_s
\left(\log p(y_i\mid\theta^{(s)})\right),
$$

$$
\operatorname{elpd}_{WAIC}=\sum_i
(\operatorname{lppd}_i-p_{WAIC,i}),\qquad
WAIC=-2\operatorname{elpd}_{WAIC}.
$$

WAIC aproxima desempeño predictivo fuera de muestra y, a diferencia de DIC,
usa toda la distribución de la log-likelihood pointwise. Fue formulado para
modelos regulares y singulares; bajo condiciones asintóticas se relaciona con
leave-one-out cross-validation.

#### Problemas de WAIC

1. **La unidad predictiva define el estimando.** Una fila, un grupo completo y
   una serie temporal no son unidades intercambiables.
2. **En HLRM el WAIC actual sigue siendo condicional.** Incluso `unit="group"`
   suma la log-likelihood dentro del grupo, pero usa efectos locales o
   asignaciones aprendidos con ese mismo grupo. No evalúa una finca nueva.
   CHLRM ya integra las etiquetas en su nuevo GOF marginal.
3. **No elimina automáticamente leakage jerárquico.** Para generalizar a nuevos
   grupos se necesita LOGO-CV o una densidad marginal que integre efectos y
   asignaciones todavía no observados.
4. **Puede ser inestable con observaciones influyentes o priors débiles.** Una
   varianza grande de la log-likelihood pointwise indica que la aproximación
   puede ser poco fiable. La implementación advierte cuando
   $p_{WAIC,i}>0.4$, pero esto es una señal diagnóstica, no una corrección.
5. **VI/SVI pueden subestimar $p_{WAIC}$.** Si $q$ subestima varianzas
   posteriores, también puede reducir artificialmente la varianza de la
   log-likelihood y producir un WAIC demasiado optimista.
6. **Comparabilidad limitada.** Sólo pueden compararse valores calculados sobre
   exactamente los mismos datos, la misma unidad predictiva y la misma
   definición de likelihood. El valor absoluto no tiene interpretación aislada.

### Qué calcula actualmente el repositorio

| Modelo | Métrica | Unidad/condicionamiento | Uso recomendado |
|---|---|---|---|
| LRM | `waic` | observación | comparación predictiva sobre las mismas filas |
| HLRM | `waic` | observación condicional en efectos de grupo | diagnóstico in-sample condicional |
| HLRM | `waic_group` | grupo observado, con efecto aprendido | dependencia intragrupo; no grupo nuevo |
| CHLRM | `waic` | observación condicional en asignación/efecto | diagnóstico condicional |
| CHLRM | `waic_group` | grupo observado, asignación aprendida | no equivale a LOGO-CV |
| Todos | `dic` | devianza total en medias posteriores | referencia histórica secundaria |

`waic_lrm()` debe delegar en el helper común `waic_from_loglik()` para devolver
pointwise, error estándar y warnings igual que HLRM y CHLRM.

### Contraste con `Computacion_revision_paper_Camilo.pdf`

Las secciones 22 y 23 del documento de revisión precisan el estimando que debe
usarse para Farms. La unidad predictiva es la finca completa y la asignación
latente de la finca debe marginalizarse. Para cada draw $s$ se requiere

$$
\ell_j^{(s)}=
\log\left\{
\sum_{k=1}^K \omega_k^{(s)}
\prod_{i=1}^{n_j}
\mathcal N\left(y_{ij}\mid x_{ij}^{\top}\beta_k^{(s)},
\sigma_k^{2(s)}\right)
\right\}.
$$

Equivalentemente, y de forma numéricamente estable,

$$
\ell_j^{(s)}=
\operatorname{logsumexp}_{k}
\left[
\log\omega_k^{(s)}+
\sum_{i=1}^{n_j}\log
\mathcal N\left(y_{ij}\mid x_{ij}^{\top}\beta_k^{(s)},
\sigma_k^{2(s)}\right)
\right].
$$

Esto **no era lo que calculaba** `waic_chlrm(..., unit="group")`, función ya
retirada. Aquella implementación tomaba `Betaj` y `sigma2j`, es decir, los parámetros del
componente asignado a cada finca en cada draw, y después suma las
log-verosimilitudes de sus observaciones. Calculaba una medida
condicional a $z_j$; no usaba `omega` ni marginalizaba la asignación. La nueva
ruta sí aplica la fórmula anterior. Agrupar las filas con `rowsum()` cambia la unidad de
WAIC, pero no convierte una likelihood condicional en una likelihood marginal.

| Aspecto | Documento de revisión | Código actual | Acción |
|---|---|---|---|
| Unidad WAIC | finca completa | finca completa | integrar en runner final |
| Asignación $z_j$ | marginalizada | marginalizada | probar en estudio final |
| Parámetros usados | `omega`, `Betak`, `sigma2k` | mismo contrato compacto | mantener pruebas |
| VI/SVI | 5000 draws de $q(\Theta)$ | GOF marginal integrado y S configurable | configurar 5000 draws en el runner |
| DIC | secundario, misma likelihood por finca | DIC marginal separado y opcional | no seleccionar K sólo con DIC |
| Media para DIC MCMC | después de relabeling | componentes alineados iterativamente | usar draws combinados tras implementar multichain |

#### Implementación propuesta para CHLRM

1. Crear `loglik_group_marginal_chlrm(y_list, X_list, draws)`, que devuelva una
   matriz $m\times S$ y use `logsumexp` sobre los $K$ componentes.
2. Exigir en `draws` las matrices/arreglos `omega` ($S\times K$), `Betak`
   ($S\times K\times p$) y `sigma2k` ($S\times K$).
3. Extender `vi_sample_draws_chlrm()` para generar
   $\omega^{(s)}\sim Dirichlet(\alpha)$. Para WAIC no generar ni utilizar
   $z_j^{(s)}\sim Categorical(\rho_j)$.
4. Calcular `waic_group_marginal_chlrm()` pasando esa matriz al helper común
   `waic_from_loglik()` y devolver `elpd`, error estándar, contribuciones por
   finca y advertencias.
5. Mantener temporalmente el cálculo existente con el nombre explícito
   `waic_group_conditional_chlrm()` para reproducibilidad, sin presentarlo como
   el WAIC principal del artículo.
6. Implementar `dic_group_marginal_chlrm()` con la misma likelihood marginal.
   En MCMC, calcular $\bar\Theta$ únicamente desde cadenas relabeladas; en
   VI/SVI usar las medias de los factores variacionales y documentar que el
   resultado depende de $q$, no de la posterior exacta.
7. Retirar el argumento público `nj`: los tamaños deben derivarse de
   `lengths(y_list)` y validarse contra las filas de `X_list`.

#### Pruebas mínimas de corrección

- Para $K=1$, la likelihood marginal debe coincidir con la likelihood agrupada
  del único componente.
- Una permutación conjunta de `omega`, `Betak` y `sigma2k` no debe cambiar
  ninguna contribución por finca, WAIC ni DIC.
- La implementación debe permanecer finita con probabilidades de mezcla muy
  pequeñas; la suma debe hacerse en escala logarítmica.
- MCMC, VI y SVI deben pasar exactamente la misma estructura $m\times S$ al
  cálculo de WAIC.
- Un ejemplo pequeño debe contrastarse contra un cálculo directo de la mezcla
  en escala natural.
- El DIC marginal debe comprobarse antes y después del relabeling: la devianza
  de cada draw es invariante, pero $D(\bar\Theta)$ sólo es interpretable cuando
  las etiquetas se alinearon.

#### Decisión para el artículo

Adoptar el WAIC marginal por finca definido en la revisión para la tabla de
Farms. Reportar además `lppd`, $p_{WAIC}$, `elpd_waic`, su error estándar y el
número de fincas con $p_{WAIC,j}>0.4$. Conservar DIC en una columna secundaria
y acompañarlo de $\bar D$ y $p_D$; no usarlo para escoger $K$ ni para concluir
superioridad predictiva. La selección de $K$ debe basarse principalmente en
MLPD de validación cruzada y estabilidad del clustering, tal como propone la
sección 21 del documento. También debe reportar componentes ocupados y
sensibilidad a `alpha0`: WAIC puede preferir una mezcla predictivamente adecuada
pero redundante, lo cual no constituye un error del criterio ni demuestra que
ese $K$ sea identificable.

El WAIC por finca y el MLPD-CV de la sección 21.4 tampoco responden exactamente
la misma pregunta. El primero toma una finca completa como unidad; la CV
propuesta deja una observación de cada finca en cada fold y predice dentro de
fincas ya observadas. Ambos pueden reportarse, pero deben llevar etiquetas
distintas y no interpretarse como dos estimadores intercambiables del mismo
riesgo predictivo. Si interesa una finca completamente nueva, hace falta además
LOGO-CV por finca.

### Recomendación de evaluación predictiva

1. Reportar `elpd` además de la escala deviance.
2. Comparar diferencias pointwise con su error estándar, no sólo dos números.
3. Mantener `waic_observation_conditional` y
   `waic_group_conditional` con nombres explícitos.
4. Para una nueva finca usar LOGO-CV marginalizando parámetros locales y, en
   CHLRM, la asignación.
5. Preferir PSIS-LOO cuando las unidades sean aproximadamente independientes y
   revisar sus diagnósticos Pareto-$k$.
6. Usar K-fold/LOGO refit cuando PSIS-LOO no sea fiable o la estructura de
   dependencia requiera dejar bloques completos fuera.
7. No usar DIC como criterio principal de selección de $K$.

## Plan de ejecución histórico (sustituido)

Este apartado conserva el plan de trabajo acordado y lo actualiza con las
decisiones anteriores. Las fases están ordenadas por dependencia: primero se
fijan y verifican los contratos estadísticos; después se reorganiza y optimiza
el código.

### Fase 0: congelar comportamiento y notación

- Crear pruebas de caracterización con semillas fijas para LRM, HLRM y CHLRM.
- Guardar dimensiones, nombres, parámetros y métricas canónicas, no objetos
  completos dependientes de detalles internos.
- Crear una tabla símbolo teórico–nombre de código–dimensión–distribución.
- Añadir pruebas específicas de LRM, actualmente ausentes.
- Probar componentes vacíos, $K=1$, $K=m$, colinealidad y datos no finitos.

Criterio de aceptación: resultados base reproducibles y mapa de notación
revisado antes de mover funciones.

### Fase 1: seguridad y contratos

- Reemplazar `logdet()` por una implementación Cholesky que exija SPD.
- Implementar validadores comunes de datos, hiperparámetros, MCMC y VI/SVI.
- Validar el estado requerido por `fit()`, `summarize_posterior()` y GOF.
- Separar media ajustada, predicción de respuestas y grupos nuevos.
  Implementado mediante `predict()`/`posterior_predict()` y `group_ids`;
  queda pendiente decidir si se expone además un verbo `fitted()` de conveniencia.
- Definir un contrato RNG común; la inicialización VI no debe consumir
  aleatoriedad de manera opaca.
- Usar políticas distintas de parada CAVI/SVI. Implementado para CHLRM;
  pendiente armonizar LRM/HLRM y los runners sin confundir épocas y checkpoints.

Criterio de aceptación: inputs inválidos fallan antes de los loops con mensajes
claros; ninguna métrica depende de una llamada previa no documentada.

### Fase 2: derivación y coherencia bayesiana

- Derivar cada condicional MCMC y cada factor CAVI/SVI desde la misma
  factorización.
- Verificar el ELBO término por término.
- Demostrar o redefinir el esquema global de SVI CHLRM.
- Separar prior de inicialización OLS.
- Introducir `ols_empirical_bayes`, `weakly_informative` y prior explícita.
- Añadir análisis de sensibilidad de prior.
- Corregir la relación entre Dahl y parámetros relabelados.

Criterio de aceptación: reducción SVI a CAVI con minibatch completo y
`learning_rate=1`, más tests numéricos de cada identidad.

### Fase 3: criterios predictivos e intervalos

- Unificar la construcción de log-likelihood pointwise.
- Hacer que LRM use `waic_from_loglik()`.
- Implementar para CHLRM la log-likelihood marginal por finca con `omega`,
  `Betak` y `sigma2k`, y usarla tanto en WAIC como en DIC marginal.
- Extender los draws VI/SVI con muestras Dirichlet de `omega`.
- Renombrar WAIC según unidad, condicionamiento y marginalización.
- Mantener el DIC condicional sólo para reproducibilidad; reportar como
  secundario el DIC marginal definido en el documento de revisión.
- Corregir el warning CHLRM: la devianza marginal por draw es label-invariant,
  pero la media paramétrica de DIC requiere relabeling.
- Implementar LOGO-CV para grupos nuevos.
- Separar LOGO-CV de la CV dentro de finca descrita en la sección 21.4.
- Unificar PPC y añadir discrepancias residuales y por grupo. El generador y el
  consumo de `y_rep` ya son comunes; faltan las discrepancias T1--T5.
- Implementar intervalos creíbles variacionales mediante un módulo común.
  Implementado para predicciones puntuales de la media.
- Añadir intervalos predictivos mediante draws de MCMC o $q$. Implementado con
  `posterior_predict()`; faltan estudios de cobertura/calibración.

Criterio de aceptación: cada métrica documenta datos, unidad predictiva,
condicionamiento y distribución usada para integrar incertidumbre.

### Fase 4: reorganización y nomenclatura

La reorganización tendrá dos capas con responsabilidades distintas:

1. **Paquete R:** implementación reusable, documentada y probada del método.
2. **Suite de validación científica:** E1--E5, benchmarks y construcción del
   manuscrito, consumiendo exclusivamente la API pública del paquete.

No todo el experimento pesado debe instalarse con la librería. La librería debe
hacer posible el experimento; la suite debe demostrar que sus capacidades
funcionan y caracterizar estadísticamente el método.

Arquitectura objetivo compatible con un paquete R:

```text
.
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── R/
│   ├── core-numerics.R
│   ├── core-validation.R
│   ├── core-diagnostics.R
│   ├── core-summaries.R
│   ├── lrm-model.R
│   ├── lrm-simulation.R
│   ├── lrm-inference.R
│   ├── hlrm-model.R
│   ├── hlrm-simulation.R
│   ├── hlrm-inference.R
│   ├── chlrm-model.R
│   ├── chlrm-simulation.R
│   ├── chlrm-priors.R
│   ├── chlrm-mcmc.R
│   ├── chlrm-vi.R
│   ├── chlrm-partitions.R
│   └── predictive-assessment.R
├── man/
├── tests/testthat/
├── vignettes/
├── inst/extdata/
└── thesis/final/code/protocol/
    ├── config/
    ├── E1-monte-carlo/
    ├── E2-select-k/
    ├── E3-svi-sensitivity/
    ├── E4-scalability/
    └── E5-farms/
```

- Eliminar `source()` y `library()` de archivos internos.
- Convertir el núcleo en paquete R; no mantener un loader paralelo como API.
- Mantener funciones puras para kernels y clases únicamente como
  coordinadoras de estado.
- Introducir nombres snake_case con wrappers temporales para la API anterior.
- Renombrar `sample_rhoj_chlrm()` como `sample_allocations_chlrm()`.
- Calcular `n_group` una vez y eliminar argumentos `nj` redundantes.
- Documentar la API exportada con roxygen2 y mantener kernels internos sin
  exportar.
- Declarar dependencias mínimas en `Imports` y herramientas experimentales en
  `Suggests`; ningún ejemplo instalado debe iniciar cómputo pesado.
- Añadir vignettes pequeñas para simulación, ajuste, diagnóstico, predicción y
  selección de K, todas con perfiles seguros.
- Establecer una política de compatibilidad y deprecación para nombres actuales.

Criterio de aceptación: `R CMD check` pasa sin errores ni warnings relevantes;
los tests y ejemplos no dependen del directorio de trabajo; el protocolo smoke
se ejecuta instalando el paquete y usando sólo funciones exportadas.

### Fase 5: performance y validación experimental

- Sustituir inversas explícitas por Cholesky y `solve(A, b)`.
- Cachear factorizaciones y estadísticas suficientes constantes.
- Procesar WAIC y posterior predictiva por bloques para evitar matrices `n x S`
  innecesarias.
- Parametrizar completamente los simuladores y retirar diseños favorables
  hardcodeados.
- Ejecutar un diseño factorial con separación, desbalance, componentes vacíos,
  distintos $m,n_j,p,K$, covariables correlacionadas y múltiples semillas.
- Comparar recuperación, cobertura, predicción, ELBO, ARI/coclustering, tiempo
  y memoria con errores Monte Carlo.

Criterio de aceptación: igualdad numérica dentro de tolerancia antes/después de
optimizar y benchmarks reproducibles con configuración y entorno registrados.

## Capacidad objetivo histórica (sustituida)

Al terminar el refactor, una persona debe poder ejecutar:

```r
run_protocol(
  experiment = "E1",
  profile = "smoke",
  workers = 2,
  memory_gb = 6,
  resume = TRUE
)
```

sin editar las clases estadísticas ni copiar funciones en scripts. La misma
interfaz debe aceptar `E1` a `E5` y los perfiles:

- `smoke`: contratos y extremo a extremo con tamaños mínimos;
- `pilot`: todos los escenarios con pocas réplicas para estimar costo;
- `full`: tamaños exactos del protocolo, sólo por solicitud explícita.

Cada experimento será una configuración declarativa. Los módulos de `R/` deben
proveer generación, prior, inicialización, ajuste, draws, relabeling, predicción,
métricas y diagnósticos. Los runners sólo orquestarán esos módulos y escribirán
resultados; no contendrán nuevas fórmulas estadísticas.

### Definición de terminado

El código estará listo para la revisión cuando:

1. todas las filas bloqueantes de
   `paper/EVALUATION_PLAN.md` estén en `Implementado`;
2. E1--E5 completen el perfil `smoke` desde un entorno limpio;
3. una ejecución interrumpida pueda reanudarse sin repetir tareas terminadas;
4. cada número de una tabla o figura pueda rastrearse hasta configuración,
   semilla, resultado crudo y versión del código;
5. `full` no pueda iniciarse accidentalmente sin estimación previa de recursos y
   confirmación explícita;
6. los estimandos marginales/condicionales y las unidades predictivas aparezcan
   en nombres y metadatos, no sólo en documentación narrativa.

## Presupuesto y seguridad computacional

El protocolo completo contiene decenas de miles de ajustes. La eficiencia no
se resolverá únicamente vectorizando kernels: también requiere controlar el
trabajo total, la memoria simultánea y la persistencia.

### Unidad de trabajo y reanudación

La tarea mínima será
`experimento/escenario/replica/metodo/inicio-o-cadena`. Cada tarea:

- recibe una semilla derivada determinísticamente de su identificador;
- escribe primero a un archivo temporal y lo publica de forma atómica;
- registra estado `pending/running/complete/failed`;
- puede repetirse aisladamente;
- no se vuelve a ejecutar con `resume=TRUE` si su salida pasó validación.

El agregador sólo leerá resultados terminados. Nunca volverá a ajustar modelos.

### Límites por defecto

- `workers=1` por defecto y paralelismo a nivel de tareas, no dentro del mismo
  objeto mutable.
- Límite configurable de memoria; el scheduler no inicia una tarea si la suma
  estimada excede el presupuesto.
- BLAS/OpenMP fijados a un thread por worker para evitar sobre-suscripción.
- Draws y log-likelihood se procesan por bloques; no se materializan arreglos
  `replicas x draws x groups x components x observations`.
- PSM `m x m` sólo se construye cuando es requerida y nunca para los tamaños de
  escalabilidad si basta una métrica streaming.
- MCMC se almacena por bloques en disco; las cuatro cadenas no se duplican al
  combinarlas o relabelarlas.
- Los objetos de ajuste completos son opcionales; siempre se guardan resúmenes,
  diagnósticos y draws mínimos necesarios para auditoría.

### Estimación previa y corte seguro

El perfil `pilot` debe medir tiempo y memoria máximos por tipo de tarea. Antes de
`full`, el runner proyectará horas de CPU, tiempo de pared y almacenamiento. Si
supera los límites configurados, abortará antes de ajustar y mostrará qué parte
del diseño causa el exceso. También debe aceptar límites como:

```r
limits = list(
  max_wall_hours = 8,
  max_memory_gb = 6,
  max_storage_gb = 20,
  max_workers = 2
)
```

Al alcanzar un límite durante la ejecución, se completa la tarea activa, se
guarda checkpoint y se detiene limpiamente. Esto permite dividir `full` en
varias sesiones sin perder resultados.

## Correspondencia entre fases y protocolo

| Fase del plan | Filas de revisión que debe cerrar | Entregable verificable |
|---|---|---|
| 0 | infraestructura de tests y notación | tests de caracterización y schema versionado |
| 1 | contratos, RNG, estado y predicción | validadores, estados serializables y API train/test |
| 2 | prior, CAVI, SVI y relabeling | derivaciones, prior del protocolo y equivalencias numéricas |
| 3 | métricas, WAIC, DIC, PPC y CV | módulo único de evaluación con tests manuales |
| 4 | arquitectura y nombres | paquete/loader sin `source()` internos y runners declarativos |
| 5 | E1--E5, performance y reporte | perfiles smoke/pilot/full, checkpoints, tablas y figuras |

El detalle de escenarios, tamaños y salidas permanece en la matriz de
trazabilidad; no se duplicará como constantes dispersas en los scripts.

## Orden inmediato histórico (sustituido)

1. Completar estado serializable/reanudación MCMC y extensión incremental por
   bloques sobre el multichain ya diagnosticado.
2. Cerrar evaluación pura: cobertura/fidelidad/CPE/MLPD/regret, T1--T5,
   MCSE dependiente y LOGO/CV con ejemplos manuales.
3. Uniformar contratos de ajuste LRM/HLRM con CHLRM y completar validadores,
   estado serializable y reanudación.
4. Convertir el núcleo en estructura de paquete, retirar `source()` internos y
   hacer que la suite científica consuma sólo la API pública.
5. Construir runner declarativo con perfiles `smoke/pilot/full`, presupuesto,
   checkpoints, procedencia y agregación automática.
6. Ejecutar smoke, luego piloto; congelar código/configuración y sólo entonces
   ejecutar E3, E1, E2, E4 y E5.

No conviene empezar por mover archivos o renombrar todos los parámetros: eso
dificultaría distinguir cambios semánticos de cambios puramente estructurales.

## Séptimo incremento: MCMC multichain diagnosticado

- Se añadió `run_mcmc_chains()` para LRM, HLRM y CHLRM, con número de cadenas
  configurable, semillas distintas y ejecución secuencial por seguridad local.
- La semilla abarca inicialización y ajuste completos; esto corrigió la falta
  de reproducibilidad bit a bit causada por inicializaciones aleatorias en los
  constructores LRM/HLRM.
- `posterior` calcula R-hat rank-normalizado, bulk ESS y tail ESS; ESS/segundo
  usa la suma de tiempo de las cadenas ejecutadas secuencialmente.
- CHLRM aprende una sola referencia de etiquetas con los draws de todas las
  cadenas y alinea cada cadena antes de diagnosticar.
- `combine_mcmc_chains()` es una operación separada y rechaza por defecto las
  cadenas que incumplen los umbrales declarados. Los draws combinados conservan
  cadena, iteración y semilla de origen.
- El runner exige `thin=1` y estima memoria antes de comenzar. No promete un
  límite duro de RAM ni extensión incremental: checkpoint/reanudación requieren
  primero uniformar el estado mutable de las tres clases.
- `tests/run-validation.R` deja reportes locales livianos bajo `tests/output/`;
  nunca persiste cadenas completas automáticamente.

## Octavo incremento: VI/SVI multistart comparable

- LRM y HLRM adoptaron el mismo contrato de estado, RNG y monitor ELBO de
  CHLRM: cambio relativo, paciencia, razón de parada y error explícito.
- `make_vi_initial_states()` genera estados serializables independientes de la
  aleatoriedad posterior del algoritmo. Un mismo conjunto puede reutilizarse
  entre VI y SVI.
- `run_vi_multistart()` ejecuta secuencialmente, conserva todos los resultados,
  errores, ELBO, iteraciones, tiempos y semillas, y no selecciona ajustes no
  convergentes salvo solicitud explícita exploratoria.
- La selección es el mayor ELBO terminal finito entre intentos elegibles; los
  empates se resuelven por el primer inicio para reproducibilidad.
- CHLRM alinea todos los inicios a una referencia suministrada —preferiblemente
  la MCMC del dataset— o aprende una referencia común entre inicios cuando no se
  suministra una.
- Validación de enteros y derivación de semillas se movieron a
  `inference-control.R`; VI no depende de helpers internos de MCMC.

## Noveno incremento: checkpoint y extensión MCMC exacta

- LRM, HLRM y CHLRM MCMC comparten estados
  `initialized/running/fitted/failed`, RNG aislado y metadatos de ejecución.
- `resume_mcmc()` continúa desde los parámetros Gibbs y estado RNG terminales,
  no repite warm-up y sólo publica el nuevo bloque cuando termina completo.
- Reanudar y ejecutar de manera continua producen draws y parámetros terminales
  idénticos en los tres modelos; esta identidad tiene pruebas de regresión.
- `mcmc_checkpoint()` usa una lista plana sin entornos R6;
  `save_mcmc_checkpoint()` publica mediante archivo temporal y rename atómico.
- Resúmenes, GOF y relabeling se invalidan al añadir draws para impedir el uso
  accidental de resultados obsoletos.
- `run_mcmc_chains()` acepta `block_size/max_samples`, recalcula R-hat y ESS por
  bloque y registra la trayectoria diagnóstica y razón de parada.
- El paquete no escribe checkpoints multichain automáticamente: el scheduler
  del paper asignará una ruta por cadena/tarea para controlar reanudación y
  procedencia sin imponer políticas de archivos a la API estadística.

## Décimo incremento: contratos matriciales y Cholesky

- Se retiró `.sym()`: promediar siempre una matriz con su traspuesta podía
  ocultar errores algebraicos. `symmetrize_numeric()` sólo corrige discrepancias
  compatibles con redondeo y rechaza asimetría material.
- `chol_spd()`, `solve_spd()` y `logdet_spd()` centralizan los sistemas,
  inversas y log-determinantes de covarianzas, precisiones y escalas. No quedan
  llamadas directas a `solve()`, `logdet()` o `.sym()` dentro de `R/`.
- Los constructores de LRM, HLRM y CHLRM, tanto MCMC como VI/SVI, validan sus
  hiperparámetros antes de inicializar el algoritmo y guardan la versión
  canónica validada.
- Las inversas opcionales se contrastan contra su matriz original en lugar de
  confiar en dos entradas potencialmente contradictorias. `n0 > p - 1` exige
  una inverse-Wishart propia sin afirmar que todos sus momentos existan.
- Los presets OLS agregan regularización mínima y mantienen varianza positiva
  aun con ajuste perfecto. Esto garantiza un objeto numérico válido, pero no
  elimina su carácter empirical Bayes.
- Las pruebas comparan la nueva álgebra con cálculos base, prueban matrices
  indefinidas/asimétricas e hiperparámetros contradictorios y cubren el caso de
  ajuste OLS perfecto.

## Undécimo incremento: contratos públicos uniformes

- `fit()` usa `fit-result-1` en LRM, HLRM y CHLRM para MCMC, VI y SVI. Sus
  bloques fijos son parámetros, variables latentes, objetivo, diagnósticos,
  tiempos y metadatos. ELBO/log-likelihood/log-joint dejaron de mezclarse con
  parámetros.
- `summarize_posterior()` usa `posterior-summary-1` y nombres paramétricos
  canónicos específicos de cada modelo. Cada resumen usa `mean`, `sd`, `lower`, `upper`,
  `level`, `distribution`, `mean_exists` y `contains_zero`.
- `predict()` y `posterior_predict()` usan `posterior-prediction-2`, con
  `model`, `method`, `target`, intervalo y resumen comunes. Asignaciones y
  draws sólo aparecen cuando aplican o se solicitaron.
- La API de evaluación inicialmente usó un contenedor único. El incremento
  posterior lo retiró para que WAIC, MSE/R² y PPC no compartan estado ni
  aparenten medir la misma propiedad.
- No se añadieron alias para los esquemas retirados. Tests y consumidores de
  la tesis se migraron al contrato vigente y una prueba transversal exige
  igualdad de nombres entre los tres modelos.

## Duodécimo incremento: simulación y estado privado

- LRM, HLRM y CHLRM devuelven `simulation-result-1`: `model`, `train`,
  `truth` y `configuration`; `test` sólo aparece cuando el generador realmente
  lo produce. La verdad usa nombres semánticos y no `real_parameters` ni
  `cluster_info`.
- Se retiraron sin aliases `sim_lrm`, `sim_hlrm`, `sim_hlrm_hierarchical`,
  `sim_chlrm` y `sim_chlrm_hierarchical`. Sus consumidores activos usan
  `simulate_lrm()`, `simulate_hlrm()`, `simulate_hlrm_hierarchy()`,
  `simulate_chlrm_hierarchy()` o el generador del protocolo.
- Toda simulación con semilla restaura el RNG del llamador. El protocolo CHLRM
  permite asignaciones fijas para escenarios balanceados sin confundirlas con
  pesos poblacionales.
- El simulador HLRM condicional ya no codifica medias `30 + 5*j`, redondeo ni
  probabilidades binarias dependientes del número de columna: probabilidades,
  medias y desviaciones son entradas explícitas y quedan en la configuración.
- La notación matemática compacta permanece dentro de kernels para mantener
  correspondencia con las ecuaciones. La API pública usa nombres descriptivos;
  `rho` queda reservado para responsabilidades y `learning_rate` para el paso
  SVI.
- CHLRM Gibbs ya no mantiene `Betaj`, `sigma2j` ni `nk` como estado duplicado:
  se derivan de asignaciones y parámetros de componente al almacenar cada draw.
  Esto elimina posibles estados internamente incompatibles sin cambiar la
  distribución objetivo ni la salida de cadenas.

## Decimotercer incremento: contrato de log-likelihood

- `R/log-likelihood.R` define `log-likelihood-1` con valores pointwise,
  identificadores de unidad predictiva, modelo, condicionamiento y variables
  integradas. Una matriz sin metadatos ya no es la única fuente de verdad del
  estimando.
- LRM admite unidades observación o grupos suministrados; HLRM distingue
  parámetros locales condicionados de la marginalización de efectos y varianzas
  de grupo; CHLRM usa la likelihood marginal sobre una asignación compartida
  por grupo.
- GOF de los tres modelos y WAIC consumen el contrato y conservan su
  procedencia. DIC CHLRM/HLRM reutiliza la misma ruta marginal. Las likelihoods
  escalares condicionales del Gibbs no se eliminaron porque pertenecen al
  objetivo interno del sampler, no a la comparación predictiva.
- Las pruebas manuales comprueban que agregar observaciones reproduce la
  likelihood conjunta por grupo, que las rutas marginales coinciden con los
  cálculos especializados y que estimandos incompatibles producen error.

## Decimocuarto incremento: evaluación por responsabilidad

- Se retiró el contenedor monolítico de evaluación y el estado mutable asociado.
  `evaluate_in_sample_fit()` describe ajuste observado,
  `posterior_predictive_check()` evalúa adecuación generativa y
  `predictive_criteria()` calcula WAIC. Ninguna operación modifica el
  objeto ajustado.
- `log_likelihood()` permanece como infraestructura pointwise. El usuario elige
  un objetivo predictivo de alto nivel y el backend usa condicionamiento local
  para grupos HLRM observados, marginalización para grupos HLRM nuevos y para
  asignaciones CHLRM, y likelihood gaussiana ordinaria en LRM.
- En LRM agrupar filas sólo cambia la unidad predictiva y produce `waic_group`;
  no se denomina marginal ni se presenta como predicción de grupos nuevos.
- Los nueve motores implementados exponen nombres iguales por responsabilidad y
  contratos `in-sample-fit-1`, `posterior-predictive-check-2` y
  `predictive-criteria-2`. Tests y análisis activos de la tesis fueron migrados
  sin wrappers de compatibilidad.

## Decisión registrada: simulación, priors y escala

La revisión posterior distinguió tres operaciones que antes aparecían cercanas:

1. simulación con verdad fija, donde coeficientes, varianzas y pesos son
   entradas y se generan `X`, variables latentes y `y`;
2. simulación jerárquica, donde se fijan parámetros poblacionales y se muestrean
   parámetros locales/de componente antes de generar datos;
3. predictiva posterior, que parte de un ajuste y permanece en
   `posterior_predict()`.

Las dos primeras son capacidades generales del paquete y deben devolver el
mismo contrato versionado. La geometría trigonométrica que garantiza distancia
`d`, los escenarios S1--S7c y la prior exacta del review son composición de la
suite experimental. Por ello los nombres `simulate_chlrm_protocol()`,
`make_component_coefficients()` y `prior_chlrm_protocol()` quedaron marcados
para separación antes de congelar el namespace. Esta decisión no invalida sus
cálculos actuales; cambia su propiedad y nombre final.

Los simuladores y las priors deben usar la misma familia matemática, soporte,
unidades y convención `shape/rate`, pero no se conectan copiando valores. La
verdad simulada es evidencia para evaluar el ajuste, no información disponible
para construir su prior. Un único objeto prior validado debe alimentar MCMC,
VI y SVI. Las especificaciones explícitas, las priors calibradas desde OLS y la
inicialización OLS son responsabilidades distintas y no compartirán un rótulo
ambiguo como `NON_INFORMATIVE`.

También se ratificó que no habrá estandarización silenciosa ni argumentos
`standardize`/`standarize` dentro de los modelos. Una decisión posterior añadió
al paquete un contrato objetivo de preprocesamiento explícito y reversible; el
paper decide cuándo aplicarlo. Farms aprende sus estadísticas únicamente en
entrenamiento y transforma parámetros/predicciones de vuelta cuando se reporta
la escala original. Los experimentos simulados no se reestandarizan
empíricamente porque su escala es parte de la verdad que se pretende recuperar.

## Decisión registrada: escalado reusable y retiro de DIC

Estandarizar puede mejorar el condicionamiento numérico, hacer comparables las
escalas de los coeficientes y reducir la sensibilidad de VI/SVI a la tasa de
aprendizaje. También modifica interceptos, pendientes, varianzas residuales,
covarianzas jerárquicas, priors y densidades. Por ello el paquete ofrecerá un
objeto `data_scaling_spec` explícito, serializable y reversible, no una bandera
oculta dentro de `fit()`. El objeto se aprende sólo en entrenamiento, registra
columnas/escala y permite transformar datos, priors, parámetros y predicciones;
los log scores incluyen el jacobiano cuando se convierten de escala.

La decisión anterior de conservar DIC como diagnóstico secundario queda
sustituida. DIC se retirará de la API y del estudio final porque usa la
likelihood en un punto posterior, depende de la parametrización y, en CHLRM,
también de una representación relabelada de una posterior singular/multimodal.
El protocolo del review lo exige sólo si se conserva. La evidencia de
comparación será WAIC con estimando explícito, PSIS-LOO cuando sus diagnósticos
sean aceptables, CV/LOGO y MLPD de test. WAIC tampoco es universalmente
comparable: requiere los mismos datos, escala de respuesta, unidad predictiva y
condicionamiento/marginalización.

## Decisión registrada: cierre de evaluación predictiva reusable

P3 queda cerrado en el paquete. `predictive_performance()` separa las métricas
fuera de muestra —RMSE, MAE, R², cobertura y ancho de intervalo, ELPD y MLPD—
del ajuste descriptivo, PPC y WAIC. CV/LOGO reajusta el modelo y aprende el
escalado dentro de cada fold; PSIS-LOO se reserva para MCMC convergente y no se
aplica como si los draws de CAVI/SVI fueran draws exactos de la posterior.

El contrato `posterior-predictive-check-2` reemplaza el anterior sin alias. Su
preset estándar implementa T1--T5 cuando son matemáticamente aplicables; T5 usa
como referencia fija el OLS estimado en los datos observados. Reporta el MCSE
IID y, en MCMC, un MCSE basado en el ESS de la secuencia
$I\{T(y^{rep})\geq T(y)\}$ ordenada por cadena. Si hay menos de cuatro draws
seleccionados por cadena, la incertidumbre dependiente queda explícitamente no
estimable en vez de presentar el MCSE IID como corregido.

Los valores editoriales —5000 replicados, escenarios, folds y agregación entre
réplicas— siguen siendo responsabilidad de `paper/evaluation/`.

## Referencias metodológicas

- Spiegelhalter, Best, Carlin y van der Linde (2002),
  [Bayesian measures of model complexity and fit](https://doi.org/10.1111/1467-9868.00353).
- Watanabe (2010),
  [Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular Learning Theory](https://www.jmlr.org/papers/v11/watanabe10a.html).
- Watanabe (2021),
  [WAIC and WBIC for mixture models](https://doi.org/10.1007/s41237-021-00133-z).
- Vehtari, Gelman y Gabry (2017),
  [Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC](https://doi.org/10.1007/s11222-016-9696-4).
