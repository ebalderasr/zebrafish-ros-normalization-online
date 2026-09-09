# zebrafish-ros-normalization-online

*Versión en español. El [README en inglés](README.md) es el principal.*

Pipeline reproducible para normalizar y analizar mediciones de **fluorescencia
DCF** de especies reactivas de oxígeno en embriones de pez cebra. Viene en dos
formas que calculan lo mismo:

- **[Aplicación web](https://ebalderasr.github.io/zebrafish-ros-normalization-online/)**,
  que corre por completo en el navegador con Pyodide. Sin instalación, sin
  subir archivos, sin servidor. Los datos no salen de la máquina.
- **Paquete offline** (`zebrafish_ros`), una CLI y biblioteca de Python para
  análisis programable y versionado.

Un conjunto de pruebas automáticas verifica que ambos produzcan idénticos
log₂FC por embrión, anclas de control, marcas de outlier y tablas por sesión.
Ver [Paridad](#paridad-entre-las-dos-implementaciones).

El problema que resuelve es concreto. DCF reporta un solo canal, así que su
intensidad absoluta no trae corrección interna por carga de sonda, tamaño del
embrión ni iluminación. En los datos de ejemplo, la mediana del control por
sesión abarca 2.6 veces, mientras que los efectos de los fármacos son del orden
de 15 a 30 %. Sin corregir, la deriva entre sesiones supera al efecto que se
quiere medir.

```
CSV crudos  ->  tidy  ->  marcado de outliers  ->  normalización intra-sesión
            ->  estadística anidada  ->  tablas Prism  ->  figuras
```

> El texto de Methods listo para manuscrito está en **[METHODS.md](METHODS.md)**,
> en inglés, tipografiado en **[docs/methods.pdf](docs/methods.pdf)**.

### Repositorio hermano

Las mediciones ratiométricas de HyPer del mismo estudio se procesan en
[hyper-normalizer](https://github.com/ebalderasr/hyper-normalizer). Ambos
repositorios aceptan el mismo formato de entrada, aplican normalización
intra-sesión contra el control y corren la misma estadística con la sesión como
réplica, de modo que los dos ensayos se pueden reportar juntos. Las diferencias
propias de DCF son el ancla por mediana y la rama de outliers, explicadas
abajo.

---

## La matemática

### 1. El dato de entrada es una intensidad de un solo canal

HyPer es ratiométrico y divide por sí mismo la concentración de sonda. DCF no.
Lo único que se mide es una intensidad por embrión, así que el número crudo no
corrige cuánta sonda captó el embrión ni qué tan brillante estaba la lámpara
esa mañana. Toda comparación entre sesiones depende por tanto del ancla de la
sección 3.

La unidad de medición es un embrión. Las filas no están pareadas: dos valores
en la misma fila del CSV son dos embriones distintos medidos el mismo día.

### 2. Formato tidy y el control compartido

Los CSV crudos son matrices anchas: una fila por embrión, una columna por
condición, celdas vacías donde no hubo medición. Se convierten a formato largo
`(GRUPO, FECHA, TRATAMIENTO, INTENSIDAD)`.

Los nombres de archivo pueden seguir la convención `{GRUPO} DCF {PANEL}.csv`,
la misma que usa `hyper-normalizer`. Los archivos que comparten grupo se
analizan juntos y su control común se de-duplica: cuando un experimento se
reparte en varios paneles de fármacos, el control se mide una vez y se copia a
la hoja de cada panel, así que contarlo por archivo multiplicaría su n y
sesgaría el ancla. Un archivo cuyo nombre no siga la convención se vuelve su
propio grupo, que es como se ha comportado siempre la aplicación web.

El lector es deliberadamente permisivo: acepta comas decimales, encabezados con
acentos y mayúsculas inconsistentes, números guardados como texto y columnas en
blanco.

### 3. Marcado de outliers

Dentro de cada `(grupo, fecha, tratamiento)` se marcan los embriones fuera de
los límites de Tukey:

$$L_{\text{inf}} = Q_1 - 1.5\,\mathrm{IQR}, \qquad L_{\text{sup}} = Q_3 + 1.5\,\mathrm{IQR}$$

Los grupos con menos de cuatro embriones quedan sin marcar, porque sus
cuartiles se interpolan a partir de muy pocos puntos.

Se producen dos ramas. `with_outliers/` conserva todos los embriones;
`without_outliers/` elimina los marcados **y recalcula la normalización**.
Recalcular es el punto: quitar un outlier de un grupo control cambia el ancla
de esa sesión, así que filtrar después de normalizar dejaría a los embriones
restantes divididos por un denominador que ya no existe. Una prueba verifica
que al menos una sesión cambie de ancla entre ramas.

### 4. Normalización intra-sesión

El modelo implícito es multiplicativo:

$$I_i \;=\; \mu_t \cdot \beta_{g,d} \cdot \varepsilon_i$$

$\mu_t$ es el efecto de la condición, $\beta_{g,d}$ el factor de lote del grupo
$g$ en la fecha $d$, y $\varepsilon_i$ el ruido entre embriones. Un resumen del
control del mismo grupo y la misma fecha estima $\beta_{g,d}$, así que dividir
por él cancela el término:

$$\mathrm{ratio}_i \;=\; \frac{I_i}{\tilde{I}^{\,\mathrm{ctrl}}_{g,d}}, \qquad
\tilde{I}^{\,\mathrm{ctrl}}_{g,d} \;=\; \operatorname{mediana}\left(I_j : j \in \mathrm{ctrl}(g,d)\right)$$

La **mediana** es el ancla por omisión, que es la que usa la aplicación web y
la que conviene a intensidades de un canal con valores extremos ocasionales.
`--anchor mean` cambia a la media aritmética, igual que `hyper-normalizer`. La
elección mueve un poco los cambios de veces que se muestran y deja intacto el
test estadístico, cosa que una prueba verifica con precisión de 1 parte en
10¹².

Se siguen dos consecuencias exactas, y ambas condicionan la estadística:

1. El control normalizado queda centrado en 1 en cada sesión, por
   construcción. Ese valor es la línea de referencia de las figuras. Es también
   la razón por la que el control no puede tratarse como muestra libre en un
   test: perdió un grado de libertad por sesión.
2. Todos los embriones de una sesión comparten denominador, así que no son
   independientes entre sí.

La propiedad que justifica el método está cubierta por una prueba,
`test_normalization_cancels_the_batch_factor`: multiplicar por 7 todas las
mediciones de una sesión, que es lo que hace un cambio de ganancia, deja los
valores normalizados idénticos con precisión de 1 parte en 10¹².

Que funcione también se mide directamente. `variation.csv` reporta el CV entre
fechas de las medianas diarias, antes y después de normalizar. En los datos de
ejemplo baja de 0.29–0.45 a 0.04–0.15, para todas las condiciones.

### 5. Estadística: la sesión es la réplica

Un test estándar sobre los embriones agrupados choca con las dos consecuencias
anteriores. La corrección es trabajar en escala logarítmica y usar la sesión de
adquisición como unidad de réplica:

$$\delta_d \;=\; \overline{\log_2 I}\big|_{t,d} \;-\; \overline{\log_2 I}\big|_{\mathrm{ctrl},d}$$

La diferencia se toma dentro de la sesión y sobre las intensidades crudas, así
que $\log_2\beta_{g,d}$ se cancela algebraicamente. El test no depende entonces
ni de la normalización ni de la elección de ancla. Los $\delta_d$, uno por
sesión, se contrastan contra 0 con una prueba t de una muestra de dos colas,
con $n$ igual al número de sesiones. El efecto reportado es $2^{\bar\delta}$,
el cambio geométrico de veces, con su IC del 95 %. Los p-valores se corrigen
por Holm-Bonferroni dentro de cada grupo.

Cada sesión pesa igual, sin importar cuántos embriones contenga.
`--min-embryos 3` descarta las sesiones por debajo de ese número y conviene
correrlo como prueba de sensibilidad.

Con `--mixed-model` se ajusta además la versión completa,

$$\log_2 I_i \;=\; \alpha_t + b_d + \varepsilon_i, \qquad b_d \sim \mathcal{N}(0, \tau^2)$$

que estima el efecto de sesión y el de condición a la vez, en lugar de
normalizar primero y contrastar después.

### 6. Cuánto importa

Sobre los datos sintéticos de ejemplo, con outliers eliminados:

| Grupo | Cond. | sesiones | fold | IC 95 % | *p* Holm | *p* ingenuo |
|-------|-------|---------:|-----:|---------|---------:|------------:|
| WT | NAC | 5 | 0.719 | [0.672, 0.770] | 0.0011 | 2.6×10⁻¹⁷ |
| WT | APO | 5 | 0.840 | [0.799, 0.883] | 0.0032 | 9.1×10⁻⁸ |
| MUT | EUK | 5 | 0.920 | [0.812, 1.043] | 0.1376 | 6.1×10⁻² |

La columna `p_naive_pooled` trata cada embrión como réplica independiente.
Queda hasta trece órdenes de magnitud por debajo de `p_holm`. El pipeline la
calcula para cuantificar esa inflación; el valor que se reporta es `p_holm`.

### 7. Tablas de Prism

Dos formas, y elegir entre ellas es una decisión metodológica, no de formato.

`prism_by_date.csv` da una fila por sesión de adquisición y una columna por
condición, con la mediana del log₂FC de esa sesión en cada celda. N es el
número de sesiones, que son las réplicas independientes del diseño. **Esta es
la tabla que se pega en Prism.**

`prism_by_embryo.csv` da una fila por embrión con todas las sesiones juntas. N
pasa a ser el número de embriones, típicamente de 4 a 10 veces mayor.
Alimentarla a una prueba t o a un ANOVA sin modelar el anidamiento infla la
significancia. Se exporta para inspección y para diseños que tomen en cuenta
esa estructura de forma explícita.

### 8. Figuras

La figura principal es un **SuperPlot** ([Lord et al., *J. Cell Biol.* 2020](https://doi.org/10.1083/jcb.202001064)).
Los embriones individuales van como puntos grises tenues de fondo, con la media
de cada sesión encima como marcador grande de color. Así se ve cuántas sesiones
sostienen cada caja, y si un efecto se mantiene entre sesiones o proviene de un
solo día. Un boxplot sobre embriones agrupados no muestra ninguna de las dos
cosas.

Una figura complementaria grafica el ancla de control cruda de cada sesión, que
es la deriva que la normalización elimina.

---

## Uso offline

```bash
git clone https://github.com/ebalderasr/zebrafish-ros-normalization-online.git
cd zebrafish-ros-normalization-online
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

```bash
# Datos sintéticos incluidos
python -m zebrafish_ros --input-dir data/example --output-dir results

# Tus datos, con las convenciones de hyper-normalizer
python -m zebrafish_ros \
    --input-dir data/raw \
    --output-dir results \
    --anchor mean --min-embryos 3 --mixed-model
```

Opciones principales:

| Opción | Efecto |
|---|---|
| `--control NOMBRE` | Columna que hace de control (por omisión `DMSO`). |
| `--anchor median\|mean` | Estadístico que resume el control de cada sesión. `median` iguala a la app web, `mean` a `hyper-normalizer`. El test no cambia. |
| `--outliers keep\|drop\|both` | Qué ramas escribir (por omisión `both`). |
| `--min-embryos N` | Descarta del test las sesiones con menos de N embriones. |
| `--strict` | Falla si algún `(grupo, fecha)` no tiene control, en vez de descartarlo. |
| `--mixed-model` | Ajusta `log2(I) ~ tratamiento + (1 \| fecha)`. Requiere `statsmodels`. |
| `--no-plots` | Omite las figuras. |

Como biblioteca:

```python
from pathlib import Path
from zebrafish_ros import run

result = run(Path("data/example"), Path("results"))
for test in result.branches["drop"].tests:
    print(test.group, test.treatment, round(test.geo_fold, 3), round(test.p_holm, 4))
```

## Uso en el navegador

Abre la [aplicación](https://ebalderasr.github.io/zebrafish-ros-normalization-online/),
suelta uno o más CSV, elige la columna de control de cada archivo y analiza. La
app dibuja gráficas interactivas y exporta todas las tablas en un ZIP.

Para servirla localmente:

```bash
python3 -m http.server 8000   # luego abre http://localhost:8000
```

La aplicación web descarga Pyodide, NumPy y Plotly de CDN públicos, así que
necesita conexión en la primera carga. Para trabajar sin red, usa el paquete
offline.

## Formato de entrada

Un archivo por combinación de grupo y panel, nombrado `{GRUPO} DCF {PANEL}.csv`,
o cualquier nombre de CSV para que el archivo se trate como su propio grupo:

```csv
FECHA,DMSO,VAS,DPI,APO
240115,108.17,79.12,91.29,100.31
240115,89.69,117.49,88.20,117.20
240115,97.87,78.72,,75.07
```

Primera columna la fecha de adquisición (YYMMDD), el resto condiciones, una
fila por embrión, celdas vacías donde no hubo medición. Es el mismo formato que
acepta `hyper-normalizer`, así que un experimento se exporta una vez y lo
procesa cualquiera de los dos pipelines. Cuando un grupo se reparte en varios
paneles, repite la columna del control en cada uno; el pipeline la de-duplica.

## Salidas

En el nivel superior:

| Archivo | Contenido |
|---|---|
| `tidy.csv` | Una fila por embrión, con el control ya de-duplicado. |
| `outliers_flagged.csv` | Cada embrión marcado, con los límites que lo excluyeron. |

Y luego una vez por rama, bajo `with_outliers/` y `without_outliers/`:

| Archivo | Contenido |
|---|---|
| `normalized.csv` | Añade `ratio_norm`, `log2_norm` y el ancla usada. |
| `control_anchors.csv` | Resumen del control por sesión: n, media, mediana, DE, estado. |
| `summary.csv` | Descriptivos por grupo y condición, aritméticos y geométricos. |
| `date_folds.csv` | Una fila por `(grupo, fecha, tratamiento)`, las réplicas reales. |
| `tests.csv` | Contraste vs control con la sesión como réplica, IC y p corregido. |
| `variation.csv` | CV entre fechas de las medianas diarias, antes y después. |
| `prism_by_date.csv` | Tabla ancha, una fila por sesión. Esta es la de inferencia. |
| `prism_by_embryo.csv` | Tabla ancha, una fila por embrión. Lee la sección 7 antes. |
| `figure_groups.png` | SuperPlot por grupo experimental. |
| `figure_control_drift.png` | Ancla de control cruda por sesión. |

Los flotantes se escriben con 12 cifras significativas, de modo que releer un
archivo de salida para un análisis posterior no vuelve a redondear.

## Paridad entre las dos implementaciones

`tests/test_engine_parity.py` corre `zebrafish_ros_engine.analyze_one_file`,
que es el módulo que Pyodide carga en el navegador, contra el paquete offline
sobre los mismos archivos, y verifica que coincidan los log₂FC por embrión, las
anclas de control, las marcas de outlier y las celdas de la tabla Prism por
sesión.

Hay una diferencia documentada en lugar de corregida. El motor del navegador
ordena sus filas largas por `(fecha, número de fila, condición)` antes de armar
la tabla Prism, así que sus columnas salen alfabéticas. El paquete offline
conserva el orden que las condiciones tienen en el encabezado del CSV, con el
control primero. Todos los valores son idénticos; solo cambia el orden de las
columnas.

## Datos

Este repositorio no incluye los datos experimentales reales. `data/example/`
contiene un conjunto sintético generado por `scripts/make_example_data.py`, con
la misma estructura que los datos del laboratorio: dos grupos, dos paneles,
control compartido, sesiones desbalanceadas, celdas faltantes, deriva fuerte
entre sesiones y embriones extremos plantados.

Para regenerarlo:

```bash
python scripts/make_example_data.py --seed 20260908
```

Los CSV reales van en `data/raw/`, que está en `.gitignore`.

## Pruebas

```bash
pip install -r requirements-dev.txt
pip install -e .
pytest -q
```

Las pruebas cubren la paridad con el motor del navegador y las invariantes
matemáticas: que el control normalizado quede centrado en 1, que escalar una
sesión entera deje sin cambio todos los valores normalizados, que la elección
de ancla no mueva ni un p-valor, que la normalización baje el CV entre fechas
de todas las condiciones, y que la *n* de los tests cuente sesiones y no
embriones.

## Estructura

```
index.html, style.css, app.js       aplicación web
zebrafish_ros_engine.py             motor de análisis que carga Pyodide
zebrafish_ros/                      paquete offline
    tidy.py        lectura de CSV anchos y de-duplicación del control
    outliers.py    marcado de Tukey dentro de fecha y condición
    normalize.py   normalización intra-sesión
    stats.py       descriptivos, sesión como réplica, Holm, modelo mixto
    prism.py       tablas anchas para GraphPad Prism
    plots.py       SuperPlots y la figura de deriva del control
    pipeline.py    orquestación y escritura de salidas
    __main__.py    CLI
scripts/
    make_example_data.py
    build_methods_pdf.py
data/example/                       datos sintéticos versionados
tests/
docs/methods.pdf                    versión tipografiada de METHODS.md
METHODS.md                          texto de Methods para manuscrito
```

## Cita

Si este código contribuye a un trabajo publicado, cita el repositorio y el
método de SuperPlot:

> Lord SJ, Velle KB, Mullins RD, Fritz-Laylin LK. SuperPlots: Communicating
> reproducibility and variability in cell biology. *J Cell Biol.*
> 2020;219(6):e202001064. doi:10.1083/jcb.202001064

## Autor

**Emiliano Balderas Ramírez**
Bioingeniero, candidato a Doctor en Ciencias Bioquímicas
Instituto de Biotecnología (IBt), UNAM

## Licencia

MIT. Ver [LICENSE](LICENSE).
