# zebrafish-ros-normalization-online

[![tests](https://github.com/ebalderasr/zebrafish-ros-normalization-online/actions/workflows/tests.yml/badge.svg)](https://github.com/ebalderasr/zebrafish-ros-normalization-online/actions/workflows/tests.yml)

Reproducible pipeline for normalizing and analysing **DCF fluorescence**
measurements of reactive oxygen species in zebrafish embryos. It ships in two
forms that compute the same thing:

- **[Browser app](https://ebalderasr.github.io/zebrafish-ros-normalization-online/)**,
  which runs entirely client-side through Pyodide. No installation, no upload,
  no backend. Data never leaves the machine.
- **Offline package** (`zebrafish_ros`), a Python CLI and library for scripted,
  version-controlled analysis.

An automated test suite checks that the two produce identical per-embryo
log₂FC, control anchors, outlier flags and per-session tables. See
[Parity](#parity-between-the-two-implementations).

The problem it solves is specific. DCF reports a single channel, so its
absolute intensity carries no internal correction for probe loading, embryo
size or illumination. In the bundled example data the session median of the
vehicle control spans 2.6-fold, while the drug effects are on the order of 15
to 30 %. Without correction, the between-session drift is larger than the
effect under study.

```
raw CSVs  ->  tidy  ->  outlier flagging  ->  within-session normalization
          ->  nested statistics  ->  Prism tables  ->  figures
```

> Manuscript-ready Methods text is in **[METHODS.md](METHODS.md)**, typeset to
> **[docs/methods.pdf](docs/methods.pdf)**.
> *Este README también está disponible [en español](README.es.md).*

### Companion repository

The ratiometric HyPer measurements of the same study are processed by
[hyper-normalizer](https://github.com/ebalderasr/hyper-normalizer). Both
repositories accept the same input layout, apply within-session normalization
against the control, and run the same session-as-replicate statistics, so the
two assays can be reported side by side. The differences specific to DCF are
the median anchor and the outlier branch, both explained below.

---

## The mathematics

### 1. The input is a single-channel intensity

HyPer is ratiometric and divides out probe concentration on its own. DCF does
not. One intensity per embryo is all that is measured, so nothing in the raw
number corrects for how much probe the embryo took up or how bright the lamp
was that morning. Every comparison across sessions therefore depends on the
anchor described in section 3.

The unit of measurement is one embryo. Rows are not paired: two values on the
same CSV row are two different embryos imaged on the same date.

### 2. Tidy format and the shared control

The raw CSVs are wide matrices: one row per embryo, one column per condition,
empty cells where no measurement was taken. These are converted to long format,
`(GROUP, DATE, TREATMENT, INTENSITY)`.

File names may follow the `{GROUP} DCF {PANEL}.csv` convention, the same shape
`hyper-normalizer` uses. Files sharing a group are analysed together and their
common control is de-duplicated: when one experiment is split across drug
panels, the vehicle control is measured once and copied into each panel's
sheet, so counting it per file would multiply its n and bias the anchor. A file
whose name does not match the convention becomes its own group, which is how
the browser app has always behaved.

The parser is deliberately permissive and accepts decimal commas, accented and
inconsistently capitalised headers, numbers stored as text, and blank columns.

### 3. Outlier flagging

Within each `(group, date, treatment)`, embryos outside Tukey's bounds are
flagged:

$$L_{\text{low}} = Q_1 - 1.5\,\mathrm{IQR}, \qquad L_{\text{high}} = Q_3 + 1.5\,\mathrm{IQR}$$

Groups of fewer than four embryos are left unflagged, since their quartiles are
interpolated from too few points.

Two branches are produced. `with_outliers/` keeps every embryo;
`without_outliers/` removes the flagged ones **and recomputes the
normalization**. Recomputing is the point: dropping an outlier from a control
group changes that session's anchor, so filtering after normalizing would leave
the remaining embryos divided by a denominator that no longer exists. A test
checks that at least one session's anchor actually moves between branches.

### 4. Within-session normalization

The implicit model is multiplicative:

$$I_i \;=\; \mu_t \cdot \beta_{g,d} \cdot \varepsilon_i$$

$\mu_t$ is the condition effect, $\beta_{g,d}$ the batch factor of group $g$ on
date $d$, and $\varepsilon_i$ the between-embryo noise. A summary of the
control from the same group and date estimates $\beta_{g,d}$, so dividing by it
cancels the term:

$$\mathrm{ratio}_i \;=\; \frac{I_i}{\tilde{I}^{\,\mathrm{ctrl}}_{g,d}}, \qquad
\tilde{I}^{\,\mathrm{ctrl}}_{g,d} \;=\; \operatorname{median}\left(I_j : j \in \mathrm{ctrl}(g,d)\right)$$

The **median** is the default anchor, which is what the browser app uses and
what suits single-channel intensities with occasional extreme values.
`--anchor mean` switches to the arithmetic mean, matching `hyper-normalizer`.
The choice moves the displayed fold changes slightly and leaves the statistical
test untouched, which a test asserts to 1 part in 10¹².

Two exact consequences follow, and both constrain the statistics:

1. The normalized control is centred on 1 in every session, by construction.
   That value is the reference line in the figures. It is also why the control
   group cannot be treated as a freely varying sample in a test: it lost a
   degree of freedom per session.
2. Every embryo from a given session shares the same denominator, so embryos
   within a session are not independent of one another.

The property that justifies the method is covered by a test,
`test_normalization_cancels_the_batch_factor`: multiplying every measurement of
one session by 7, which is what a change in detector gain does, leaves the
normalized values identical to 1 part in 10¹².

That it works is also measured directly. `variation.csv` reports the across-date
CV of the daily medians before and after normalization. On the example data it
falls from 0.29–0.45 to 0.04–0.15, for every condition.

### 5. Statistics: the session is the replicate

A standard test over pooled embryos runs into both consequences above. The
correction is to work on a logarithmic scale and use the acquisition session as
the unit of replication:

$$\delta_d \;=\; \overline{\log_2 I}\big|_{t,d} \;-\; \overline{\log_2 I}\big|_{\mathrm{ctrl},d}$$

The difference is taken within a session and on the raw intensities, so
$\log_2\beta_{g,d}$ cancels algebraically. The test therefore depends neither
on the normalization nor on the anchor choice. The $\delta_d$ values, one per
session, are tested against zero with a two-tailed one-sample *t* test, with
$n$ equal to the number of sessions. The reported effect is $2^{\bar\delta}$,
the geometric fold change, with its 95 % CI. *P* values are Holm-Bonferroni
corrected within each group.

Every session carries equal weight, regardless of how many embryos it holds.
`--min-embryos 3` drops sessions below that count and is worth running as a
sensitivity check.

With `--mixed-model`, the full version is also fitted,

$$\log_2 I_i \;=\; \alpha_t + b_d + \varepsilon_i, \qquad b_d \sim \mathcal{N}(0, \tau^2)$$

which estimates session and condition effects jointly instead of normalizing
first and testing afterwards.

### 6. How much this matters

On the synthetic example data, outliers removed:

| Group | Cond. | sessions | fold | 95 % CI | *p* Holm | *p* naive |
|-------|-------|---------:|-----:|---------|---------:|----------:|
| WT | NAC | 5 | 0.719 | [0.672, 0.770] | 0.0011 | 2.6×10⁻¹⁷ |
| WT | APO | 5 | 0.840 | [0.799, 0.883] | 0.0032 | 9.1×10⁻⁸ |
| MUT | EUK | 5 | 0.920 | [0.812, 1.043] | 0.1376 | 6.1×10⁻² |

The `p_naive_pooled` column treats every embryo as an independent replicate. It
sits up to thirteen orders of magnitude below `p_holm`. The pipeline computes
and prints it to quantify that inflation; the value to report is `p_holm`.

### 7. Prism tables

Two shapes, and the choice between them is methodological rather than
cosmetic.

`prism_by_date.csv` gives one row per acquisition session and one column per
condition, each cell the median log₂FC of that session. N is the number of
sessions, which is the design's number of independent replicates. **This is the
table to paste into Prism.**

`prism_by_embryo.csv` gives one row per embryo with all sessions pooled. N is
then the number of embryos, typically 4 to 10 times larger. Feeding it to a
*t* test or an ANOVA without modelling the nesting inflates significance. It is
exported for inspection and for designs that account for the structure
explicitly.

### 8. Figures

The main figure is a **SuperPlot** ([Lord et al., *J. Cell Biol.* 2020](https://doi.org/10.1083/jcb.202001064)).
Individual embryos appear as faint grey background points, with the mean of
each acquisition session drawn on top as a large colour-coded marker. This
shows how many sessions support each box, and whether an effect holds across
sessions or comes from a single date. A boxplot over pooled embryos shows
neither.

A companion figure plots the raw control anchor per session, which is the drift
the normalization removes.

---

## Offline use

```bash
git clone https://github.com/ebalderasr/zebrafish-ros-normalization-online.git
cd zebrafish-ros-normalization-online
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

```bash
# Bundled synthetic data
python -m zebrafish_ros --input-dir data/example --output-dir results

# Your data, matching hyper-normalizer's conventions
python -m zebrafish_ros \
    --input-dir data/raw \
    --output-dir results \
    --anchor mean --min-embryos 3 --mixed-model
```

Main options:

| Option | Effect |
|---|---|
| `--control NAME` | Column acting as the control (default `DMSO`). |
| `--anchor median\|mean` | Statistic summarising each session's control. `median` matches the browser app, `mean` matches `hyper-normalizer`. The test is unaffected. |
| `--outliers keep\|drop\|both` | Which branches to write (default `both`). |
| `--min-embryos N` | Drop sessions measured in fewer than N embryos from the test. |
| `--strict` | Fail if any `(group, date)` lacks a control, instead of dropping it. |
| `--mixed-model` | Fit `log2(I) ~ treatment + (1 \| date)`. Requires `statsmodels`. |
| `--no-plots` | Skip figure generation. |

As a library:

```python
from pathlib import Path
from zebrafish_ros import run

result = run(Path("data/example"), Path("results"))
for test in result.branches["drop"].tests:
    print(test.group, test.treatment, round(test.geo_fold, 3), round(test.p_holm, 4))
```

## Browser use

Open the [live app](https://ebalderasr.github.io/zebrafish-ros-normalization-online/),
drop one or more CSVs, choose the control column for each file, and analyse.
The app renders interactive plots and exports every table as a ZIP.

To serve it locally:

```bash
python3 -m http.server 8000   # then open http://localhost:8000
```

The browser app fetches Pyodide, NumPy and Plotly from public CDNs, so it needs
a network connection on first load. For work without network access, use the
offline package.

## Input format

One file per group-and-panel combination, named `{GROUP} DCF {PANEL}.csv`, or
any CSV name to have the file treated as its own group:

```csv
FECHA,DMSO,VAS,DPI,APO
240115,108.17,79.12,91.29,100.31
240115,89.69,117.49,88.20,117.20
240115,97.87,78.72,,75.07
```

First column the acquisition date (YYMMDD), the rest conditions, one row per
embryo, empty cells where no measurement was taken. This is the same layout
`hyper-normalizer` accepts, so one experiment can be exported once and
processed by either pipeline. When a group is split across panels, repeat the
control column in each; the pipeline de-duplicates it.

## Outputs

At the top level:

| File | Contents |
|---|---|
| `tidy.csv` | One row per embryo, control already de-duplicated. |
| `outliers_flagged.csv` | Every flagged embryo with the bounds that excluded it. |

Then once per branch, under `with_outliers/` and `without_outliers/`:

| File | Contents |
|---|---|
| `normalized.csv` | Adds `ratio_norm`, `log2_norm` and the anchor used. |
| `control_anchors.csv` | Control summary per session: n, mean, median, SD, status. |
| `summary.csv` | Descriptives per group and condition, arithmetic and geometric. |
| `date_folds.csv` | One row per `(group, date, treatment)`, the real replicates. |
| `tests.csv` | Contrast vs control with session as replicate, CI, corrected *p*. |
| `variation.csv` | Across-date CV of the daily medians, before and after normalization. |
| `prism_by_date.csv` | Wide table, one row per session. Use this one for inference. |
| `prism_by_embryo.csv` | Wide table, one row per embryo. Read section 7 first. |
| `figure_groups.png` | SuperPlot per experimental group. |
| `figure_control_drift.png` | Raw control anchor per session. |

Floats are written with 12 significant digits, so re-reading an output file for
a later analysis does not round a second time.

## Parity between the two implementations

`tests/test_engine_parity.py` runs `zebrafish_ros_engine.analyze_one_file`,
which is the module Pyodide loads in the browser, against the offline package
on the same files, and asserts that the per-embryo log₂FC, the control anchors,
the outlier flags and the per-session Prism cells all match.

One difference is documented rather than fixed. The browser engine sorts its
long rows by `(date, row number, condition)` before building the Prism table,
so its columns come out alphabetical. The offline package keeps the order the
conditions have in the CSV header, with the control first. Every cell value is
identical; only the columns sit in a different sequence.

## Data

This repository does not include the real experimental data. `data/example/`
holds a synthetic dataset generated by `scripts/make_example_data.py`, with the
same structure as the laboratory data: two groups, two panels, a shared
control, unbalanced sessions, missing cells, strong session drift and planted
extreme embryos.

To regenerate it:

```bash
python scripts/make_example_data.py --seed 20260908
```

Real CSVs belong in `data/raw/`, which is in `.gitignore`.

## Tests

```bash
pip install -r requirements-dev.txt
pip install -e .
pytest -q
```

The tests cover parity with the browser engine and the mathematical invariants:
that the normalized control is centred on 1, that rescaling a whole session
leaves every normalized value unchanged, that the anchor choice does not move a
single *p* value, that normalization lowers the across-date CV of every
condition, and that the *n* of the tests counts sessions and not embryos.

## Layout

```
index.html, style.css, app.js       browser app
zebrafish_ros_engine.py             analysis engine loaded by Pyodide
zebrafish_ros/                      offline package
    tidy.py        reading wide CSVs and de-duplicating the control
    outliers.py    Tukey flagging within date and condition
    normalize.py   within-session normalization
    stats.py       descriptives, session as replicate, Holm, mixed model
    prism.py       wide tables for GraphPad Prism
    plots.py       SuperPlots and the control drift figure
    pipeline.py    orchestration and output writing
    __main__.py    CLI
scripts/
    make_example_data.py
    build_methods_pdf.py
data/example/                       versioned synthetic data
tests/
docs/methods.pdf                    typeset version of METHODS.md
METHODS.md                          manuscript-ready Methods text
```

## Citation

If this code contributes to published work, please cite the repository and the
SuperPlot method:

> Lord SJ, Velle KB, Mullins RD, Fritz-Laylin LK. SuperPlots: Communicating
> reproducibility and variability in cell biology. *J Cell Biol.*
> 2020;219(6):e202001064. doi:10.1083/jcb.202001064

## Author

**Emiliano Balderas Ramírez**
Bioengineer, PhD candidate in Biochemical Sciences
Instituto de Biotecnología (IBt), UNAM

## Licence

MIT. See [LICENSE](LICENSE).
