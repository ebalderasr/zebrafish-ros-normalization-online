<div align="center">

# Zebrafish ROS Normalizer Online

### From raw DCF fluorescence CSVs to normalized, date-aware ROS analysis — entirely in your browser

<br>

**[→ Open the live app](https://ebalderasr.github.io/zebrafish-ros-normalization-online/)**

<br>

[![Stack](https://img.shields.io/badge/Stack-Pyodide_·_NumPy_·_Plotly-4A90D9?style=for-the-badge)]()
[![Focus](https://img.shields.io/badge/Focus-Zebrafish_Embryos_·_DCF_Analysis-34C759?style=for-the-badge)]()
[![License](https://img.shields.io/badge/License-MIT-blue?style=for-the-badge)](./LICENSE)
[![Part of](https://img.shields.io/badge/Part_of-HostCell_Lab_Suite-5856D6?style=for-the-badge)](https://github.com/ebalderasr)

</div>

---

## What is Zebrafish ROS Normalizer Online?

Zebrafish ROS Normalizer Online is a **browser-based analysis app** for DCF fluorescence experiments in zebrafish embryos. It takes one or more raw CSV files, interprets each non-empty cell as one independent embryo, normalizes intensities by acquisition date using an internal control from the same day, flags outliers, renders interactive plots, and packages the processed tables into a downloadable ZIP.

It runs entirely in the browser through [Pyodide](https://pyodide.org), so **no Python installation, no local scripts, and no data upload are required**. Your experimental data stays on your machine.

---

## Why it matters

DCF fluorescence is highly sensitive to acquisition settings and handling conditions. Across different experiment days, changes in laser power, detector gain, embryo timing, and operator workflow can shift raw intensities enough to make direct cross-date comparisons misleading.

Without a normalization strategy, a treatment can appear different simply because it was acquired under different microscope settings.

This app solves that by anchoring each date to an internal control measured on the same day. What used to require repeated spreadsheet cleanup, manual per-day normalization, and multiple plotting rounds can now be done in one browser session.

---

## How it works

Four steps. No command line. No backend.

<br>

**Step 1 — Upload one or more CSV files**

Drag and drop your files or click to browse. The app accepts multiple `.csv` files at once and ignores non-CSV inputs automatically. Each file is treated as a separate experiment and is analyzed independently.

<br>

**Step 2 — The app detects your columns**

The engine inspects each CSV to identify:

- the date column;
- the condition columns;
- and a default candidate control.

If a `DMSO`-like column exists, it is preselected automatically.

<br>

**Step 3 — Choose one control per file**

When you click **Analyze files**, the app opens a modal asking for the control column of **each uploaded CSV**. This step is mandatory. `DMSO` is selected by default when present, but the user can choose any other condition to serve as the internal reference.

<br>

**Step 4 — Analyze and explore**

The app processes every file in the browser, generates normalized tables and warnings, and displays interactive plots for:

- raw distribution by condition;
- raw median intensity across dates;
- control-anchor drift across dates;
- normalized distributions;
- and variation reduction after normalization.

Each file also includes downloadable per-file outputs, and the full session can be exported as a ZIP.

---

## Methods

### Input model

Each CSV is expected to contain:

- one date column;
- one column per experimental condition;
- zero or more empty cells.

Each non-empty numeric cell is treated as **one independent embryo**. If the same row contains values under multiple conditions, they are **not** interpreted as matched embryos; they are simply embryos acquired on the same date under different conditions.

---

### Normalization by date

For each file and each acquisition date \( d \), the user-selected control column is used as the internal anchor for that same day:

$$
\widetilde{I}_{\mathrm{control},d} = \mathrm{median}\left(I_{\mathrm{control},d}\right)
$$

Then each embryo intensity is normalized as:

$$
\mathrm{ratio}_{i,d} = \frac{I_{i,d}}{\widetilde{I}_{\mathrm{control},d}}
$$

and expressed on a log base 2 scale as:

$$
\log_2 \mathrm{FC}_{i,d} = \log_2 \left(\mathrm{ratio}_{i,d}\right)
$$

This design keeps normalization local to each file and each date. The engine does **not** use:

- a global reference across all dates;
- a pooled median across conditions;
- or a shared anchor across uploaded files.

---

### Outlier handling

Outliers are detected within each `date × condition` group using the standard 1.5×IQR rule:

$$
\mathrm{IQR} = Q_3 - Q_1
$$

$$
L_{\mathrm{inf}} = Q_1 - 1.5 \cdot \mathrm{IQR}
$$

$$
L_{\mathrm{sup}} = Q_3 + 1.5 \cdot \mathrm{IQR}
$$

Any embryo with intensity outside those bounds is flagged as a potential outlier.

The app returns two analysis branches for every file:

- `with_outliers`: all embryos retained;
- `without_outliers`: flagged outliers removed and normalization recalculated from the cleaned data.

This is important because removing outliers can also change the control median of a given date, which means the normalized values must be recomputed rather than simply filtered.

---

### Statistical unit and pseudo-replication

The normalization produces a log₂FC value for every individual embryo. The app exports two Prism table formats and **the researcher decides which one fits their experimental design**.

**Per-date table (recommended for standard statistical inference)**

Each acquisition date represents one independent experimental run. Embryos measured on the same day under the same condition are biological sub-replicates within that run, not independent experiments. In this design, the correct unit for statistical testing is the date, and N equals the number of dates with valid data:

| Date | DMSO | Treatment A | Treatment B |
|---|---|---|---|
| 2026-03-17 | 0.00 | 1.14 | −0.43 |
| 2026-03-18 | 0.00 | 0.98 | −0.61 |
| 2026-03-19 | 0.00 | 1.02 | −0.55 |

Each cell is the **median log₂FC** of all embryos in that date-condition group.

**Per-embryo table (use with care)**

Each row is one individual embryo's log₂FC, with all dates pooled. Columns are not paired. This format gives N = number of embryos, which can be 4–10× larger than the number of independent experimental days. Using it directly in t-tests or ANOVAs without accounting for the nested structure of the data inflates statistical power and increases the risk of false positives. It may be appropriate under specific designs (e.g., if the researcher treats embryos as the unit of inference within a single date, or applies a mixed model), but this should be a deliberate methodological choice.

---

### Warnings and traceability

The engine does not fail silently. It records warnings and notes for issues such as:

- missing or ambiguous date columns;
- non-numeric cells;
- empty condition columns;
- inconsistent or nearly duplicated condition names;
- dates without valid control measurements;
- low control sample size for a given date;
- and rows discarded during parsing.

All warnings are shown in the interface and are also exported as CSV.

---

## Features

| | |
|---|---|
| **Zero installation** | Runs fully client-side via Pyodide |
| **Multiple CSV support** | Each uploaded file is treated as an independent experiment |
| **Mandatory control selection** | One control column must be chosen for each file before analysis |
| **DMSO by default** | If a `DMSO` column is detected, it is preselected automatically |
| **Date-aware normalization** | Every date is normalized against the control measured on that same date |
| **Prism table by date** | Wide table with one row per experimental day — recommended unit for statistical inference |
| **Prism table by embryo** | Wide table with one row per embryo — use with awareness of pseudo-replication risk |
| **Dual output branches** | `sin_outliers` and `con_outliers`, both fully exported |
| **Interactive plots** | Plotly-based browser plots for raw, normalized and control-anchor views |
| **Explicit warnings** | Parsing, normalization, raw-data assumption, and statistical notes remain visible |
| **ZIP export** | All per-file outputs can be downloaded in one click |
| **No data upload** | Everything runs locally in the browser |

---

## Input format

### Required structure

Each CSV must contain:

| Element | Requirement |
|---|---|
| Date column | A date-like column detectable as `YYMMDD` rows, for example `260315` |
| Condition columns | One or more treatment/control columns |
| Values | Numeric fluorescence intensities, optionally using decimal commas |

### Supported quirks

The parser is intentionally permissive. It can handle:

- leading/trailing spaces in headers;
- accents and special characters in column names;
- inconsistent capitalization;
- empty cells;
- completely empty columns;
- numeric values stored as text;
- decimal commas such as `123,4`;
- and sparse condition coverage across dates.

### Important interpretation rule

Rows are **not** treated as paired embryos across conditions.

If a row contains:

| Fecha | DMSO | Treatment A |
|---|---:|---:|
| 260315 | 120 | 145 |

that means:

- one embryo under `DMSO` acquired on `2026-03-15`;
- one different embryo under `Treatment A` acquired on `2026-03-15`.

---

## Outputs

For each uploaded CSV, the app generates 14 files organized in two groups.

**Prism GraphPad tables** (primary deliverables):

| File | Contents |
|---|---|
| `*_PRISM_por_fecha_sin_outliers.csv` | Wide table: one row per date, one column per condition, value = median log₂FC. **Recommended for statistical inference.** |
| `*_PRISM_por_fecha_con_outliers.csv` | Same, retaining flagged outliers. |
| `*_PRISM_por_embrion_sin_outliers.csv` | Wide table: one column per condition, one row per embryo (all dates pooled), value = log₂FC. |
| `*_PRISM_por_embrion_con_outliers.csv` | Same, retaining flagged outliers. |

**Detailed tables** (for traceability and deeper inspection):

| File | Contents |
|---|---|
| `*_datos_normalizados_sin_outliers.csv` | Long-format table of all embryo measurements after normalization (outliers removed). |
| `*_datos_normalizados_con_outliers.csv` | Same, retaining flagged outliers. |
| `*_outliers_eliminados.csv` | Records of embryos flagged and removed as outliers. |
| `*_resumen_por_fecha_condicion_sin_outliers.csv` | Summary statistics (n, mean, median, SD) per date × condition. |
| `*_resumen_por_fecha_condicion_con_outliers.csv` | Same, retaining flagged outliers. |
| `*_ancla_control_por_fecha_sin_outliers.csv` | Control anchor statistics used for normalization, per date. |
| `*_ancla_control_por_fecha_con_outliers.csv` | Same, retaining flagged outliers. |
| `*_variacion_por_condicion_sin_outliers.csv` | Across-date CV of daily medians, before and after normalization. |
| `*_variacion_por_condicion_con_outliers.csv` | Same, retaining flagged outliers. |
| `*_advertencias.csv` | All warnings and processing notes generated during analysis. |

At the interface level, each file also includes:

- summary cards;
- warnings and processing notes;
- interactive plots for both branches;
- per-file CSV downloads.

The global ZIP export packages the outputs of all uploaded files without mixing experiments.

---

## Tech stack

**Analysis (in-browser via WebAssembly)**

![Pyodide](https://img.shields.io/badge/Pyodide-3776AB?style=flat-square&logo=python&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-013243?style=flat-square&logo=numpy&logoColor=white)
![Plotly](https://img.shields.io/badge/Plotly-3F4F75?style=flat-square&logo=plotly&logoColor=white)

**Frontend**

![HTML5](https://img.shields.io/badge/HTML5-E34F26?style=flat-square&logo=html5&logoColor=white)
![CSS3](https://img.shields.io/badge/CSS3-1572B6?style=flat-square&logo=css3&logoColor=white)
![JavaScript](https://img.shields.io/badge/JavaScript-F7DF1E?style=flat-square&logo=javascript&logoColor=black)
![Bootstrap](https://img.shields.io/badge/Bootstrap-7952B3?style=flat-square&logo=bootstrap&logoColor=white)

**Packaging**

![JSZip](https://img.shields.io/badge/JSZip-555555?style=flat-square)
![FileSaver](https://img.shields.io/badge/FileSaver.js-555555?style=flat-square)

---

## Project structure

```text
zebrafish-ros-normalization-online/
├── index.html                ← main UI
├── style.css                 ← app styling
├── app.js                    ← Pyodide init, UI logic, Plotly rendering, ZIP export
├── zebrafish_ros_engine.py   ← analysis engine running in-browser
├── .nojekyll                 ← GitHub Pages compatibility
└── LICENSE                   ← MIT license
```

---

## Local development

To preview locally:

```bash
python3 -m http.server 8000
```

Then open:

```text
http://localhost:8000
```

---

## GitHub Pages

The project is designed to be deployed as a static GitHub Pages site.

Typical setup:

1. push `main` to GitHub;
2. enable GitHub Pages from the `main` branch and root folder;
3. wait for GitHub to publish the site.

---

## Notes

- The first Pyodide load may take several seconds.
- The app downloads runtime dependencies from public CDNs.
- A detectable date column and at least one condition column are required per file.
- The control does not need to be named `DMSO`, but one control must always be selected before analysis starts.

---

## Author

**Emiliano Balderas Ramírez**  
Bioengineer · PhD Candidate in Biochemical Sciences  
Instituto de Biotecnología (IBt), UNAM

[![LinkedIn](https://img.shields.io/badge/LinkedIn-emilianobalderas-0A66C2?style=flat-square&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/emilianobalderas/)
[![Email](https://img.shields.io/badge/Email-ebalderas%40live.com.mx-D14836?style=flat-square&logo=gmail&logoColor=white)](mailto:ebalderas@live.com.mx)

---

<div align="center"><i>Zebrafish ROS Normalizer Online — upload your CSVs, choose your controls, compare dates correctly.</i></div>
