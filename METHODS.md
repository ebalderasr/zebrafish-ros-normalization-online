# Methods

Manuscript-ready description of the data processing, normalization, and
statistical analysis implemented in this repository. Three versions are
provided: an extended version for the Methods section of a research article, a
compact version for short-format journals, and a single paragraph for a strict
word limit.

The analysis mirrors the one in
[hyper-normalizer](https://github.com/ebalderasr/hyper-normalizer), which
handles the ratiometric HyPer measurements of the same study. Both assays are
normalized against the control acquired in the same session and analysed with
the session as the unit of replication, so the two sets of results can be
reported side by side. The differences specific to DCF, namely the median
anchor and the outlier branch, are stated where they arise.

Text in square brackets marks acquisition details that the code does not
determine and that the authors must supply. Sample sizes refer to the dataset
processed at the time of writing and should be checked against the final
dataset.

---

## Extended version

### Fluorescence quantification and data structure

Intracellular reactive oxygen species were measured with the fluorogenic probe
2',7'-dichlorodihydrofluorescein diacetate (H<sub>2</sub>DCFDA), which yields
the fluorescent product DCF upon oxidation. For each embryo, a single
fluorescence intensity was quantified [*add probe concentration and incubation
time, microscope, objective, excitation and emission settings, and
region-of-interest definition*].

Unlike a ratiometric sensor, DCF reports a single channel, so its absolute
intensity carries no internal correction for probe loading, embryo size, or
illumination. Values are therefore comparable only among embryos imaged in the
same session, which is what makes the within-session anchor described below
necessary rather than merely convenient.

The individual embryo was the unit of measurement. The independent acquisition
session, defined as a cohort of embryos stained, treated, and imaged on the
same day, was the unit of replication. Values were exported as one matrix per
experimental group and per pharmacological panel, with rows corresponding to
individual embryos and columns to conditions, and converted to long format in
which each record comprises the group, acquisition date, condition, and
intensity of a single embryo. Rows are not paired: two values on the same row
are two different embryos imaged on the same date, not one embryo measured
twice. Not every compound was tested in every session, so absent measurements
were missing by design; they were excluded and not imputed. Where one group was
split across panels, the shared vehicle control was recorded in each panel's
sheet and duplicated records were collapsed to a single series per group and
session before analysis. The final dataset comprised [568] embryos across
[2] experimental groups and [7] acquisition sessions, with each compound tested
in [4 to 5] sessions.

### Outlier handling

Within each group, acquisition date, and condition, embryos falling outside
Tukey's bounds were flagged:

$$L_{\text{low}} = Q_1 - 1.5\,\mathrm{IQR}, \qquad
L_{\text{high}} = Q_3 + 1.5\,\mathrm{IQR}, \qquad
\mathrm{IQR} = Q_3 - Q_1.$$

Groups of fewer than four embryos were left unflagged, since their quartiles
are interpolated from too few points for the rule to be informative. Two
parallel analyses were carried out and are reported in full: one retaining
every embryo and one excluding the flagged ones. Because removing an embryo
from a control group changes that session's anchor, the normalization was
recomputed from the filtered data rather than applied first and filtered
afterwards. [*State which branch the reported figures use, and confirm that the
conclusions hold in both.*] In the present dataset [28] of [568] embryos
([4.9] %) were flagged.

### Within-session normalization

Absolute DCF intensity varied between acquisition sessions as a consequence of
probe loading, illumination power, detector gain, and developmental stage.
Across the present dataset the session median of the vehicle control varied by
up to [2.6]-fold, which exceeds the treatment effects of interest. Uncorrected,
this variation would dominate the between-group comparison. The variation is
multiplicative, and the raw intensity of embryo *i* was modelled as

$$I_i = \mu_t \cdot \beta_{g,d} \cdot \varepsilon_i,$$

where $\mu_t$ denotes the effect of condition *t*, $\beta_{g,d}$ a
session-specific gain factor for group *g* on date *d*, and $\varepsilon_i$ the
residual between-embryo variability. Each measurement was expressed relative to
the median of the vehicle (DMSO) controls acquired in the same group and the
same session,

$$\mathrm{ratio}_i = \frac{I_i}{\tilde{I}^{\,\text{veh}}_{g,d}}, \qquad
\tilde{I}^{\,\text{veh}}_{g,d} = \operatorname{median}\left(I_j : j \in \text{veh}(g,d)\right),$$

and reported on a log base 2 scale as $\log_2 \mathrm{FC}_i$. The median was
preferred over the mean because single-channel intensities carry occasional
extreme values; the mean gives the same conclusions and is available as an
option in the software. The denominator was always the concurrent control of
the same group and session. No pooled or global control mean was used, and no
anchor was shared between groups.

This transformation is invariant to any rescaling of a complete session, so
session-level differences in gain are eliminated rather than reduced. Its
effect was verified directly: the coefficient of variation of the daily medians
across sessions fell from [0.29 to 0.45] before normalization to [0.04 to 0.15]
after it, for every condition. Sessions with fewer than three control embryos
were flagged as having a poorly determined anchor.

Two consequences constrained the subsequent analysis. All embryos from a given
session share the same denominator and are therefore not mutually independent.
The vehicle group is also centred on 1 in every session by construction and
cannot be treated as a freely varying sample.

### Statistical analysis

All analyses were performed on log<sub>2</sub>-transformed values. The
logarithm is the appropriate scale for a strictly positive multiplicative
quantity. It renders the model above additive, makes symmetric intervals
meaningful, since a doubling and a halving are then equidistant from zero, and
removes the upward bias affecting the arithmetic mean of ratios.

The independent acquisition session was used as the unit of replication, which
avoids pseudoreplication. For each group *g*, condition *t*, and session *d*,
the within-session effect was computed as the difference between the mean
log<sub>2</sub> intensity of the treated embryos and that of the vehicle
controls acquired in the same session,

$$\delta_{t,d} = \overline{\log_2 I}\big|_{t,d} - \overline{\log_2 I}\big|_{\text{veh},d}.$$

This difference is taken within a session and on the raw intensities, so
$\log_2 \beta_{g,d}$ cancels algebraically. The estimate therefore depends
neither on the normalization step, which serves presentation, nor on whether
the median or the mean was used as the anchor. The resulting values, one per
session, were tested against zero by a two-tailed one-sample *t* test, with *n*
equal to the number of independent sessions in which the condition and its
concurrent vehicle control were both acquired (*n* = [4 to 5] per compound).
Effects are reported as the geometric fold change $2^{\bar\delta}$ with its
95 % confidence interval, and *p* values were adjusted for multiple comparisons
within each experimental group using the Holm-Bonferroni procedure.

Each session carries equal weight in this test, independently of how many
embryos it contains. As a sensitivity analysis, the test was repeated after
excluding sessions in which fewer than three embryos were measured
[*state whether the conclusions were unchanged*]. An equivalent linear
mixed-effects model,
$\log_2 I \sim \text{condition} + (1\,|\,\text{session})$, was also fitted for
each group with session as a random intercept, estimating the session and
condition effects jointly rather than normalizing before testing. Tests
treating individual embryos as independent replicates were not used. They
disregard the nested structure of the design and yield anti-conservative
*p* values, which for the present dataset were optimistic by up to thirteen
orders of magnitude.

### Data presentation

Results are displayed as SuperPlots (Lord et al., 2020). Individual embryo
values are plotted as small background points, and the mean of each independent
session is superimposed as a large colour-coded marker, so that the number of
independent replicates and the between-session consistency of each effect are
visible. Boxes indicate the median and interquartile range of the individual
embryo values, and the dashed line marks the normalized vehicle control at 1. A
companion figure plots the raw control anchor of each session, which is the
between-session drift that the normalization removes.

### Software and data availability

The analysis is implemented twice from a single specification: as a
browser-based application that runs entirely client-side through Pyodide, with
no data leaving the user's machine, and as an offline Python package for
scripted and reproducible use. Both are available at
<https://github.com/ebalderasr/zebrafish-ros-normalization-online> (version
1.0.0, MIT licence), and the application is deployed at
<https://ebalderasr.github.io/zebrafish-ros-normalization-online/>. The offline
package uses Python 3.10 with NumPy, SciPy, pandas, Matplotlib, seaborn, and
statsmodels. An automated test suite verifies that the two implementations
produce identical per-embryo log<sub>2</sub>FC values, control anchors, outlier
flags, and per-session tables, and checks the defining invariants of the
normalization, including invariance to session-level rescaling.

**Reference.** Lord SJ, Velle KB, Mullins RD, Fritz-Laylin LK. SuperPlots:
Communicating reproducibility and variability in cell biology. *J Cell Biol.*
2020;219(6):e202001064. doi:10.1083/jcb.202001064

---

## Compact version

Intracellular reactive oxygen species were measured with H<sub>2</sub>DCFDA,
and one DCF fluorescence intensity was quantified per embryo. Within each
group, acquisition date, and condition, embryos outside Tukey's 1.5 IQR bounds
were flagged, and the analysis was carried out both retaining and excluding
them. Because single-channel DCF intensity carries no internal correction and
drifts between acquisition sessions, each measurement was normalized to the
median of the vehicle (DMSO) controls acquired in the same group and the same
session, yielding a log<sub>2</sub> fold change relative to vehicle. The anchor
was recomputed after outlier removal, since dropping a control embryo changes
that session's denominator.

Statistical analysis used the independent acquisition session, rather than the
individual embryo, as the unit of replication. For each session, the difference
between the mean log<sub>2</sub> intensity of treated embryos and that of the
concurrent vehicle controls was computed; these session-level differences were
tested against zero by a two-tailed one-sample *t* test (*n* = [4 to 5]
independent sessions per compound) with Holm-Bonferroni correction within each
group. Effects are reported as geometric fold changes with 95 % confidence
intervals and displayed as SuperPlots showing per-session means superimposed on
individual embryo values. Analyses were performed in Python, in a browser
application and an equivalent offline package whose agreement is covered by
automated tests; both are available at
<https://github.com/ebalderasr/zebrafish-ros-normalization-online>.

---

## Single-paragraph version

For contexts with a strict word limit.

DCF fluorescence intensity was quantified per embryo, embryos outside Tukey's
1.5 IQR bounds within their own group, date and condition were flagged, and
each measurement was normalized to the median of the vehicle (DMSO) controls
acquired in the same group and the same acquisition session, which removes the
multiplicative session-to-session drift in absolute signal and yields a
log<sub>2</sub> fold change relative to vehicle. Inference used the independent
acquisition session as the unit of replication: per-session differences between
treated embryos and concurrent vehicle controls were tested against zero by a
two-tailed one-sample *t* test (*n* = [4 to 5] sessions per compound) with
Holm-Bonferroni correction within each group, and effects are reported as
geometric fold changes with 95 % confidence intervals. Analysis code is
available at
<https://github.com/ebalderasr/zebrafish-ros-normalization-online>.
