"""Statistics for a nested design: embryos within acquisition sessions.

This mirrors :mod:`hyper_normalizer.stats` so that both assays in the same
study are analysed the same way and their results can be reported side by side.

Why this module exists
----------------------
Within-session normalization (:mod:`zebrafish_ros.normalize`) removes the batch
effect, but it leaves two problems if a standard test is then applied to the
pooled embryos:

1. Pseudoreplication. Every embryo from one session shares the same anchor, so
   they are not independent replicates. Treating 40 embryos from 5 sessions as
   n = 40 overstates the available information and produces optimistic
   p-values, often by orders of magnitude.
2. Deflated control variance. The normalized control is centred on 1 in every
   session by construction, so it has lost a degree of freedom per session and
   comparing it as a freely varying sample is anti-conservative.

Both are avoided by working on a log scale and using the acquisition session as
the unit of replication::

    delta_d = mean(log2 I for treatment on date d)
            - mean(log2 I for control on date d)

The difference is taken within a session and on the raw intensities, so the
batch factor cancels algebraically and neither the normalization nor the choice
of anchor affects the test. The ``delta_d`` values, one per session, are tested
against zero with a one-sample t test where n is the number of sessions. Their
exponential, ``2 ** mean(delta)``, is the geometric fold change.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, asdict
from statistics import fmean, stdev

from scipy import stats as sps

from .normalize import Normalized
from .tidy import DEFAULT_CONTROL


@dataclass(frozen=True)
class Summary:
    """Descriptives per ``(group, treatment)``, arithmetic and geometric."""

    group: str
    treatment: str
    n_embryos: int
    n_dates: int
    mean: float
    sd: float
    sem: float
    geo_mean: float
    geo_sd: float
    raw_mean: float          # mean of the unnormalized intensity, for reference
    raw_median: float


@dataclass(frozen=True)
class DateFold:
    """Summary of one ``(group, date, treatment)``, which is the real replicate."""

    group: str
    date: str
    treatment: str
    n_embryos: int
    mean_norm: float          # arithmetic mean of ratio_norm that date
    median_log2fc: float      # median log2FC, the value the Prism table carries
    geo_fold: float           # geometric change against the same date's control
    delta_log2: float         # log2(geo_fold), the quantity entering the test


@dataclass(frozen=True)
class TestResult:
    """Treatment against control, with the session as the unit of replication."""

    group: str
    treatment: str
    n_dates: int
    n_embryos: int
    geo_fold: float
    ci95_low: float
    ci95_high: float
    t: float
    df: int
    p_value: float
    p_holm: float
    p_naive_pooled: float


def _geometric(values: list[float]) -> tuple[float, float]:
    """Geometric mean and standard deviation of strictly positive values."""
    logs = [math.log(v) for v in values if v > 0]
    if not logs:
        return math.nan, math.nan
    mean_log = fmean(logs)
    sd_log = stdev(logs) if len(logs) > 1 else 0.0
    return math.exp(mean_log), math.exp(sd_log)


def describe(rows: list[Normalized]) -> list[Summary]:
    """Descriptives per group and treatment.

    Both the arithmetic mean (comparable with the browser app's tables) and the
    geometric mean (appropriate for a ratio) are reported. The wider the gap
    between them, the more skewed the group's distribution.
    """
    grouped: dict[tuple[str, str], list[Normalized]] = defaultdict(list)
    for row in rows:
        grouped[(row.group, row.treatment)].append(row)

    out: list[Summary] = []
    for (group, treatment), rows_in_group in sorted(grouped.items()):
        values = [r.ratio_norm for r in rows_in_group]
        raw = sorted(r.intensity for r in rows_in_group)
        n = len(values)
        mean = fmean(values)
        sd = stdev(values) if n > 1 else 0.0
        geo_mean, geo_sd = _geometric(values)
        middle = raw[n // 2] if n % 2 else (raw[n // 2 - 1] + raw[n // 2]) / 2
        out.append(
            Summary(
                group=group,
                treatment=treatment,
                n_embryos=n,
                n_dates=len({r.date for r in rows_in_group}),
                mean=mean,
                sd=sd,
                sem=sd / math.sqrt(n) if n > 1 else 0.0,
                geo_mean=geo_mean,
                geo_sd=geo_sd,
                raw_mean=fmean(raw),
                raw_median=middle,
            )
        )
    return out


def date_folds(
    rows: list[Normalized], control: str = DEFAULT_CONTROL
) -> list[DateFold]:
    """Collapse each ``(group, date, treatment)`` to a single value.

    This is the level at which the experiment has true replicates. If a
    condition was repeated on five acquisition dates, n = 5, however many
    embryos were measured on each of them.
    """
    grouped: dict[tuple[str, str, str], list[Normalized]] = defaultdict(list)
    for row in rows:
        grouped[(row.group, row.date, row.treatment)].append(row)

    # Mean log2 of each (group, date) control: the paired reference point.
    control_log2: dict[tuple[str, str], float] = {
        (group, date): fmean([math.log2(r.intensity) for r in rows_in_group])
        for (group, date, treatment), rows_in_group in grouped.items()
        if treatment == control
    }

    out: list[DateFold] = []
    for (group, date, treatment), rows_in_group in sorted(grouped.items()):
        reference = control_log2.get((group, date))
        if reference is None:
            continue
        delta = fmean([math.log2(r.intensity) for r in rows_in_group]) - reference
        log2fcs = sorted(r.log2_norm for r in rows_in_group)
        n = len(log2fcs)
        median = log2fcs[n // 2] if n % 2 else (log2fcs[n // 2 - 1] + log2fcs[n // 2]) / 2
        out.append(
            DateFold(
                group=group,
                date=date,
                treatment=treatment,
                n_embryos=n,
                mean_norm=fmean([r.ratio_norm for r in rows_in_group]),
                median_log2fc=median,
                geo_fold=2.0**delta,
                delta_log2=delta,
            )
        )
    return out


def holm(p_values: list[float]) -> list[float]:
    """Holm-Bonferroni correction, preserving the input order."""
    m = len(p_values)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: p_values[i])
    adjusted = [0.0] * m
    running = 0.0
    for rank, index in enumerate(order):
        candidate = (m - rank) * p_values[index]
        running = max(running, candidate)
        adjusted[index] = min(1.0, running)
    return adjusted


def paired_tests(
    rows: list[Normalized],
    control: str = DEFAULT_CONTROL,
    min_embryos: int = 1,
) -> list[TestResult]:
    """One-sample t test on the per-session ``delta_log2`` values, per treatment.

    Requires at least two sessions holding both the treatment and the control.
    P-values are Holm-corrected within each group, treating each group as one
    family of tests. ``p_naive_pooled`` repeats the contrast with every embryo
    as an independent replicate. It is computed to show how much significance
    inflates when the nested structure is ignored, not to be reported.

    ``min_embryos`` drops sessions in which the treatment or its control was
    measured in fewer than that many embryos. Every session carries the same
    weight, so a session resting on one or two embryos counts as much as one
    resting on twenty. Raising the threshold to 3 is a useful sensitivity check.
    """
    folds = date_folds(rows, control=control)

    control_n = {
        (f.group, f.date): f.n_embryos for f in folds if f.treatment == control
    }
    by_treatment: dict[tuple[str, str], list[DateFold]] = defaultdict(list)
    for fold in folds:
        if fold.treatment == control:
            continue
        if fold.n_embryos < min_embryos:
            continue
        if control_n.get((fold.group, fold.date), 0) < min_embryos:
            continue
        by_treatment[(fold.group, fold.treatment)].append(fold)

    pooled: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        pooled[(row.group, row.treatment)].append(row.log2_norm)

    results: list[TestResult] = []
    for (group, treatment), folds_in_group in sorted(by_treatment.items()):
        deltas = [f.delta_log2 for f in folds_in_group]
        n = len(deltas)
        if n < 2:
            continue

        t_stat, p_value = sps.ttest_1samp(deltas, popmean=0.0)
        mean_delta = fmean(deltas)
        sem = stdev(deltas) / math.sqrt(n)
        half_width = sps.t.ppf(0.975, df=n - 1) * sem

        naive = sps.ttest_ind(
            pooled[(group, treatment)], pooled[(group, control)], equal_var=False
        )

        results.append(
            TestResult(
                group=group,
                treatment=treatment,
                n_dates=n,
                n_embryos=sum(f.n_embryos for f in folds_in_group),
                geo_fold=2.0**mean_delta,
                ci95_low=2.0 ** (mean_delta - half_width),
                ci95_high=2.0 ** (mean_delta + half_width),
                t=float(t_stat),
                df=n - 1,
                p_value=float(p_value),
                p_holm=math.nan,  # filled in below, per family
                p_naive_pooled=float(naive.pvalue),
            )
        )

    # Holm correction within each group.
    corrected: list[TestResult] = []
    by_group: dict[str, list[TestResult]] = defaultdict(list)
    for result in results:
        by_group[result.group].append(result)
    for group_results in by_group.values():
        adjusted = holm([r.p_value for r in group_results])
        for result, p_holm in zip(group_results, adjusted):
            corrected.append(TestResult(**{**asdict(result), "p_holm": p_holm}))

    return sorted(corrected, key=lambda r: (r.group, r.treatment))


def variation(rows: list[Normalized], control: str = DEFAULT_CONTROL) -> list[dict]:
    """Across-date coefficient of variation, before and after normalization.

    This is the direct evidence that the normalization did its job: the CV of
    the daily medians should drop for every condition.
    """
    by_date: dict[tuple[str, str, str], list[Normalized]] = defaultdict(list)
    for row in rows:
        by_date[(row.group, row.treatment, row.date)].append(row)

    def _median(values: list[float]) -> float:
        ordered = sorted(values)
        n = len(ordered)
        return ordered[n // 2] if n % 2 else (ordered[n // 2 - 1] + ordered[n // 2]) / 2

    per_condition: dict[tuple[str, str], list[tuple[float, float]]] = defaultdict(list)
    for (group, treatment, _date), rows_in_date in by_date.items():
        per_condition[(group, treatment)].append(
            (
                _median([r.intensity for r in rows_in_date]),
                _median([r.ratio_norm for r in rows_in_date]),
            )
        )

    out: list[dict] = []
    for (group, treatment), pairs in sorted(per_condition.items()):
        if len(pairs) < 2:
            continue
        raw = [p[0] for p in pairs]
        norm = [p[1] for p in pairs]
        raw_cv = stdev(raw) / fmean(raw) if fmean(raw) else math.nan
        norm_cv = stdev(norm) / fmean(norm) if fmean(norm) else math.nan
        out.append(
            {
                "group": group,
                "treatment": treatment,
                "n_dates": len(pairs),
                "raw_daily_median_cv": raw_cv,
                "normalized_daily_median_cv": norm_cv,
                "cv_reduction": raw_cv - norm_cv,
            }
        )
    return out


def mixed_model(
    rows: list[Normalized], control: str = DEFAULT_CONTROL
) -> dict[str, object]:
    """Mixed model ``log2(I) ~ treatment + (1 | date)``, fitted per group.

    This is the full version of the analysis. Instead of normalizing first and
    testing afterwards, it estimates the session effect and the treatment
    effect at once, and uses every embryo without assuming independence.
    Requires ``statsmodels``. If it is not installed, a note is returned rather
    than raising.
    """
    try:
        import pandas as pd
        import statsmodels.formula.api as smf
    except ImportError:
        return {
            "available": False,
            "message": (
                "Install 'statsmodels' and 'pandas' for the mixed model: "
                "pip install statsmodels pandas"
            ),
        }

    frame = pd.DataFrame(
        {
            "group": [r.group for r in rows],
            "date": [r.date for r in rows],
            "treatment": [r.treatment for r in rows],
            "log2_intensity": [math.log2(r.intensity) for r in rows],
        }
    )

    fits: dict[str, object] = {"available": True}
    for group, subset in frame.groupby("group"):
        if subset["date"].nunique() < 3:
            fits[group] = "Fewer than 3 acquisition dates: mixed model not fitted."
            continue
        formula = f'log2_intensity ~ C(treatment, Treatment(reference="{control}"))'
        model = smf.mixedlm(formula, subset, groups=subset["date"])
        try:
            fit = model.fit(reml=True, method="lbfgs")
        except Exception as error:  # report convergence failures, keep the pipeline alive
            fits[group] = f"Model did not converge: {error}"
            continue
        fits[group] = fit.summary().as_text()
    return fits
