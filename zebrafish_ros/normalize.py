"""Within-session normalization of DCF fluorescence intensities.

The implicit model is multiplicative::

    I_i = mu_t * beta_(group,date) * eps_i

``beta`` absorbs everything that changes between acquisition sessions (laser
power, detector gain, dye loading, embryo stage, operator) and ``mu_t`` is the
treatment effect. A summary of the control from the same group and date
estimates ``beta``, so dividing by it cancels the term and leaves a
dimensionless fold change::

    ratio_i = I_i / anchor(group, date)
    log2fc_i = log2(ratio_i)

The anchor is the **median** of the control by default, which is what the
browser app uses and what suits raw intensities with occasional extreme values.
``--anchor mean`` switches to the arithmetic mean, which is what
`hyper-normalizer` uses on ratiometric data. The choice moves the displayed
fold changes slightly and does not affect the statistical test, which works on
within-session differences of raw logs (see :mod:`zebrafish_ros.stats`).

Two consequences follow, and both constrain the analysis downstream:

1. With the median anchor, the normalized control has a median of 1 in every
   session. With the mean anchor it has a mean of exactly 1. Either way that
   value is the reference line in the figures, and the control group has lost a
   degree of freedom per session.
2. Every embryo from one session shares the same denominator, so embryos within
   a session are not independent of one another.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from statistics import fmean

import numpy as np

from .tidy import DEFAULT_CONTROL, Measurement

#: Below this many control embryos the anchor is reported as poorly determined.
MIN_CONTROL_N_RECOMMENDED = 3

ANCHOR_CHOICES = ("median", "mean")


@dataclass(frozen=True)
class Normalized:
    """A measurement with its fold change against the control of its own date."""

    group: str
    panel: str
    date: str
    treatment: str
    intensity: float         # raw DCF intensity
    ratio_norm: float        # intensity / control anchor of (group, date)
    log2_norm: float         # log2(ratio_norm), a symmetric and additive scale
    control_n: int           # control replicates behind the anchor
    control_anchor: float    # the denominator, kept for traceability
    anchor_status: str       # ok | low_control_n
    source: str
    row_number: int


@dataclass(frozen=True)
class Anchor:
    """The control summary of one ``(group, date)``."""

    group: str
    date: str
    control_n: int
    anchor: float
    control_mean: float
    control_median: float
    control_sd: float
    status: str


class NormalizationError(ValueError):
    """Some ``(group, date)`` group has no usable control."""


def control_anchors(
    measurements: list[Measurement],
    control: str = DEFAULT_CONTROL,
    statistic: str = "median",
) -> tuple[dict[tuple[str, str], Anchor], list[str]]:
    """Anchor and control summary for each ``(group, date)``."""
    if statistic not in ANCHOR_CHOICES:
        raise ValueError(f"anchor must be one of {ANCHOR_CHOICES}, got {statistic!r}")

    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for m in measurements:
        if m.treatment == control:
            grouped[(m.group, m.date)].append(m.intensity)

    anchors: dict[tuple[str, str], Anchor] = {}
    warnings: list[str] = []
    for (group, date), values in sorted(grouped.items()):
        if not values:
            continue
        array = np.asarray(values, dtype=float)
        mean = float(np.mean(array))
        median = float(np.median(array))
        sd = float(np.std(array, ddof=1)) if len(values) > 1 else 0.0
        anchor = median if statistic == "median" else mean

        if anchor <= 0:
            raise NormalizationError(
                f"The control anchor for {group}/{date} is {anchor}, "
                "which cannot serve as a denominator."
            )

        status = "ok"
        if len(values) < MIN_CONTROL_N_RECOMMENDED:
            status = "low_control_n"
            warnings.append(
                f"{group}/{date} has only {len(values)} control embryos; "
                f"the anchor is poorly determined (recommended minimum "
                f"{MIN_CONTROL_N_RECOMMENDED})."
            )

        anchors[(group, date)] = Anchor(
            group=group, date=date, control_n=len(values), anchor=anchor,
            control_mean=mean, control_median=median, control_sd=sd, status=status,
        )
    return anchors, warnings


def normalize(
    measurements: list[Measurement],
    control: str = DEFAULT_CONTROL,
    statistic: str = "median",
    strict: bool = False,
) -> tuple[list[Normalized], list[Anchor], list[str]]:
    """Divide each measurement by the control anchor of its ``(group, date)``.

    Sessions without a control are dropped and reported as a warning. With
    ``strict=True`` they raise instead, which is the right setting for a final
    analysis, where a session without a control is a design failure rather than
    a few data points less.
    """
    anchors, warnings = control_anchors(measurements, control=control, statistic=statistic)

    out: list[Normalized] = []
    missing: set[tuple[str, str]] = set()
    for m in measurements:
        key = (m.group, m.date)
        anchor = anchors.get(key)
        if anchor is None:
            missing.add(key)
            continue
        ratio = m.intensity / anchor.anchor
        out.append(
            Normalized(
                group=m.group,
                panel=m.panel,
                date=m.date,
                treatment=m.treatment,
                intensity=m.intensity,
                ratio_norm=ratio,
                log2_norm=math.log2(ratio),
                control_n=anchor.control_n,
                control_anchor=anchor.anchor,
                anchor_status=anchor.status,
                source=m.source,
                row_number=m.row_number,
            )
        )

    missing_warnings = [
        f"No '{control}' control in {group}/{date}: "
        f"{sum(1 for m in measurements if (m.group, m.date) == (group, date))} "
        "measurements dropped."
        for group, date in sorted(missing)
    ]
    if missing and strict:
        raise NormalizationError("; ".join(missing_warnings))

    return out, list(anchors.values()), warnings + missing_warnings
