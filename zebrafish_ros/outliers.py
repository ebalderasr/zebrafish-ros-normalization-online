"""Outlier flagging within each ``(group, date, treatment)``.

Raw DCF intensity is more prone to isolated extreme values than a ratiometric
readout, because nothing in the measurement divides out embryo-to-embryo
differences in dye loading. The browser app therefore reports two branches, one
keeping every embryo and one dropping the flagged ones, and this module
reproduces that behaviour.

Flagging uses Tukey's rule on the raw intensities of each group::

    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR

Groups of fewer than four embryos are left unflagged: with three points the
quartiles are interpolated from too little information for the rule to mean
anything.

Removing outliers changes the control anchor of the affected date, so the
normalization has to be recomputed on the filtered data rather than applied and
then filtered. :func:`zebrafish_ros.pipeline.run` does that.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np

from .tidy import Measurement

#: Tukey's multiplier. 1.5 is the convention behind the boxplot whisker.
IQR_MULTIPLIER = 1.5

#: Below this many embryos the quartiles carry too little information.
MIN_N_FOR_FLAGGING = 4


@dataclass(frozen=True)
class Outlier:
    """One flagged embryo, with the bounds that excluded it."""

    group: str
    date: str
    treatment: str
    intensity: float
    lower_bound: float
    upper_bound: float
    source: str
    row_number: int


def iqr_bounds(values: list[float]) -> tuple[float, float] | None:
    """Tukey bounds for a group, or ``None`` when the group is too small."""
    if len(values) < MIN_N_FOR_FLAGGING:
        return None
    array = np.asarray(values, dtype=float)
    q1 = float(np.percentile(array, 25))
    q3 = float(np.percentile(array, 75))
    spread = q3 - q1
    return q1 - IQR_MULTIPLIER * spread, q3 + IQR_MULTIPLIER * spread


def flag_outliers(
    measurements: list[Measurement],
) -> tuple[list[Measurement], list[Outlier]]:
    """Split the measurements into those kept and those flagged."""
    grouped: dict[tuple[str, str, str], list[Measurement]] = defaultdict(list)
    for m in measurements:
        grouped[(m.group, m.date, m.treatment)].append(m)

    bounds = {
        key: iqr_bounds([m.intensity for m in group])
        for key, group in grouped.items()
    }

    kept: list[Measurement] = []
    flagged: list[Outlier] = []
    for m in measurements:
        limits = bounds[(m.group, m.date, m.treatment)]
        if limits is None:
            kept.append(m)
            continue
        lower, upper = limits
        if m.intensity < lower or m.intensity > upper:
            flagged.append(
                Outlier(
                    group=m.group,
                    date=m.date,
                    treatment=m.treatment,
                    intensity=m.intensity,
                    lower_bound=lower,
                    upper_bound=upper,
                    source=m.source,
                    row_number=m.row_number,
                )
            )
        else:
            kept.append(m)

    return kept, flagged
