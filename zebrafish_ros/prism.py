"""Wide tables for GraphPad Prism.

Two shapes, and the choice between them is a methodological decision rather
than a formatting one.

``by_date`` gives one row per acquisition session and one column per condition,
each cell the median log2FC of that session. N is the number of sessions, which
is the design's true number of independent replicates. This is the table to
paste into Prism.

``by_embryo`` gives one row per embryo with all sessions pooled. N is then the
number of embryos, typically 4 to 10 times larger. Feeding it to a t test or an
ANOVA without modelling the nesting inflates significance. It is exported for
inspection and for designs that account for the structure explicitly.
"""

from __future__ import annotations

from collections import defaultdict

from .normalize import Normalized
from .stats import DateFold
from .tidy import DEFAULT_CONTROL


def _ordered_treatments(rows: list[Normalized], control: str) -> list[str]:
    """Control first, then the remaining conditions in order of appearance."""
    seen: list[str] = []
    for row in rows:
        if row.treatment not in seen:
            seen.append(row.treatment)
    if control in seen:
        seen.remove(control)
        seen.insert(0, control)
    return seen


def by_date(
    rows: list[Normalized],
    folds: list[DateFold],
    control: str = DEFAULT_CONTROL,
) -> tuple[list[str], list[dict]]:
    """One row per acquisition date, cells are the median log2FC."""
    treatments = _ordered_treatments(rows, control)
    headers = ["DATE"] + treatments

    lookup = {(f.date, f.treatment): f.median_log2fc for f in folds}
    dates = sorted({f.date for f in folds})

    table = []
    for date in dates:
        record: dict[str, object] = {"DATE": date}
        for treatment in treatments:
            value = lookup.get((date, treatment))
            record[treatment] = "" if value is None else value
        table.append(record)
    return headers, table


def by_embryo(
    rows: list[Normalized], control: str = DEFAULT_CONTROL
) -> tuple[list[str], list[dict]]:
    """One row per embryo, columns padded to the longest condition."""
    treatments = _ordered_treatments(rows, control)

    values: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        values[row.treatment].append(row.log2_norm)

    longest = max((len(v) for v in values.values()), default=0)
    table = []
    for index in range(longest):
        record: dict[str, object] = {}
        for treatment in treatments:
            column = values.get(treatment, [])
            record[treatment] = column[index] if index < len(column) else ""
        table.append(record)
    return treatments, table
