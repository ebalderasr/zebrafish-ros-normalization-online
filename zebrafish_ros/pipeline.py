"""Orchestration: from raw CSVs to output tables and figures.

The pipeline runs one or two branches. ``keep`` retains every embryo; ``drop``
removes the ones flagged by Tukey's rule and recomputes the normalization from
the filtered data. Recomputing matters: dropping an outlier from a control
group changes that session's anchor, so filtering after normalizing would leave
the remaining embryos divided by a denominator that no longer exists.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from .normalize import Anchor, Normalized, normalize
from .outliers import Outlier, flag_outliers
from .prism import by_date, by_embryo
from .stats import (
    DateFold, Summary, TestResult, date_folds, describe, paired_tests, variation,
)
from .tidy import DEFAULT_CONTROL, Measurement, build_tidy, discover_inputs

BRANCHES = ("keep", "drop", "both")

#: Directory name for each branch inside the output directory.
BRANCH_DIRS = {"keep": "with_outliers", "drop": "without_outliers"}


@dataclass
class Branch:
    """One analysis branch: every embryo, or outliers removed."""

    label: str
    normalized: list[Normalized]
    anchors: list[Anchor]
    summary: list[Summary]
    folds: list[DateFold]
    tests: list[TestResult]
    variation: list[dict]


@dataclass
class PipelineResult:
    """Everything one run produces, in memory and already written to disk."""

    tidy: list[Measurement]
    outliers: list[Outlier]
    branches: dict[str, Branch]
    warnings: list[str]
    outputs: list[Path]


def _write_csv(path: Path, rows: list, float_fmt: str = "{:.12g}") -> Path:
    """Write a list of dataclasses to CSV, using their fields as the header.

    Floats are written with 12 significant digits. Rounding to the decimals
    used for display would make any later re-analysis of these files round a
    second time, which shifts values sitting on a rounding boundary.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return path

    columns = [f.name for f in fields(rows[0])]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(col.upper() for col in columns)
        for row in rows:
            record = asdict(row)
            writer.writerow(
                float_fmt.format(record[c]) if isinstance(record[c], float) else record[c]
                for c in columns
            )
    return path


def _write_table(path: Path, headers: list[str], table: list[dict],
                 float_fmt: str = "{:.12g}") -> Path:
    """Write a wide table given explicit headers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for record in table:
            writer.writerow(
                float_fmt.format(record[h]) if isinstance(record.get(h), float)
                else record.get(h, "")
                for h in headers
            )
    return path


def _write_dicts(path: Path, rows: list[dict], float_fmt: str = "{:.12g}") -> Path:
    """Write a list of plain dicts, using the first record's keys as header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return path
    headers = list(rows[0].keys())
    return _write_table(path, [h.upper() for h in headers],
                        [{h.upper(): r[h] for h in headers} for r in rows], float_fmt)


def _build_branch(
    label: str,
    measurements: list[Measurement],
    control: str,
    anchor: str,
    strict: bool,
    min_embryos: int,
) -> tuple[Branch, list[str]]:
    normalized, anchors, warnings = normalize(
        measurements, control=control, statistic=anchor, strict=strict
    )
    branch = Branch(
        label=label,
        normalized=normalized,
        anchors=anchors,
        summary=describe(normalized),
        folds=date_folds(normalized, control=control),
        tests=paired_tests(normalized, control=control, min_embryos=min_embryos),
        variation=variation(normalized, control=control),
    )
    return branch, [f"[{label}] {w}" for w in warnings]


def run(
    input_dir: Path,
    output_dir: Path,
    control: str = DEFAULT_CONTROL,
    anchor: str = "median",
    outliers: str = "both",
    make_plots: bool = True,
    strict: bool = False,
    min_embryos: int = 1,
) -> PipelineResult:
    """Run the full pipeline over a directory of raw CSVs.

    ``outliers`` selects which branches to produce: ``keep`` retains every
    embryo, ``drop`` removes those flagged within their own date and condition,
    and ``both`` writes each into its own subdirectory.
    """
    if outliers not in BRANCHES:
        raise ValueError(f"outliers must be one of {BRANCHES}, got {outliers!r}")

    paths = discover_inputs(input_dir)
    tidy, warnings = build_tidy(paths, control=control)
    kept, flagged = flag_outliers(tidy)

    wanted = ("keep", "drop") if outliers == "both" else (outliers,)
    branches: dict[str, Branch] = {}
    outputs: list[Path] = [
        _write_csv(output_dir / "tidy.csv", tidy),
        _write_csv(output_dir / "outliers_flagged.csv", flagged),
    ]

    for label in wanted:
        rows = tidy if label == "keep" else kept
        branch, branch_warnings = _build_branch(
            label, rows, control, anchor, strict, min_embryos
        )
        branches[label] = branch
        warnings.extend(branch_warnings)

        base = output_dir / BRANCH_DIRS[label] if outliers == "both" else output_dir
        outputs += [
            _write_csv(base / "normalized.csv", branch.normalized),
            _write_csv(base / "control_anchors.csv", branch.anchors),
            _write_csv(base / "summary.csv", branch.summary),
            _write_csv(base / "date_folds.csv", branch.folds),
            _write_csv(base / "tests.csv", branch.tests),
            _write_dicts(base / "variation.csv", branch.variation),
        ]

        headers, table = by_date(branch.normalized, branch.folds, control=control)
        outputs.append(_write_table(base / "prism_by_date.csv", headers, table))
        headers, table = by_embryo(branch.normalized, control=control)
        outputs.append(_write_table(base / "prism_by_embryo.csv", headers, table))

        if make_plots:
            from .plots import plot_control_drift, plot_groups

            outputs.append(plot_groups(branch.normalized, base / "figure_groups.png", control))
            outputs.append(plot_control_drift(branch.anchors, base / "figure_control_drift.png"))

    if flagged:
        warnings.append(
            f"{len(flagged)} of {len(tidy)} embryos "
            f"({100 * len(flagged) / len(tidy):.1f} %) fall outside the 1.5x IQR "
            "bounds of their own date and condition."
        )

    return PipelineResult(
        tidy=tidy, outliers=flagged, branches=branches,
        warnings=warnings, outputs=outputs,
    )
