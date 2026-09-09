"""The offline package must agree with the browser engine, value by value.

The app published on GitHub Pages was used to process data for publication, so
the offline pipeline is only useful if it reproduces it. These tests run
`zebrafish_ros_engine.analyze_one_file`, which is the module Pyodide loads in
the browser, and compare its output against `zebrafish_ros` on the same file.

The comparison covers the quantities that reach a figure or a manuscript: the
per-embryo log2FC, the control anchor of each date, and the outlier flags.
"""

from __future__ import annotations

import shutil
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE_DIR = ROOT / "data" / "example"

sys.path.insert(0, str(ROOT))  # the browser engine lives at the repository root

engine = pytest.importorskip("zebrafish_ros_engine")

from zebrafish_ros.normalize import normalize
from zebrafish_ros.outliers import flag_outliers
from zebrafish_ros.tidy import build_tidy

EXAMPLE_FILES = sorted(p.name for p in EXAMPLE_DIR.glob("*.csv"))

#: Branch names differ between the two implementations.
BRANCHES = {"keep": "with_outliers", "drop": "without_outliers"}


def _offline(path: Path, drop_outliers: bool):
    """Run the offline package over a directory holding only this file."""
    measurements, _ = build_tidy([path])
    if drop_outliers:
        measurements, _ = flag_outliers(measurements)
    rows, anchors, _ = normalize(measurements, statistic="median")
    return rows, anchors


def _engine(path: Path, branch: str):
    result = engine.analyze_one_file(path.name, path.read_text(encoding="utf-8"), "DMSO")
    return result["branches"][branch]


@pytest.mark.parametrize("filename", EXAMPLE_FILES)
@pytest.mark.parametrize("offline_branch", ["keep", "drop"])
def test_log2fc_matches_the_browser_engine(filename, offline_branch):
    path = EXAMPLE_DIR / filename
    rows, _ = _offline(path, drop_outliers=(offline_branch == "drop"))
    branch = _engine(path, BRANCHES[offline_branch])

    mine = Counter(
        (r.date, r.treatment, round(r.log2_norm, 10)) for r in rows
    )
    theirs = Counter(
        (r["acquisition_date"], r["condition_original"], round(r["log2fc_vs_control"], 10))
        for r in branch["long_rows"]
        if r["log2fc_vs_control"] is not None
    )
    assert mine == theirs, f"{filename}/{offline_branch}: log2FC values differ"


@pytest.mark.parametrize("filename", EXAMPLE_FILES)
def test_control_anchor_matches_the_browser_engine(filename):
    path = EXAMPLE_DIR / filename
    _, anchors = _offline(path, drop_outliers=False)
    branch = _engine(path, "with_outliers")

    mine = {a.date: (a.control_n, round(a.anchor, 10)) for a in anchors}
    theirs = {
        r["acquisition_date"]: (r["control_n"], round(r["control_median"], 10))
        for r in branch["control_rows"]
        if r["control_median"] is not None
    }
    assert mine == theirs, f"{filename}: control anchors differ"


@pytest.mark.parametrize("filename", EXAMPLE_FILES)
def test_outlier_flags_match_the_browser_engine(filename):
    path = EXAMPLE_DIR / filename
    measurements, _ = build_tidy([path])
    _, flagged = flag_outliers(measurements)

    result = engine.analyze_one_file(path.name, path.read_text(encoding="utf-8"), "DMSO")
    theirs = Counter(
        (r["acquisition_date"], r["condition_original"], round(r["intensity"], 10))
        for r in result["branches"]["with_outliers"]["long_rows"]
        if r["is_iqr_outlier_within_date_condition"]
    )
    mine = Counter(
        (o.date, o.treatment, round(o.intensity, 10)) for o in flagged
    )
    assert mine == theirs, f"{filename}: outlier flags differ"


def test_prism_by_date_matches_the_browser_engine():
    """The table meant for statistical inference must carry the same values.

    Column order is the one documented difference. The browser engine sorts its
    long rows by ``(date, row number, condition)`` before building the table, so
    its Prism columns come out alphabetical. The offline package keeps the order
    the conditions have in the CSV header, with the control first, which is the
    order the experimenter wrote and the one `hyper-normalizer` uses. Every cell
    value is identical; only the columns sit in a different sequence.
    """
    from zebrafish_ros.prism import by_date
    from zebrafish_ros.stats import date_folds

    for filename in EXAMPLE_FILES:
        path = EXAMPLE_DIR / filename
        rows, _ = _offline(path, drop_outliers=True)
        headers, table = by_date(rows, date_folds(rows))

        branch = _engine(path, "without_outliers")
        their_rows, their_headers = engine.build_prism_by_date(branch["long_rows"])

        assert set(headers[1:]) == set(their_headers[1:]), f"{filename}: columns differ"
        assert headers[1] == their_headers[1] == "DMSO", (
            f"{filename}: the control must lead both tables"
        )

        mine = {
            (row["DATE"], column): row[column]
            for row in table for column in headers[1:]
        }
        theirs = {
            (row["acquisition_date"], column): row[column]
            for row in their_rows for column in their_headers[1:]
        }
        assert set(mine) == set(theirs), f"{filename}: date/condition cells differ"
        for key, value in mine.items():
            other = theirs[key]
            if value == "" or other == "":
                assert value == other == "", f"{filename}/{key}: one side is empty"
            else:
                assert value == pytest.approx(other, rel=1e-12), f"{filename}/{key}"
