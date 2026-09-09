"""Tests of the pipeline's mathematical properties.

These do not check experimental numbers. They check the invariants the
normalization claims to satisfy. If one of them breaks, the analysis no longer
means what the README says it means.
"""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path

import pytest

from zebrafish_ros.normalize import NormalizationError, control_anchors, normalize
from zebrafish_ros.outliers import flag_outliers, iqr_bounds
from zebrafish_ros.stats import date_folds, describe, holm, paired_tests, variation
from zebrafish_ros.tidy import (
    InputError,
    Measurement,
    build_tidy,
    deduplicate_controls,
    discover_inputs,
    parse_date,
    parse_filename,
    parse_measurement,
    read_wide,
)

EXAMPLE_DIR = Path(__file__).resolve().parent.parent / "data" / "example"


def _measurement(**kwargs) -> Measurement:
    base = dict(
        group="WT", panel="ANTINOX", date="2024-01-15", date_raw="240115",
        treatment="DMSO", intensity=100.0, source="X.csv", row_number=2,
    )
    return Measurement(**{**base, **kwargs})


@pytest.fixture(scope="module")
def tidy():
    measurements, warnings = build_tidy(discover_inputs(EXAMPLE_DIR))
    assert warnings == [], f"The example data should produce no warnings: {warnings}"
    return measurements


@pytest.fixture(scope="module")
def normalized(tidy):
    rows, _anchors, _warnings = normalize(tidy, statistic="median")
    return rows


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def test_parse_measurement_accepts_decimal_commas():
    assert parse_measurement("123,4") == pytest.approx(123.4)
    assert parse_measurement("1.234,5") == pytest.approx(1234.5)
    assert parse_measurement("1,234.5") == pytest.approx(1234.5)
    assert parse_measurement(" 98.6 ") == pytest.approx(98.6)


def test_parse_measurement_rejects_blanks_and_null_tokens():
    for value in ("", "   ", "NA", "n/a", "nd", "S/D", None, "abc"):
        assert parse_measurement(value) is None


def test_parse_date_converts_yymmdd_to_iso():
    assert parse_date("240115") == "2024-01-15"
    assert parse_date("240115.0") == "2024-01-15"  # a spreadsheet made it a float
    assert parse_date("24-01-15") == "2024-01-15"


def test_parse_date_rejects_unparseable_values():
    for value in ("", "2024", "abc", "999999", None):
        assert parse_date(value) is None


def test_parse_filename_reads_the_group_and_panel():
    assert parse_filename(Path("WT DCF ANTINOX.csv")) == ("WT", "ANTINOX")
    assert parse_filename(Path("MUT DCF ANTIox.csv")) == ("MUT", "ANTIOX")


def test_parse_filename_falls_back_to_the_stem():
    """A file off the convention is its own group, which is the app's behaviour."""
    assert parse_filename(Path("experiment_1.csv")) == ("experiment_1", "")


def test_read_wide_skips_empty_cells():
    rows = read_wide(EXAMPLE_DIR / "WT DCF ANTINOX.csv")
    assert rows
    assert all(math.isfinite(r.intensity) for r in rows)
    counts = defaultdict(int)
    for r in rows:
        counts[r.treatment] += 1
    assert len(set(counts.values())) > 1, "The design is unbalanced by construction"


# --------------------------------------------------------------------------
# Control de-duplication
# --------------------------------------------------------------------------

def test_control_is_not_counted_twice(tidy):
    """The control appears in both panels of a group and must be counted once."""
    deduped = defaultdict(int)
    for m in tidy:
        if m.treatment == "DMSO":
            deduped[(m.group, m.date)] += 1

    raw = defaultdict(int)
    for path in discover_inputs(EXAMPLE_DIR):
        for m in read_wide(path):
            if m.treatment == "DMSO":
                raw[(m.group, m.date)] += 1

    shared = [k for k in raw if raw[k] > deduped[k]]
    assert shared, "The example must include sessions with a control in both panels"
    for key in shared:
        assert deduped[key] < raw[key]


def test_dedup_warns_when_controls_disagree():
    measurements = [
        _measurement(intensity=100.0, source="A.csv"),
        _measurement(intensity=110.0, source="A.csv"),
        _measurement(intensity=999.0, source="B.csv"),
    ]
    kept, warnings = deduplicate_controls(measurements)
    assert len(kept) == 2
    assert len(warnings) == 1
    assert "differs" in warnings[0]


# --------------------------------------------------------------------------
# Outliers
# --------------------------------------------------------------------------

def test_iqr_bounds_need_at_least_four_values():
    assert iqr_bounds([1.0, 2.0, 3.0]) is None
    assert iqr_bounds([1.0, 2.0, 3.0, 4.0]) is not None


def test_flag_outliers_splits_without_losing_embryos(tidy):
    kept, flagged = flag_outliers(tidy)
    assert len(kept) + len(flagged) == len(tidy)
    assert flagged, "The example data plants extreme embryos on purpose"


def test_flagged_embryos_sit_outside_their_own_bounds(tidy):
    _kept, flagged = flag_outliers(tidy)
    for outlier in flagged:
        assert outlier.intensity < outlier.lower_bound or outlier.intensity > outlier.upper_bound


# --------------------------------------------------------------------------
# Normalization
# --------------------------------------------------------------------------

def test_median_anchor_centres_the_control_on_one(tidy):
    rows, _anchors, _ = normalize(tidy, statistic="median")
    grouped = defaultdict(list)
    for r in rows:
        if r.treatment == "DMSO":
            grouped[(r.group, r.date)].append(r.ratio_norm)

    assert grouped
    for key, values in grouped.items():
        ordered = sorted(values)
        n = len(ordered)
        median = ordered[n // 2] if n % 2 else (ordered[n // 2 - 1] + ordered[n // 2]) / 2
        assert median == pytest.approx(1.0, abs=1e-12), key


def test_mean_anchor_gives_the_control_a_mean_of_one(tidy):
    """With --anchor mean the pipeline matches hyper-normalizer's convention."""
    rows, _anchors, _ = normalize(tidy, statistic="mean")
    grouped = defaultdict(list)
    for r in rows:
        if r.treatment == "DMSO":
            grouped[(r.group, r.date)].append(r.ratio_norm)

    for key, values in grouped.items():
        assert sum(values) / len(values) == pytest.approx(1.0, abs=1e-12), key


def test_normalization_cancels_the_batch_factor(tidy):
    """The central property: rescaling a whole session changes nothing.

    Multiplying every measurement of one session by 7, which is what a change
    in detector gain does, must leave that session's normalized values
    identical.
    """
    original, _, _ = normalize(tidy, statistic="median")

    session = ("WT", "2024-02-05")
    rescaled = [
        Measurement(
            group=m.group, panel=m.panel, date=m.date, date_raw=m.date_raw,
            treatment=m.treatment,
            intensity=m.intensity * 7.0 if (m.group, m.date) == session else m.intensity,
            source=m.source, row_number=m.row_number,
        )
        for m in tidy
    ]
    perturbed, _, _ = normalize(rescaled, statistic="median")

    assert len(original) == len(perturbed)
    for a, b in zip(original, perturbed):
        assert a.ratio_norm == pytest.approx(b.ratio_norm, rel=1e-12)


def test_log2_matches_the_ratio(normalized):
    for r in normalized:
        assert r.log2_norm == pytest.approx(math.log2(r.ratio_norm), rel=1e-12)


def test_normalize_drops_sessions_without_control():
    orphan = [_measurement(treatment="NAC", date="2024-09-09")]
    rows, _anchors, warnings = normalize(orphan)
    assert rows == []
    assert any("No 'DMSO' control" in w for w in warnings)

    with pytest.raises(NormalizationError):
        normalize(orphan, strict=True)


def test_low_control_n_is_reported(tidy):
    """The thin session in the example data must raise a warning, not pass silently."""
    _anchors, warnings = control_anchors(tidy)
    assert any("only 2 control embryos" in w for w in warnings)


def test_anchor_rejects_unknown_statistic(tidy):
    with pytest.raises(ValueError):
        control_anchors(tidy, statistic="mode")


# --------------------------------------------------------------------------
# Nested statistics
# --------------------------------------------------------------------------

def test_control_folds_to_one_every_session(normalized):
    for fold in date_folds(normalized):
        if fold.treatment == "DMSO":
            assert fold.delta_log2 == pytest.approx(0.0, abs=1e-12)
            assert fold.geo_fold == pytest.approx(1.0, abs=1e-12)


def test_test_n_counts_sessions_not_embryos(normalized):
    for result in paired_tests(normalized):
        assert result.df == result.n_dates - 1
        assert result.n_dates < result.n_embryos


def test_naive_test_is_more_optimistic(normalized):
    results = paired_tests(normalized)
    assert results
    inflated = [r for r in results if r.p_naive_pooled < r.p_value]
    assert len(inflated) >= len(results) * 0.7


def test_the_test_ignores_the_anchor_choice(tidy):
    """The contrast uses within-session log differences of raw intensities.

    Median or mean anchor moves the displayed fold changes, and must leave the
    p-values untouched.
    """
    median_rows, _, _ = normalize(tidy, statistic="median")
    mean_rows, _, _ = normalize(tidy, statistic="mean")

    by_median = {(r.group, r.treatment): r for r in paired_tests(median_rows)}
    by_mean = {(r.group, r.treatment): r for r in paired_tests(mean_rows)}

    assert set(by_median) == set(by_mean)
    for key, result in by_median.items():
        assert result.p_value == pytest.approx(by_mean[key].p_value, rel=1e-12), key
        assert result.geo_fold == pytest.approx(by_mean[key].geo_fold, rel=1e-12), key


def test_normalization_reduces_across_date_variation(normalized):
    """The direct evidence that the anchor works."""
    rows = [v for v in variation(normalized) if v["treatment"] != "DMSO"]
    assert rows
    improved = [v for v in rows if v["normalized_daily_median_cv"] < v["raw_daily_median_cv"]]
    assert len(improved) == len(rows), "Normalization should lower the CV of every condition"


def test_geometric_mean_does_not_exceed_arithmetic(normalized):
    for row in describe(normalized):
        assert row.geo_mean <= row.mean + 1e-12


def test_holm_is_monotone_and_conservative():
    raw = [0.001, 0.02, 0.04, 0.5]
    adjusted = holm(raw)
    assert all(a >= b for a, b in zip(adjusted, raw))
    ordered = [adjusted[i] for i in sorted(range(4), key=lambda i: raw[i])]
    assert ordered == sorted(ordered)
    assert all(p <= 1.0 for p in adjusted)


def test_holm_with_empty_list():
    assert holm([]) == []


def test_min_embryos_drops_thin_sessions(normalized):
    loose = {(r.group, r.treatment): r.n_dates for r in paired_tests(normalized)}
    strict = {(r.group, r.treatment): r.n_dates for r in paired_tests(normalized, min_embryos=3)}
    assert strict
    for key, dates in strict.items():
        assert dates <= loose[key]
    assert any(strict[k] < loose[k] for k in strict)


# --------------------------------------------------------------------------
# Pipeline outputs
# --------------------------------------------------------------------------

def test_both_branches_are_written_without_plots(tmp_path):
    from zebrafish_ros.pipeline import run

    result = run(EXAMPLE_DIR, tmp_path, make_plots=False, outliers="both")
    written = {p.relative_to(tmp_path).as_posix() for p in result.outputs}

    assert "tidy.csv" in written
    assert "outliers_flagged.csv" in written
    for directory in ("with_outliers", "without_outliers"):
        for name in ("normalized.csv", "control_anchors.csv", "summary.csv",
                     "date_folds.csv", "tests.csv", "variation.csv",
                     "prism_by_date.csv", "prism_by_embryo.csv"):
            assert f"{directory}/{name}" in written
    assert all(p.exists() for p in result.outputs)


def test_dropping_outliers_changes_the_anchor(tmp_path):
    """Filtering after normalizing would be wrong, so the pipeline recomputes."""
    from zebrafish_ros.pipeline import run

    result = run(EXAMPLE_DIR, tmp_path, make_plots=False, outliers="both")
    keep = {(a.group, a.date): a.anchor for a in result.branches["keep"].anchors}
    drop = {(a.group, a.date): a.anchor for a in result.branches["drop"].anchors}

    assert set(keep) == set(drop)
    assert any(keep[k] != drop[k] for k in keep), (
        "At least one session should have a control outlier that moves its anchor"
    )


def test_rejects_an_unknown_branch(tmp_path):
    from zebrafish_ros.pipeline import run

    with pytest.raises(ValueError):
        run(EXAMPLE_DIR, tmp_path, outliers="sometimes")


def test_empty_directory_is_reported(tmp_path):
    with pytest.raises(InputError):
        discover_inputs(tmp_path)
