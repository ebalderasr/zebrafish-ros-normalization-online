"""Normalization and analysis pipeline for DCF fluorescence in zebrafish embryos.

Flow::

    raw CSVs -> tidy -> outlier flagging -> within-session normalization
             -> statistics -> Prism tables -> figures

Quick use::

    from pathlib import Path
    from zebrafish_ros import run

    result = run(Path("data/example"), Path("results"))
    for test in result.branches["drop"].tests:
        print(test.treatment, round(test.geo_fold, 3), round(test.p_holm, 4))

This package is the offline counterpart of the browser app in this repository.
Both implement the same normalization; the test suite checks that they agree
value by value.
"""

from .normalize import Anchor, Normalized, normalize
from .outliers import Outlier, flag_outliers
from .pipeline import Branch, PipelineResult, run
from .stats import (
    DateFold, Summary, TestResult, date_folds, describe, paired_tests, variation,
)
from .tidy import Measurement, build_tidy, discover_inputs

__version__ = "1.0.0"

__all__ = [
    "Measurement", "Normalized", "Anchor", "Outlier",
    "Summary", "DateFold", "TestResult", "Branch", "PipelineResult",
    "build_tidy", "discover_inputs", "flag_outliers", "normalize",
    "describe", "date_folds", "paired_tests", "variation", "run",
    "__version__",
]
