#!/usr/bin/env python3
"""Generate the SYNTHETIC CSVs in `data/example/`.

The real experimental data is not distributed with this repository. This script
builds a dataset with the same structure as the laboratory files, and with the
same layout `hyper-normalizer` expects, so one experiment can be exported once
and processed by either pipeline.

The simulator uses the same multiplicative model the normalization assumes::

    I = mu_condition * beta_session * lognormal(0, sigma)

It includes three features on purpose:

* strong between-session drift in the raw intensity, which is what the
  within-session anchor is meant to remove;
* a control shared across the two drug panels of each group, which exercises
  the de-duplication;
* a handful of extreme embryos, which exercise the 1.5x IQR branch.

Usage::

    python scripts/make_example_data.py [--output-dir data/example] [--seed 20260908]
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

# Fictional dates, in the same YYMMDD format the laboratory uses.
DATES = ["240115", "240122", "240129", "240205", "240212", "240219", "240226"]

PANELS = {
    "ANTINOX": ["DMSO", "VAS", "DPI", "APO"],
    "ANTIox": ["DMSO", "EUK", "TEMPOL", "NAC"],
}

#: Embryos per (session, condition). Unbalanced on purpose.
DESIGN: dict[str, dict[str, int]] = {
    "240115": {"DMSO": 8, "VAS": 8, "DPI": 8, "APO": 8},
    "240122": {"DMSO": 10, "VAS": 9, "DPI": 10, "APO": 10, "EUK": 8, "TEMPOL": 8},
    "240129": {"DMSO": 9, "VAS": 8, "APO": 9, "EUK": 9, "TEMPOL": 9, "NAC": 7},
    "240205": {"DMSO": 10, "VAS": 10, "DPI": 9, "APO": 10, "EUK": 10, "TEMPOL": 10, "NAC": 8},
    "240212": {"DMSO": 8, "DPI": 8, "APO": 8, "EUK": 8, "TEMPOL": 8, "NAC": 8},
    "240219": {"DMSO": 9, "NAC": 9, "EUK": 9},
    "240226": {"DMSO": 2, "NAC": 2},          # a thin session, on purpose
}

#: True effect of each compound, as a fraction of the control signal.
TRUE_EFFECT = {
    "DMSO": 1.00, "VAS": 0.86, "DPI": 0.74, "APO": 0.80,
    "EUK": 0.88, "TEMPOL": 0.83, "NAC": 0.69,
}

#: Session gain. DCF intensity drifts far more than a ratiometric readout.
SESSION_EFFECT = {
    "240115": 1.00, "240122": 1.62, "240129": 0.78, "240205": 1.21,
    "240212": 2.05, "240219": 0.91, "240226": 1.44,
}

GROUPS = {"WT": 118.0, "MUT": 141.0}  # baseline DCF intensity of each group
NOISE_SD = 0.16                        # sigma of the per-embryo lognormal noise

#: Extreme embryos planted so the outlier branch has something to remove.
OUTLIERS = {
    ("WT", "240122", "VAS"): 3.1,
    ("MUT", "240205", "APO"): 0.28,
    ("WT", "240212", "NAC"): 2.7,
}


def simulate(rng: random.Random) -> dict[tuple[str, str, str], list[float]]:
    """Return the simulated intensities per ``(group, session, condition)``."""
    values: dict[tuple[str, str, str], list[float]] = {}
    for group, baseline in GROUPS.items():
        for date, conditions in DESIGN.items():
            for condition, n in conditions.items():
                center = baseline * TRUE_EFFECT[condition] * SESSION_EFFECT[date]
                series = [
                    round(center * rng.lognormvariate(0.0, NOISE_SD), 2)
                    for _ in range(n)
                ]
                factor = OUTLIERS.get((group, date, condition))
                if factor is not None and series:
                    series[0] = round(series[0] * factor, 2)
                values[(group, date, condition)] = series
    return values


def write_panel(
    path: Path,
    group: str,
    conditions: list[str],
    values: dict[tuple[str, str, str], list[float]],
) -> None:
    """Write one wide CSV, leaving cells empty where nothing was measured.

    The control is written into every panel that shares a session, with the
    same values. That is how the data leaves the laboratory, and it is what the
    de-duplication in `zebrafish_ros.tidy` resolves.
    """
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["FECHA", *conditions])

        for date in DATES:
            present = {c: values.get((group, date, c), []) for c in conditions}
            n_rows = max((len(v) for v in present.values()), default=0)
            if n_rows == 0:
                continue  # this panel was not run that session
            for i in range(n_rows):
                row = [date]
                for condition in conditions:
                    series = present[condition]
                    row.append(f"{series[i]:.2f}" if i < len(series) else "")
                writer.writerow(row)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output-dir", type=Path, default=Path("data/example"))
    parser.add_argument("--seed", type=int, default=20260908)
    args = parser.parse_args(argv)

    rng = random.Random(args.seed)
    values = simulate(rng)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for panel, conditions in PANELS.items():
        for group in GROUPS:
            path = args.output_dir / f"{group} DCF {panel}.csv"
            write_panel(path, group, conditions, values)
            print(f"written: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
