"""Command-line interface.

    python -m zebrafish_ros --input-dir data/example --output-dir results
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .normalize import ANCHOR_CHOICES
from .pipeline import BRANCHES, run
from .stats import mixed_model
from .tidy import DEFAULT_CONTROL, InputError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="zebrafish-ros",
        description=(
            "Normalize DCF fluorescence intensities against the same-session "
            "control and analyse the result using the acquisition session as "
            "the unit of replication. Offline counterpart of the browser app."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-dir", type=Path, default=Path("data/example"),
        help="Directory holding the raw CSV files.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results"),
        help="Directory to write tables and figures into.",
    )
    parser.add_argument(
        "--control", default=DEFAULT_CONTROL,
        help="Name of the control column, used as the normalization anchor.",
    )
    parser.add_argument(
        "--anchor", choices=ANCHOR_CHOICES, default="median",
        help=(
            "Statistic summarising the control of each session. 'median' matches "
            "the browser app; 'mean' matches hyper-normalizer. The test is "
            "unaffected either way."
        ),
    )
    parser.add_argument(
        "--outliers", choices=BRANCHES, default="both",
        help="Keep every embryo, drop those outside 1.5x IQR, or write both branches.",
    )
    parser.add_argument(
        "--min-embryos", type=int, default=1, metavar="N",
        help=(
            "Drop sessions measured in fewer than N embryos from the test. "
            "Every session carries equal weight, so try 3 as a sensitivity check."
        ),
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="Fail if any (group, date) has no control, instead of dropping it.",
    )
    parser.add_argument(
        "--mixed-model", action="store_true",
        help="Also fit log2(I) ~ treatment + (1|date). Requires statsmodels.",
    )
    parser.add_argument(
        "--no-plots", action="store_true", help="Skip figure generation.",
    )
    return parser


def _print_branch(label: str, branch) -> None:
    heading = {"keep": "all embryos", "drop": "outliers removed"}[label]
    print(f"\n{'=' * 72}\nBranch: {heading}\n{'=' * 72}")

    print("\nDescriptives (normalized intensity)")
    header = (
        f"{'Group':10s} {'Cond.':10s} {'n_emb':>6s} {'dates':>6s} "
        f"{'mean':>7s} {'SD':>7s} {'geo.mean':>9s}"
    )
    print(header)
    print("-" * len(header))
    for row in branch.summary:
        print(
            f"{row.group[:10]:10s} {row.treatment[:10]:10s} {row.n_embryos:6d} "
            f"{row.n_dates:6d} {row.mean:7.3f} {row.sd:7.3f} {row.geo_mean:9.3f}"
        )

    if branch.variation:
        print("\nAcross-date CV of the daily medians, before and after normalization")
        header = f"{'Group':10s} {'Cond.':10s} {'raw CV':>8s} {'norm CV':>8s} {'change':>8s}"
        print(header)
        print("-" * len(header))
        for row in branch.variation:
            print(
                f"{row['group'][:10]:10s} {row['treatment'][:10]:10s} "
                f"{row['raw_daily_median_cv']:8.3f} "
                f"{row['normalized_daily_median_cv']:8.3f} "
                f"{-row['cv_reduction']:+8.3f}"
            )

    print("\nContrast vs control, session as replicate")
    header = (
        f"{'Group':10s} {'Cond.':10s} {'dates':>6s} {'fold':>6s} "
        f"{'95% CI':>16s} {'p':>9s} {'p_Holm':>9s} {'p_naive':>10s}"
    )
    print(header)
    print("-" * len(header))
    for row in branch.tests:
        interval = f"[{row.ci95_low:.3f},{row.ci95_high:.3f}]"
        print(
            f"{row.group[:10]:10s} {row.treatment[:10]:10s} {row.n_dates:6d} "
            f"{row.geo_fold:6.3f} {interval:>16s} {row.p_value:9.4f} "
            f"{row.p_holm:9.4f} {row.p_naive_pooled:10.2e}"
        )


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        result = run(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            control=args.control,
            anchor=args.anchor,
            outliers=args.outliers,
            make_plots=not args.no_plots,
            strict=args.strict,
            min_embryos=args.min_embryos,
        )
    except InputError as error:
        print(f"Input error: {error}", file=sys.stderr)
        return 2

    for warning in result.warnings:
        print(f"Warning: {warning}", file=sys.stderr)

    print(
        f"\n{len(result.tidy)} embryos read from "
        f"{len({m.source for m in result.tidy})} files "
        f"({len({(m.group, m.date) for m in result.tidy})} group x session groups), "
        f"anchor = {args.anchor} of {args.control}."
    )

    for label in ("keep", "drop"):
        if label in result.branches:
            _print_branch(label, result.branches[label])

    print(
        "\n'p_naive' treats each embryo as an independent replicate. It is shown\n"
        "only to expose the inflation from pseudoreplication. Report 'p_Holm'."
    )

    if args.mixed_model:
        reference = result.branches.get("drop") or next(iter(result.branches.values()))
        print("\n=== Mixed model: log2(I) ~ treatment + (1 | date) ===")
        for key, value in mixed_model(reference.normalized, control=args.control).items():
            if key == "available":
                continue
            print(f"\n--- {key} ---\n{value}")

    print(f"\nOutputs written to {args.output_dir}/:")
    for path in result.outputs:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
