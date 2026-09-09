"""Figures for the analysis.

The main figure is a SuperPlot (Lord et al., J. Cell Biol. 2020). Individual
embryos are drawn as faint grey background points, and the mean of each
acquisition session is drawn on top as a large colour-coded marker. This shows
how many sessions support each box, and whether an effect holds across sessions
or comes from a single date. A boxplot over pooled embryos shows neither.

A second figure plots the raw control intensity per session, which is the
quantity the normalization removes. It is the fastest way to see how much
between-session drift the experiment carried.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless backend: the pipeline runs in a terminal or CI

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .normalize import Anchor, Normalized
from .stats import date_folds
from .tidy import DEFAULT_CONTROL


def treatment_order(rows: list[Normalized], control: str = DEFAULT_CONTROL) -> list[str]:
    """Plotting order: the control first, then conditions in order of appearance."""
    seen: list[str] = []
    for row in rows:
        if row.treatment not in seen:
            seen.append(row.treatment)
    if control in seen:
        seen.remove(control)
        seen.insert(0, control)
    return seen


def palette_for(order: list[str]) -> dict[str, tuple]:
    """Map treatments to colours by axis position, not alphabetically.

    Seaborn assigns a hue palette in category order, which would recolour a
    condition whenever another one is absent from a panel. Binding the colour to
    the axis position keeps it stable across panels, runs and subsets.
    """
    colours = sns.color_palette("Set2", n_colors=max(len(order), 8))
    return {name: colours[i] for i, name in enumerate(order)}


def _frames(rows: list[Normalized], control: str):
    frame = pd.DataFrame(
        {
            "GROUP": [r.group for r in rows],
            "DATE": [r.date for r in rows],
            "TREATMENT": [r.treatment for r in rows],
            "RATIO_NORM": [r.ratio_norm for r in rows],
        }
    )
    folds = date_folds(rows, control=control)
    day_frame = pd.DataFrame(
        {
            "GROUP": [f.group for f in folds],
            "DATE": [f.date for f in folds],
            "TREATMENT": [f.treatment for f in folds],
            "DATE_MEAN": [f.mean_norm for f in folds],
        }
    )
    return frame, day_frame


def _draw_panel(ax, frame, day_frame, order, title, control):
    sns.boxplot(
        data=frame, x="TREATMENT", y="RATIO_NORM", order=order, ax=ax,
        hue="TREATMENT", palette=palette_for(order), legend=False,
        fliersize=0, boxprops=dict(alpha=0.45), width=0.6,
    )
    # Individual embryos: background context, kept deliberately faint.
    sns.stripplot(
        data=frame, x="TREATMENT", y="RATIO_NORM", order=order, ax=ax,
        color="0.45", size=3.5, alpha=0.45, jitter=0.22,
    )
    # Per-session means: the true replicates of the experiment.
    sns.stripplot(
        data=day_frame, x="TREATMENT", y="DATE_MEAN", order=order, ax=ax,
        hue="DATE", palette="tab10", size=11, alpha=0.95,
        edgecolor="black", linewidth=0.9, jitter=0.12, dodge=False,
    )
    ax.axhline(1.0, color="crimson", linestyle="--", linewidth=1.4, zorder=0)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.set_xlabel("Condition", fontsize=11, labelpad=6)


def plot_groups(
    rows: list[Normalized],
    output: Path,
    control: str = DEFAULT_CONTROL,
    dpi: int = 300,
) -> Path:
    """Comparative figure with one panel per experiment group."""
    frame, day_frame = _frames(rows, control)
    groups = sorted(frame["GROUP"].unique())
    order = treatment_order(rows, control=control)

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(
        1, len(groups), figsize=(7 * len(groups), 6), sharey=True, squeeze=False
    )

    for ax, group in zip(axes[0], groups):
        subset = frame[frame["GROUP"] == group]
        day_subset = day_frame[day_frame["GROUP"] == group]
        _draw_panel(
            ax, subset, day_subset,
            [t for t in order if t in set(subset["TREATMENT"])],
            group, control,
        )
        ax.set_ylabel("")
        ax.legend_.remove() if ax.legend_ else None

    axes[0][0].set_ylabel(
        f"Normalized DCF intensity (vs same-session {control})", fontsize=11, labelpad=8
    )

    handles, labels = axes[0][-1].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles, labels, title="Acquisition date",
            loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
        )

    fig.suptitle(
        "Large points: per-session means (the replicates). "
        "Faint points: individual embryos.",
        fontsize=9, color="0.35", y=0.005, va="bottom",
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_control_drift(
    anchors: list[Anchor], output: Path, dpi: int = 300
) -> Path:
    """Raw control intensity per session: the drift the normalization removes."""
    frame = pd.DataFrame(
        {
            "GROUP": [a.group for a in anchors],
            "DATE": [a.date for a in anchors],
            "ANCHOR": [a.anchor for a in anchors],
            "N": [a.control_n for a in anchors],
        }
    ).sort_values(["GROUP", "DATE"])

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(max(7, 0.9 * len(frame)), 4.6))
    sns.lineplot(
        data=frame, x="DATE", y="ANCHOR", hue="GROUP", marker="o",
        ax=ax, linewidth=1.6, markersize=8,
    )
    ax.set_xlabel("Acquisition date", fontsize=11, labelpad=6)
    ax.set_ylabel("Raw control anchor", fontsize=11, labelpad=8)
    ax.set_title(
        "Control intensity per session, before normalization",
        fontsize=12, fontweight="bold", pad=10,
    )
    ax.tick_params(axis="x", rotation=45)
    if frame["GROUP"].nunique() < 2 and ax.legend_:
        ax.legend_.remove()

    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output
