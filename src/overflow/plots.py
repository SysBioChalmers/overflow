"""Figures summarising enzyme capacity usage and the ribosomal subunit abundances.

Each panel is one pathway and each box one condition, so a pathway that
tightens as the carbon-to-nitrogen ratio falls shows up as boxes walking
upward. Colour repeats the panel title rather than carrying identity of
its own; the numbers behind every panel are written beside the figure as
a table.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from overflow.config import CONDITION_ORDER  # noqa: E402

#: Pathways shown in the main figure and in the supplement.
SELECTED_SYSTEMS = ("Glycolysis", "TCA cycle", "ETC", "Ribosome")
SUPPLEMENT_SYSTEMS = (
    "PP shunt",
    "THF cycle",
    "Nitrogen metabolism",
    "Amino acid metabolism",
)

#: Categorical slots one to four, validated against both surfaces.
PALETTE = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")

INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e6e5e0"
SURFACE = "#fcfcfb"


def capacity_usage_figure(
    capacity: pd.DataFrame,
    systems: Sequence[str] = SELECTED_SYSTEMS,
    path: Optional[Path | str] = None,
    conditions: Optional[Sequence[str]] = None,
    width: float = 7.2,
    height: float = 2.6,
):
    """One panel per pathway, one box per condition.

    ``capacity`` is the table from
    :func:`overflow.usage.capacity_usage_by_system`.
    """
    if conditions is None:
        conditions = [c for c in CONDITION_ORDER if c in capacity.columns]

    missing = [s for s in systems if s not in set(capacity["GOterm"])]
    if missing:
        raise ValueError(f"no enzymes annotated to {missing}")

    figure, axes = plt.subplots(
        1, len(systems), figsize=(width, height), sharey=True,
        facecolor=SURFACE, layout="constrained",
    )
    axes = [axes] if len(systems) == 1 else list(axes)

    for index, (axis, system) in enumerate(zip(axes, systems)):
        rows = capacity[capacity["GOterm"] == system]
        values = [rows[c].dropna().to_numpy() for c in conditions]
        colour = PALETTE[index % len(PALETTE)]

        drawn = axis.boxplot(
            values,
            widths=0.55,
            showfliers=True,
            patch_artist=True,
            medianprops=dict(color=INK, linewidth=1.1),
            flierprops=dict(
                marker="o", markersize=2.2, markerfacecolor="none",
                markeredgecolor=colour, markeredgewidth=0.6,
            ),
        )
        for box in drawn["boxes"]:
            box.set(facecolor=colour, alpha=0.18, edgecolor=colour, linewidth=0.9)
        for part in ("whiskers", "caps"):
            for line in drawn[part]:
                line.set(color=colour, linewidth=0.9)

        axis.set_title(system, fontsize=8, color=INK, pad=6)
        axis.set_xticks(range(1, len(conditions) + 1))
        axis.set_xticklabels(conditions, rotation=90, fontsize=7, color=MUTED)
        axis.tick_params(axis="y", labelsize=7, colors=MUTED, length=3)
        axis.tick_params(axis="x", length=0)
        axis.set_facecolor(SURFACE)
        axis.set_axisbelow(True)
        axis.yaxis.grid(True, color=GRID, linewidth=0.5)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            axis.spines[side].set(color=GRID, linewidth=0.6)
        if index == 0:
            axis.set_ylabel("Capacity usage (%)", fontsize=8, color=INK)

    axes[0].set_ylim(-4, 104)

    if path is not None:
        figure.savefig(path, facecolor=SURFACE)
        plt.close(figure)
    return figure


def subunit_density(
    log_abundance: Sequence[float], bandwidth: float = 0.1, points: int = 100
):
    """Gaussian kernel density of log10 abundances on a grid.

    The grid runs three bandwidths past the data on each side and the
    kernel bandwidth is absolute, in log10 units, as in the published
    analysis.
    """
    values = np.asarray(log_abundance, float)
    grid = np.linspace(values.min() - 3 * bandwidth, values.max() + 3 * bandwidth, points)
    scaled = (grid[:, None] - values[None, :]) / bandwidth
    density = np.exp(-0.5 * scaled**2).sum(axis=1) / (len(values) * bandwidth * np.sqrt(2 * np.pi))
    return grid, density


def subunit_abundance_figure(
    means: pd.Series,
    threshold: float,
    path: Optional[Path | str] = None,
    width: float = 4.0,
    height: float = 2.6,
):
    """Distribution of the average abundance of the candidate ribosomal subunits.

    A dashed line marks the abundance below which a subunit is not part of
    the core ribosome; the count either side of it is written on the plot.
    """
    positive = means[means > 0]
    grid, density = subunit_density(np.log10(positive.to_numpy()))
    cut = np.log10(threshold)
    kept = int((positive >= threshold).sum())

    figure, axis = plt.subplots(figsize=(width, height), facecolor=SURFACE)
    axis.set_facecolor(SURFACE)
    axis.plot(grid, density, color=PALETTE[0], linewidth=1.6)
    axis.axvline(cut, color=MUTED, linewidth=0.8, linestyle="--")
    axis.set_xlim(grid.min(), grid.max() + 1.2)
    axis.text(
        cut + 0.08, axis.get_ylim()[1] * 0.96,
        f"core ribosome:\n{kept} of {len(positive)} subunits\nat or above {threshold:g}",
        fontsize=7, color=INK, va="top", ha="left",
    )
    axis.set_xlabel("Average subunit abundance, log10 (mmol/gDW)", fontsize=8, color=INK)
    axis.set_ylabel("Density", fontsize=8, color=INK)
    axis.set_title("Average ribosomal subunit abundance", fontsize=9, color=INK, loc="left")
    axis.tick_params(labelsize=7, colors=MUTED, length=2)
    axis.yaxis.grid(True, color=GRID, linewidth=0.5)
    axis.set_axisbelow(True)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set(color=GRID, linewidth=0.6)
    figure.tight_layout()

    if path is not None:
        figure.savefig(path, facecolor=SURFACE)
        plt.close(figure)
    return figure
