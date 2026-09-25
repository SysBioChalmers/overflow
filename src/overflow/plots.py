"""Figures summarising enzyme capacity usage.

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
