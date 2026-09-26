"""Figures summarising enzyme capacity usage and the ribosomal subunit abundances.

They reproduce the published figures. The usage figure follows
``legacy_matlab/code/boxplotEnzymeUsage.R`` (ggplot2, ``theme_classic``, 7 pt
text, 10 x 4.5 cm, one colour per panel), and the subunit figure is the
default MATLAB density plot ``ribosome.m`` saved. The numbers behind the usage
panels are written beside the figure as a table.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

from overflow.config import CONDITION_ORDER  # noqa: E402

#: Pathways shown in the main figure and in the supplement.
SELECTED_SYSTEMS = ("Glycolysis", "TCA cycle", "ETC", "Ribosome")
SUPPLEMENT_SYSTEMS = (
    "PP shunt",
    "THF cycle",
    "Nitrogen metabolism",
    "Amino acid metabolism",
)

#: Panel colours of the published figures, by panel position.
LEGACY_COLOURS = ("#CBBBA0", "#1D1D1B", "#1D71B8", "#878787")

#: Helvetica, or the nearest metric-compatible face that is installed.
_INSTALLED = {face.name for face in font_manager.fontManager.ttflist}
FONTS = [
    name for name in ("Helvetica", "Arial", "Nimbus Sans", "Liberation Sans")
    if name in _INSTALLED
] or ["DejaVu Sans"]

# ggplot2 sizes. A line of ``size`` s is drawn s * .pt * 0.75 points wide, and
# ``text = element_text(size = 7)`` makes axis and strip text 0.8 of that.
_POINTS_PER_SIZE = 72.27 / 25.4 * 0.75
_TEXT = 7.0
_SMALL_TEXT = 0.8 * _TEXT
_AXIS_LINE = 0.15 * _POINTS_PER_SIZE
_BOX_LINE = 0.35 * _POINTS_PER_SIZE
_TICK_LENGTH = 2.75
_PANEL_GAP = 5.5
_OUTLIER = 4.4
_GREY10, _GREY20, _GREY30 = "#1A1A1A", "#333333", "#4D4D4D"


def capacity_usage_figure(
    capacity: pd.DataFrame,
    systems: Sequence[str] = SELECTED_SYSTEMS,
    path: Optional[Path | str] = None,
    conditions: Optional[Sequence[str]] = None,
    width: float = 10 / 2.54,
    height: float = 4.5 / 2.54,
):
    """One panel per pathway, one unfilled box per condition.

    ``capacity`` is the table from
    :func:`overflow.usage.capacity_usage_by_system`.
    """
    if conditions is None:
        conditions = [c for c in CONDITION_ORDER if c in capacity.columns]

    missing = [s for s in systems if s not in set(capacity["GOterm"])]
    if missing:
        raise ValueError(f"no enzymes annotated to {missing}")

    with plt.rc_context({"font.family": FONTS, "pdf.fonttype": 42}):
        figure = plt.figure(figsize=(width, height), facecolor="white")
        left, right, bottom, top = 0.105, 0.985, 0.30, 0.855
        panel = (right - left) * width * 72 - _PANEL_GAP * (len(systems) - 1)
        panel /= len(systems)
        axes = figure.subplots(
            1, len(systems), sharey=True,
            gridspec_kw=dict(left=left, right=right, bottom=bottom, top=top,
                             wspace=_PANEL_GAP / panel),
        )
        axes = [axes] if len(systems) == 1 else list(axes)

        for index, (axis, system) in enumerate(zip(axes, systems)):
            rows = capacity[capacity["GOterm"] == system]
            values = [rows[c].dropna().to_numpy() for c in conditions]
            colour = LEGACY_COLOURS[index % len(LEGACY_COLOURS)]

            axis.boxplot(
                values,
                widths=0.75,
                whis=1.5,
                showcaps=False,
                showfliers=True,
                patch_artist=True,
                boxprops=dict(facecolor="white", edgecolor=colour, linewidth=_BOX_LINE),
                medianprops=dict(color=colour, linewidth=2 * _BOX_LINE,
                                 solid_capstyle="butt"),
                whiskerprops=dict(color=colour, linewidth=_BOX_LINE),
                flierprops=dict(marker="o", markersize=_OUTLIER, linestyle="none",
                                markerfacecolor=colour, markeredgecolor=colour,
                                markeredgewidth=0),
            )

            axis.set_title(system, fontsize=_SMALL_TEXT, color=_GREY10, pad=4.4)
            axis.set_xticks(range(1, len(conditions) + 1))
            axis.set_xticklabels(conditions, rotation=90, fontsize=_SMALL_TEXT, color=_GREY30)
            axis.tick_params(axis="x", length=_TICK_LENGTH, width=_AXIS_LINE,
                             color=_GREY20, pad=2.2)
            axis.set_xlim(0.4, len(conditions) + 0.6)
            axis.set_facecolor("white")
            for side in ("top", "right"):
                axis.spines[side].set_visible(False)
            axis.spines["bottom"].set(color="black", linewidth=_AXIS_LINE)
            if index == 0:
                axis.spines["left"].set(color="black", linewidth=_AXIS_LINE)
                axis.set_yticks([0, 25, 50, 75, 100])
                axis.tick_params(axis="y", length=_TICK_LENGTH, width=_AXIS_LINE,
                                 color=_GREY20, pad=2.2, labelsize=_SMALL_TEXT,
                                 labelcolor=_GREY30)
                axis.set_ylabel("Capacity usage (%)", fontsize=_TEXT, color="black",
                                labelpad=2.75)
            else:
                axis.spines["left"].set_visible(False)
                axis.tick_params(axis="y", left=False, labelleft=False)

        axes[0].set_ylim(-5, 105)

        if path is not None:
            figure.savefig(path, facecolor="white")
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


# MATLAB's default axes, as ``saveas`` printed them: a 341 x 247 pt box, 8 pt
# tick labels, 8.7 pt labels and bold title, and 0.4 pt lines.
_MATLAB_BLUE = "#0072BD"
_MATLAB_INK = "#262626"
_MATLAB_LINE = 0.4


def subunit_abundance_figure(
    means: pd.Series,
    path: Optional[Path | str] = None,
):
    """Distribution of the average abundance of the candidate ribosomal subunits."""
    positive = means[means > 0]
    grid, density = subunit_density(np.log10(positive.to_numpy()))

    axes_width, axes_height = 341.0, 247.0
    margin = dict(left=46.0, right=14.0, bottom=40.0, top=26.0)
    width = axes_width + margin["left"] + margin["right"]
    height = axes_height + margin["bottom"] + margin["top"]

    with plt.rc_context({"font.family": FONTS, "pdf.fonttype": 42}):
        figure = plt.figure(figsize=(width / 72, height / 72), facecolor="white")
        axis = figure.add_axes([
            margin["left"] / width, margin["bottom"] / height,
            axes_width / width, axes_height / height,
        ])
        axis.plot(grid, density, color=_MATLAB_BLUE, linewidth=_MATLAB_LINE)

        axis.set_xlim(np.floor(grid.min() * 2) / 2, np.ceil(grid.max() * 2) / 2)
        axis.set_ylim(0, np.ceil(density.max() * 10) / 10)
        axis.xaxis.set_major_locator(MultipleLocator(0.5))
        axis.yaxis.set_major_locator(MultipleLocator(0.1))
        axis.tick_params(direction="in", top=True, right=True, length=3.4,
                         width=_MATLAB_LINE, labelsize=8, color=_MATLAB_INK,
                         labelcolor=_MATLAB_INK)
        for spine in axis.spines.values():
            spine.set(color=_MATLAB_INK, linewidth=_MATLAB_LINE)
        axis.set_xlabel("Subunit abundance (log10(mmol/gDCW))", fontsize=8.7,
                        color=_MATLAB_INK, labelpad=3)
        axis.set_ylabel("Density", fontsize=8.7, color=_MATLAB_INK, labelpad=3)
        axis.set_title("Distribution of average ribosomal subunit abundances",
                       fontsize=8.7, fontweight="bold", color=_MATLAB_INK, pad=4)

        if path is not None:
            figure.savefig(path, facecolor="white")
            plt.close(figure)
    return figure
