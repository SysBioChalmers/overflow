"""Summarising enzyme usage.

The transformation from per-enzyme usage to capacity usage per system is
pure arithmetic on a table, so it can be checked against the published
output exactly: feed it the committed `enzymeUsages.txt` and it must
reproduce the committed `capUsage.txt`.
"""
import numpy as np
import pandas as pd
import pytest

from overflow.config import ROOT
from overflow.usage import (
    capacity_usage_by_system,
    combine_usage,
    read_annotation,
    usage_by_system,
)

LEGACY = ROOT / "legacy_matlab" / "results" / "enzymeUsage"


@pytest.fixture(scope="module")
def legacy_usage():
    return pd.read_csv(LEGACY / "enzymeUsages.txt", sep="\t")


@pytest.fixture(scope="module")
def legacy_capacity():
    return pd.read_csv(LEGACY / "capUsage.txt", sep="\t")


@pytest.fixture(scope="module")
def annotation():
    return read_annotation()


def test_capacity_usage_reproduces_the_published_table(
    legacy_usage, legacy_capacity, annotation
):
    ours = capacity_usage_by_system(legacy_usage, annotation)
    assert list(ours.columns) == list(legacy_capacity.columns)
    assert list(ours["protID"]) == list(legacy_capacity["protID"])
    assert list(ours["GOterm"]) == list(legacy_capacity["GOterm"])
    for condition in ("CN4", "CN22", "CN38", "CN75", "hGR"):
        assert ours[condition].to_numpy() == pytest.approx(
            legacy_capacity[condition].to_numpy()
        )


def test_enzymes_used_in_no_condition_are_dropped(annotation):
    usage = pd.DataFrame(
        {
            "protID": ["P00127", "P00128"],
            "geneID": ["YFR033C", "YDR529C"],
            "protName": ["QCR6", "QCR7"],
            "capUse_CN4": [0.0, 0.5],
            "capUse_CN22": [0.0, 0.25],
        }
    )
    result = capacity_usage_by_system(usage, annotation)
    assert list(result["protID"]) == ["P00128"]


def test_capacity_usage_is_reported_as_a_percentage(annotation):
    usage = pd.DataFrame(
        {
            "protID": ["P00127"], "geneID": ["YFR033C"], "protName": ["QCR6"],
            "capUse_CN4": [0.4213],
        }
    )
    assert capacity_usage_by_system(usage, annotation)["CN4"].iloc[0] == pytest.approx(42.13)


def test_an_unannotated_protein_is_left_out(annotation):
    usage = pd.DataFrame(
        {
            "protID": ["NOT_ANNOTATED"], "geneID": ["YXX000W"], "protName": ["X"],
            "capUse_CN4": [0.5],
        }
    )
    assert capacity_usage_by_system(usage, annotation).empty


def test_a_protein_in_two_systems_keeps_the_first(annotation):
    """ADE3 is listed under both amino acid metabolism and the THF cycle."""
    duplicated = annotation[annotation["protID"] == "P07245"]
    assert len(duplicated) == 2
    usage = pd.DataFrame(
        {
            "protID": ["P07245"], "geneID": ["YGR204W"], "protName": ["ADE3"],
            "capUse_CN4": [0.5],
        }
    )
    result = capacity_usage_by_system(usage, annotation)
    assert result["GOterm"].iloc[0] == duplicated["system"].iloc[0]


def test_combining_conditions_keeps_the_published_column_layout():
    def one(cap, abs_, ub):
        return pd.DataFrame(
            {
                "protID": ["P1"], "geneID": ["Y1"], "protName": ["N1"],
                "capUse": [cap], "absUse": [abs_], "UB": [ub],
            }
        )

    wide = combine_usage({"CN4": one(0.5, 1.0, 2.0), "hGR": one(0.25, 1.0, 4.0)})
    assert list(wide.columns) == [
        "protID", "geneID", "protName",
        "capUse_CN4", "capUse_hGR", "absUse_CN4", "absUse_hGR", "UB_CN4", "UB_hGR",
    ]
    assert wide["capUse_hGR"].iloc[0] == 0.25


def test_combining_conditions_matches_proteins_by_identifier():
    """The conditions need not list their enzymes in the same order."""
    first = pd.DataFrame(
        {"protID": ["P1", "P2"], "geneID": ["Y1", "Y2"], "protName": ["A", "B"],
         "capUse": [0.1, 0.2], "absUse": [1.0, 2.0], "UB": [10.0, 10.0]}
    )
    second = first.iloc[::-1].reset_index(drop=True).copy()
    second["capUse"] = [0.9, 0.8]
    wide = combine_usage({"CN4": first, "hGR": second})
    row = wide.set_index("protID").loc["P1"]
    assert row["capUse_CN4"] == 0.1 and row["capUse_hGR"] == 0.8


def test_system_medians_summarise_the_table(legacy_usage, annotation):
    capacity = capacity_usage_by_system(legacy_usage, annotation)
    medians = usage_by_system(capacity)
    assert set(medians["GOterm"]) == set(capacity["GOterm"])
    assert (medians[["CN4", "hGR"]] <= 100.0).all().all()


# --- figures ----------------------------------------------------------

def test_the_figure_draws_one_panel_per_system(legacy_usage, annotation, tmp_path):
    from overflow.plots import SELECTED_SYSTEMS, capacity_usage_figure

    capacity = capacity_usage_by_system(legacy_usage, annotation)
    figure = capacity_usage_figure(capacity, SELECTED_SYSTEMS, tmp_path / "f.pdf")
    assert [axis.get_title() for axis in figure.axes] == list(SELECTED_SYSTEMS)
    assert (tmp_path / "f.pdf").stat().st_size > 1000


def test_every_panel_actually_has_boxes(legacy_usage, annotation, tmp_path):
    """A panel whose system matched nothing would still draw, empty."""
    from overflow.plots import SELECTED_SYSTEMS, capacity_usage_figure

    capacity = capacity_usage_by_system(legacy_usage, annotation)
    figure = capacity_usage_figure(capacity, SELECTED_SYSTEMS, tmp_path / "f.pdf")
    for axis, system in zip(figure.axes, SELECTED_SYSTEMS):
        assert len(axis.patches) == 5, f"{system} should have one box per condition"


def test_a_system_with_no_enzymes_is_refused(legacy_usage, annotation, tmp_path):
    from overflow.plots import capacity_usage_figure

    capacity = capacity_usage_by_system(legacy_usage, annotation)
    with pytest.raises(ValueError, match="no enzymes annotated"):
        capacity_usage_figure(capacity, ["Photosynthesis"], tmp_path / "f.pdf")


def test_the_axis_covers_the_full_percentage_range(legacy_usage, annotation, tmp_path):
    """Capacity usage is a fraction of what was available, so panels have
    to share a 0-100 axis or they cannot be compared."""
    from overflow.plots import SELECTED_SYSTEMS, capacity_usage_figure

    capacity = capacity_usage_by_system(legacy_usage, annotation)
    figure = capacity_usage_figure(capacity, SELECTED_SYSTEMS, tmp_path / "f.pdf")
    for axis in figure.axes:
        low, high = axis.get_ylim()
        assert low <= 0 and high >= 100
