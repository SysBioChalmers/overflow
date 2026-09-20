"""Sampling the conventional model under the measured rates.

The ATP accounting is arithmetic on a flux vector, so it is checked
against the published summary: fed the committed sampling means it must
reproduce the committed `selectedFluxes.txt`.
"""
import numpy as np
import pandas as pd
import pytest

from overflow.atp import atp_budget, cofactor_turnover
from overflow.config import CONDITION_ORDER, ROOT, load_conditions
from overflow.sampling import (
    FORMATE_DEHYDROGENASE,
    FORMATE_EXCHANGE,
    MEASURED_EXCHANGES,
    measured_rate,
    sampling_bounds,
)

LEGACY = ROOT / "legacy_matlab" / "results" / "randomSampling"

#: yeast-GEM 8.3.4 split glucose phosphorylation over two reactions.
LEGACY_GLYCOLYSIS_CONSUMING = ("r_4235", "r_0886", "r_0534")


@pytest.fixture(scope="module")
def legacy_means():
    table = pd.read_csv(LEGACY / "allFluxes.txt", sep="\t")
    return {
        condition: dict(zip(table["ID"], table[f"{condition}_AVERAGE"]))
        for condition in CONDITION_ORDER
    }


@pytest.fixture(scope="module")
def legacy_selected():
    return pd.read_csv(LEGACY / "selectedFluxes.txt", sep="\t").set_index("Row")


# --- the bounds -------------------------------------------------------

def test_uptakes_are_negative_and_secretions_positive(conditions):
    bounds = sampling_bounds(conditions["CN75"])
    low, high = bounds["r_1714"]
    assert low < high < 0, "glucose is taken up"
    low, high = bounds["r_1761"]
    assert 0 < low < high, "ethanol is secreted"


def test_every_measured_rate_is_held_within_five_percent(conditions):
    condition = conditions["CN75"]
    for reaction_id, (field_name, sign) in MEASURED_EXCHANGES.items():
        value = measured_rate(condition, field_name)
        if value == 0:
            continue
        low, high = bounds = sampling_bounds(condition)[reaction_id]
        assert abs(low) == pytest.approx(0.95 * value if sign > 0 else 1.05 * value)
        assert abs(high) == pytest.approx(1.05 * value if sign > 0 else 0.95 * value)


def test_an_undetected_byproduct_gets_a_small_allowance(conditions):
    """Blocking it would prevent the sampler from saying whether the
    model wants to make it, which is the question the run asks."""
    low, high = sampling_bounds(conditions["CN4"])["r_1761"]
    assert low == 0.0
    assert high == pytest.approx(0.01)


def test_formate_can_be_left_out_of_the_measurement_set(conditions):
    with_formate = sampling_bounds(conditions["CN22"], include_formate=True)
    without = sampling_bounds(conditions["CN22"], include_formate=False)
    assert FORMATE_EXCHANGE in with_formate
    assert FORMATE_EXCHANGE not in without
    assert set(with_formate) - {FORMATE_EXCHANGE} == set(without)


def test_growth_is_held_at_the_dilution_rate(conditions):
    low, high = sampling_bounds(conditions["hGR"])["r_2111"]
    assert low == pytest.approx(0.95 * 0.29)
    assert high == pytest.approx(1.05 * 0.29)


# --- the ATP accounting, against the published summary ---------------

@pytest.mark.parametrize("condition", CONDITION_ORDER)
def test_atp_budget_reproduces_the_published_summary(
    legacy_means, legacy_selected, condition
):
    means = legacy_means[condition]
    published = legacy_selected[condition]
    gam = float(published["GAEC_rATP"]) / means["r_4041"]

    ours = atp_budget(
        means,
        gam=gam,
        growth_rate=load_conditions()[condition].d_rate,
        glycolysis_consuming=LEGACY_GLYCOLYSIS_CONSUMING,
    )

    # The published table is written to four significant figures, so a
    # quantity derived from several of them carries their rounding.
    for row in (
        "rGlu", "ETC_rATP", "ETC_YATP_glu", "ETC_YATP_mu",
        "glycolysis_rATP", "glycolysis_YATP_glu", "glycolysis_YATP_mu",
        "GAEC_rATP", "GAEC_YATP_glu", "NGAM_rATP", "NGAM_YATP_glu",
        "Metabolism_rATP", "GAEC+NGAM+Metabolism_rATP",
        "rPDH", "rIDH", "rMDHc", "rMDHm", "rNDE",
    ):
        assert float(ours[row]) == pytest.approx(float(published[row]), rel=2e-3), row


def test_the_budget_balances(legacy_means, legacy_selected):
    """What respiration and glycolysis make is what growth, maintenance
    and the rest of metabolism spend."""
    means = legacy_means["CN4"]
    gam = float(legacy_selected.loc["GAEC_rATP", "CN4"]) / means["r_4041"]
    budget = atp_budget(means, gam=gam, growth_rate=0.1,
                        glycolysis_consuming=LEGACY_GLYCOLYSIS_CONSUMING)
    made = budget["ETC_rATP"] + budget["glycolysis_rATP"]
    spent = budget["GAEC+NGAM+Metabolism_rATP"]
    assert made == pytest.approx(spent)


def test_yields_are_rates_divided_by_uptake_and_growth(legacy_means, legacy_selected):
    means = legacy_means["CN38"]
    gam = float(legacy_selected.loc["GAEC_rATP", "CN38"]) / means["r_4041"]
    budget = atp_budget(means, gam=gam, growth_rate=0.1,
                        glycolysis_consuming=LEGACY_GLYCOLYSIS_CONSUMING)
    assert budget["ETC_YATP_glu"] == pytest.approx(budget["ETC_rATP"] / budget["rGlu"])
    assert budget["ETC_YATP_mu"] == pytest.approx(budget["ETC_rATP"] / 0.1)


def test_a_condition_with_no_glucose_does_not_divide_by_zero():
    budget = atp_budget({"r_1166": 0.0, "r_4041": 0.1}, gam=50.0, growth_rate=0.1)
    assert np.isnan(budget["ETC_YATP_glu"])
    assert np.isfinite(budget["ETC_YATP_mu"])


# --- cofactor turnover ------------------------------------------------

@pytest.mark.slow
def test_cofactor_turnover_counts_production_only(conv_model, legacy_means):
    turnover = cofactor_turnover(conv_model, legacy_means["CN4"], "NAD", "rNAD")
    assert (turnover >= 0).all()
    assert turnover["rNAD[tot]"] == pytest.approx(
        turnover.drop("rNAD[tot]").sum()
    )


@pytest.mark.slow
def test_cofactor_turnover_is_reported_per_compartment(conv_model, legacy_means):
    turnover = cofactor_turnover(conv_model, legacy_means["hGR"], "NAD", "rNAD")
    compartments = {k for k in turnover.index if k != "rNAD[tot]"}
    assert len(compartments) >= 2, "NAD turns over in more than one compartment"
    assert all(k.startswith("rNAD[") for k in compartments)


# --- the loop-free screening cache -----------------------------------

def test_the_loop_free_screening_round_trips(tmp_path):
    """Screening every reaction for loop involvement is the dominant
    cost of a sampling run and depends only on the network, so it is
    kept rather than repeated."""
    from overflow.run_sampling import load_good_reactions, save_good_reactions

    assert load_good_reactions(tmp_path, "free") is None
    save_good_reactions(tmp_path, "free", ["r_2", "r_1"])
    assert load_good_reactions(tmp_path, "free") == ["r_1", "r_2"]
    assert load_good_reactions(tmp_path, "full") is None


def test_no_cache_directory_means_no_caching(tmp_path):
    from overflow.run_sampling import load_good_reactions, save_good_reactions

    save_good_reactions(None, "free", ["r_1"])
    assert load_good_reactions(None, "free") is None
    save_good_reactions(tmp_path, "free", [])
    assert load_good_reactions(tmp_path, "free") is None


# --- holding reactions to their loop-free range ----------------------

def test_loopless_bounds_tighten_only_what_a_cycle_inflates():
    """A reaction whose ordinary range reaches the model's arbitrary
    bound, but whose loop-free range is a few units, is held to the
    latter; a reaction already inside its loop-free range is left be."""
    import cobra
    import pandas as pd

    from overflow.sampling import apply_loopless_bounds

    model = cobra.Model("m")
    a = cobra.Reaction("a", lower_bound=-1000.0, upper_bound=1000.0)
    b = cobra.Reaction("b", lower_bound=0.0, upper_bound=3.0)
    model.add_reactions([a, b])

    ranges = pd.DataFrame(
        {"minimum": [-2.079, 0.0], "maximum": [2.908, 5.0]}, index=["a", "b"]
    )
    tightened = apply_loopless_bounds(model, ranges)

    assert tightened == 1
    assert a.lower_bound == pytest.approx(-2.079, abs=1e-6)
    assert a.upper_bound == pytest.approx(2.908, abs=1e-6)
    assert b.bounds == (0.0, 3.0), "already inside its loop-free range"


def test_loopless_bounds_never_widen_a_reaction():
    """The loop-free range of a reaction can exceed a bound the
    condition imposes; the condition wins."""
    import cobra
    import pandas as pd

    from overflow.sampling import apply_loopless_bounds

    model = cobra.Model("m")
    r = cobra.Reaction("r", lower_bound=-1.0, upper_bound=1.0)
    model.add_reactions([r])
    apply_loopless_bounds(
        model, pd.DataFrame({"minimum": [-50.0], "maximum": [50.0]}, index=["r"])
    )
    assert r.bounds == (-1.0, 1.0)


def test_a_crossed_range_leaves_the_reaction_alone():
    """Numerical noise can put a loop-free minimum above its maximum;
    that must not produce an infeasible reaction."""
    import cobra
    import pandas as pd

    from overflow.sampling import apply_loopless_bounds

    model = cobra.Model("m")
    r = cobra.Reaction("r", lower_bound=0.0, upper_bound=1.0)
    model.add_reactions([r])
    apply_loopless_bounds(
        model, pd.DataFrame({"minimum": [5.0], "maximum": [2.0]}, index=["r"])
    )
    assert r.bounds == (0.0, 1.0)
