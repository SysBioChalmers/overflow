"""Abundance reconciliation, on a model whose arithmetic is exact.

In the tiny model a unit of flux costs 10 mg of E1 or 5 mg of E2, so
every quantity below can be worked out by hand.
"""
import numpy as np
import pytest
from geckopy import constrain_enz_concs

from overflow.flexibilize import (
    apply_concentrations,
    flexibilize_proteins,
    measured_enzymes,
    minimum_usage_pass,
    write_back_concentrations,
)
from tiny import tiny_ec_model


@pytest.fixture
def model():
    """Pool of 75 mg forces both routes to carry flux at full supply.

    Ten units of product cost 10x + 5(10 - x) mg with x units on E1, so
    the pool allows at most 5 units on E1 and demands at least 5 on E2:
    E2 must supply at least 25 mg of enzyme.
    """
    return tiny_ec_model(pool=75.0, supply=10.0)


def _set_concentrations(model, **values):
    concentrations = np.full(len(model.ec.enzymes), np.nan)
    for enzyme, value in values.items():
        concentrations[model.ec.enzymes.index(enzyme)] = value
    model.ec.concs = concentrations


# --- minimum usage pass ----------------------------------------------

def test_an_abundance_below_what_the_model_needs_is_raised(model):
    _set_concentrations(model, E2=10.0)
    result = minimum_usage_pass(model, target_growth=10.0, bio_rxn="BIO")
    assert result.raised == ["E2"]
    assert result.minimum[0] == pytest.approx(1.01 * 25.0)
    assert result.concentrations[model.ec.enzymes.index("E2")] == pytest.approx(25.25)


def test_a_sufficient_abundance_is_left_alone(model):
    _set_concentrations(model, E2=40.0)
    result = minimum_usage_pass(model, target_growth=10.0, bio_rxn="BIO")
    assert result.raised == []
    assert result.concentrations[model.ec.enzymes.index("E2")] == 40.0


def test_unmeasured_enzymes_are_not_considered(model):
    _set_concentrations(model, E2=10.0)
    assert measured_enzymes(model) == [model.ec.enzymes.index("E2")]
    result = minimum_usage_pass(model, target_growth=10.0, bio_rxn="BIO")
    assert np.isnan(result.concentrations[model.ec.enzymes.index("E1")])


def test_the_pass_leaves_the_growth_bounds_as_it_found_them(model):
    _set_concentrations(model, E2=10.0)
    before = model.reactions.get_by_id("BIO").bounds
    minimum_usage_pass(model, target_growth=10.0, bio_rxn="BIO")
    assert model.reactions.get_by_id("BIO").bounds == before


def test_the_pass_does_not_constrain_the_model_it_measures(model):
    """Each minimum must be what one enzyme needs on its own, so the
    pass must not leave earlier enzymes capped behind it."""
    _set_concentrations(model, E1=1.0, E2=1.0)
    minimum_usage_pass(model, target_growth=10.0, bio_rxn="BIO")
    for enzyme in model.ec.enzymes:
        assert model.reactions.get_by_id(f"usage_prot_{enzyme}").upper_bound == 1000.0


# --- iterative flexibilization ---------------------------------------

def test_a_limiting_cap_is_released_and_growth_reaches_the_target(model):
    _set_concentrations(model, E1=5.0, E2=5.0)
    apply_concentrations(model, model.ec.concs)
    assert model.optimize().objective_value < 10.0

    result = flexibilize_proteins(model, target_growth=10.0, bio_rxn="BIO")
    assert result.reached_target
    assert result.growth == pytest.approx(10.0)
    assert set(result.released) <= {"E1", "E2"}


def test_a_released_enzyme_is_pulled_back_to_what_it_uses(model):
    """Releasing a cap should cost only what the model spends, not the
    1000 mg the release temporarily allows."""
    _set_concentrations(model, E1=5.0, E2=5.0)
    apply_concentrations(model, model.ec.concs)
    result = flexibilize_proteins(model, target_growth=10.0, bio_rxn="BIO")

    for enzyme in result.released:
        bound = model.reactions.get_by_id(f"usage_prot_{enzyme}").upper_bound
        assert bound < 1000.0
    table = result.table()
    assert (table["modified_mg_gDW"] >= table["previous_mg_gDW"]).all()


def test_the_pool_is_grown_when_no_single_enzyme_is_limiting():
    """With generous caps but too small a pool, there is nothing to
    release, so the budget itself has to give."""
    model = tiny_ec_model(pool=10.0, supply=10.0)
    _set_concentrations(model, E1=1000.0, E2=1000.0)
    apply_concentrations(model, model.ec.concs)

    result = flexibilize_proteins(model, target_growth=10.0, bio_rxn="BIO")
    assert result.pool_increases > 0
    assert result.pool_after > result.pool_before
    assert result.reached_target


def test_an_unreachable_target_is_reported_rather_than_looped_forever():
    """Supply, not protein, caps this one: no amount of enzyme helps."""
    model = tiny_ec_model(pool=1000.0, supply=1.0)
    _set_concentrations(model, E1=1.0, E2=1.0)
    apply_concentrations(model, model.ec.concs)

    result = flexibilize_proteins(model, target_growth=10.0, max_iterations=20, bio_rxn="BIO")
    assert not result.reached_target
    assert result.growth == pytest.approx(1.0)


# --- keeping the record straight -------------------------------------

def test_concentrations_are_written_back_from_the_caps(model):
    _set_concentrations(model, E1=5.0, E2=5.0)
    apply_concentrations(model, model.ec.concs)
    model.reactions.get_by_id("usage_prot_E1").upper_bound = 42.0

    write_back_concentrations(model)
    assert model.ec.concs[model.ec.enzymes.index("E1")] == 42.0
    assert model.ec.concs[model.ec.enzymes.index("E2")] == 5.0


def test_applying_concentrations_caps_the_usage_reactions(model):
    _set_concentrations(model, E2=7.0)
    apply_concentrations(model, model.ec.concs)
    assert model.reactions.get_by_id("usage_prot_E2").upper_bound == 7.0
    assert model.reactions.get_by_id("usage_prot_E1").upper_bound == 1000.0


# --- robustness of finding the limiting cap ---------------------------

def test_the_limiting_cap_is_found_even_when_shadow_prices_are_silent(model, monkeypatch):
    """On a degenerate optimum the duals can all come back zero. The
    probe has to find the limiting cap anyway, or the search concludes
    that nothing is limiting and grows the pool forever."""
    import numpy as np

    import overflow.flexibilize as flexibilize

    def silent(model, proteins=None, **kwargs):
        n = len(proteins if proteins is not None else model.ec.enzymes)
        return np.ones(n, dtype=bool), np.zeros(n, dtype=float)

    monkeypatch.setattr(flexibilize, "get_conc_control_coeffs", silent)

    _set_concentrations(model, E1=5.0, E2=5.0)
    apply_concentrations(model, model.ec.concs)
    result = flexibilize_proteins(model, target_growth=10.0, bio_rxn="BIO")

    assert result.reached_target
    assert result.released, "no cap was released despite one being limiting"
    assert result.pool_increases == 0


def test_binding_caps_ignores_enzymes_that_are_already_free(model):
    from overflow.flexibilize import binding_caps

    _set_concentrations(model, E1=5.0, E2=5.0)
    apply_concentrations(model, model.ec.concs)
    solution = model.optimize()
    assert set(binding_caps(model, solution, 1e-6)) <= {"E1", "E2"}

    model.reactions.get_by_id("usage_prot_E2").upper_bound = 1000.0
    solution = model.optimize()
    assert "E2" not in binding_caps(model, solution, 1e-6)


def test_a_stalled_search_gives_up_instead_of_inflating_the_pool(model, monkeypatch):
    """Growing the protein budget without ever gaining growth means the
    budget is not the constraint; the run should stop and say so rather
    than leave behind a model with an invented pool."""
    import overflow.flexibilize as flexibilize

    monkeypatch.setattr(flexibilize, "_most_limiting", lambda *a, **k: None)

    _set_concentrations(model, E1=5.0, E2=5.0)
    apply_concentrations(model, model.ec.concs)
    result = flexibilize_proteins(
        model, target_growth=10.0, bio_rxn="BIO", patience=5, max_iterations=500
    )
    assert not result.reached_target
    assert result.pool_increases <= 7
