import pytest
from raven_toolbox.biomass.scale import scale_biomass, sum_biomass

from overflow.biomass import GAM_COFACTORS, YEAST_BIOMASS
from overflow.config import BIO_RXN

pytestmark = pytest.mark.slow


def test_biomass_config_matches_the_model(conv_model):
    names = {r.name for r in conv_model.reactions}
    for component in YEAST_BIOMASS.components:
        assert component.pseudoreaction_name in names, component.name
    assert conv_model.metabolites.get_by_id(YEAST_BIOMASS.proton_met).name == "H+"


def test_components_sum_to_one_gram_per_gdw(conv_model):
    fractions = sum_biomass(conv_model, YEAST_BIOMASS)
    assert fractions["total"] == pytest.approx(1.0, abs=0.05)
    assert fractions["protein"] == pytest.approx(0.4, abs=0.15)


def test_scaling_protein_sets_the_requested_content(conv_model):
    model = conv_model.copy()
    scale_biomass(model, YEAST_BIOMASS, "protein", 0.27, balance_out="carbohydrate")
    fractions = sum_biomass(model, YEAST_BIOMASS)
    assert fractions["protein"] == pytest.approx(0.27, abs=1e-6)
    assert fractions["total"] == pytest.approx(1.0, abs=0.05)


def test_gam_cofactors_are_all_in_the_biomass_reaction(conv_model):
    present = {m.name for m in conv_model.reactions.get_by_id(BIO_RXN).metabolites}
    assert set(GAM_COFACTORS) <= present


def test_the_shipped_maintenance_is_read_back(conv_model):
    from overflow.biomass import current_gam

    assert current_gam(conv_model) == pytest.approx(55.3)


def test_keeping_the_model_maintenance_scales_only_the_composition(conv_model):
    from overflow.biomass import apply_protein_content, current_gam

    model = conv_model.copy()
    before = current_gam(model)
    composition = apply_protein_content(model, 0.27, gam="model")
    assert composition["protein"] == pytest.approx(0.27, abs=1e-6)
    assert composition["GAM"] == pytest.approx(before)
    assert current_gam(model) == pytest.approx(before)


def test_the_polymerization_policy_rebuilds_maintenance_from_the_composition(conv_model):
    from overflow.biomass import (
        GAM_NO_POLYMERIZATION,
        apply_protein_content,
        polymerization_cost,
    )

    model = conv_model.copy()
    composition = apply_protein_content(model, 0.27, gam="polymerization")
    expected = GAM_NO_POLYMERIZATION + polymerization_cost(model)
    assert composition["GAM"] == pytest.approx(expected)


def test_a_lower_protein_content_costs_less_to_polymerise(conv_model):
    from overflow.biomass import apply_protein_content

    lean = apply_protein_content(conv_model.copy(), 0.27, gam="polymerization")
    rich = apply_protein_content(conv_model.copy(), 0.55, gam="polymerization")
    assert lean["GAM"] < rich["GAM"]


def test_protein_scaling_can_be_turned_off(conv_model):
    from raven_toolbox.biomass.scale import sum_biomass

    from overflow.biomass import apply_protein_content

    model = conv_model.copy()
    before = sum_biomass(model, YEAST_BIOMASS)["protein"]
    composition = apply_protein_content(model, 0.27, scale_protein=False)
    assert composition["protein"] == pytest.approx(before)


def test_an_unknown_maintenance_policy_is_rejected(conv_model):
    from overflow.biomass import apply_protein_content

    with pytest.raises(ValueError, match="gam policy"):
        apply_protein_content(conv_model.copy(), 0.5, gam="guess")


def test_changing_maintenance_keeps_the_biomass_proton_offset(conv_model):
    """yeast-GEM's biomass releases three fewer protons than it
    hydrolyses ATP. Setting the proton to the GAM instead of shifting it
    adds three protons per gram of biomass, which is invisible in the
    stoichiometry and costs about a percent of growth."""
    from overflow.biomass import current_gam, set_growth_maintenance

    model = conv_model.copy()
    reaction = model.reactions.get_by_id(BIO_RXN)
    proton = next(m for m in reaction.metabolites if m.name == "H+")
    offset = reaction.metabolites[proton] - current_gam(model)

    set_growth_maintenance(model, current_gam(model) + 10.0)
    assert current_gam(model) == pytest.approx(65.3)
    assert reaction.metabolites[proton] - current_gam(model) == pytest.approx(offset)


def test_setting_the_maintenance_it_already_has_changes_nothing(conv_model):
    from overflow.biomass import current_gam, set_growth_maintenance

    model = conv_model.copy()
    reaction = model.reactions.get_by_id(BIO_RXN)
    before = {m.id: c for m, c in reaction.metabolites.items()}
    set_growth_maintenance(model, current_gam(model))
    after = {m.id: c for m, c in reaction.metabolites.items()}
    assert before == after
