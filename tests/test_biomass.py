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
