"""The vendored model is an input, not a build product: these tests pin
the properties the rest of the pipeline relies on, so that swapping the
model file in cannot silently change them."""
import numpy as np
import pytest

from overflow.config import (
    BIO_RXN,
    BYPRODUCT_RXNS,
    C_SOURCE,
    CO2_RXN,
    NGAM_RXN,
    O2_RXN,
    OXPHOS_RXNS,
    POOL_RXN,
    PROTEIN_RXN,
)

pytestmark = pytest.mark.slow


def test_model_is_a_full_ecmodel(ec_model):
    assert ec_model.ec.gecko_light is False
    assert len(ec_model.ec.enzymes) == 1144
    assert len(ec_model.ec.rxns) == 4850


def test_every_reaction_the_pipeline_names_exists(ec_model):
    required = {
        BIO_RXN, NGAM_RXN, PROTEIN_RXN, C_SOURCE, CO2_RXN, O2_RXN, POOL_RXN,
        *BYPRODUCT_RXNS.values(),
    }
    missing = {r for r in required if r not in {x.id for x in ec_model.reactions}}
    assert not missing


def test_oxphos_reactions_exist_as_split_variants(ec_model):
    ids = {r.id for r in ec_model.reactions}
    for base in OXPHOS_RXNS:
        assert any(i == base or i.startswith(base + "_") for i in ids), base


def test_every_enzyme_has_a_usage_reaction_drawing_from_the_pool(ec_model):
    ids = {r.id for r in ec_model.reactions}
    missing = [e for e in ec_model.ec.enzymes if f"usage_prot_{e}" not in ids]
    assert not missing
    usage = ec_model.reactions.get_by_id(f"usage_prot_{ec_model.ec.enzymes[0]}")
    assert any(m.id == "prot_pool" for m in usage.metabolites)


def test_molecular_weights_are_daltons(ec_model):
    """Abundances are converted mmol/gDW -> mg/gDW by multiplying with
    ec.mw, which only works if mw is in Da (== mg/mmol)."""
    mw = np.asarray(ec_model.ec.mw, dtype=float)
    assert mw.min() > 1_000       # smallest yeast protein is a few kDa
    assert mw.max() < 1_000_000
    assert np.median(mw) == pytest.approx(55_000, rel=0.3)


def test_ngam_ships_pinned_and_must_be_freed_before_fitting(ec_model):
    ngam = ec_model.reactions.get_by_id(NGAM_RXN)
    assert ngam.bounds == (0.7, 0.7)


def test_growth_associated_maintenance_is_the_expected_magnitude(ec_model):
    """yeast-GEM folds polymerization into GAM; the original analysis used
    34 mmol ATP/gDW plus polymerization, which lands in the same place."""
    atp = next(m for m in ec_model.reactions.get_by_id(BIO_RXN).metabolites
               if m.name == "ATP")
    assert ec_model.reactions.get_by_id(BIO_RXN).metabolites[atp] == pytest.approx(-55.3)


def test_kcats_are_all_positive(ec_model):
    kcat = np.asarray(ec_model.ec.kcat, dtype=float)
    assert (kcat > 0).all()
