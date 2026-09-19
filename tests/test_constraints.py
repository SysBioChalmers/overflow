import pytest

from overflow.config import BIO_RXN, BYPRODUCT_RXNS, C_SOURCE, NGAM_RXN, POOL_RXN
from overflow.constraints import (
    constrain_byproducts,
    constrain_uptake,
    free_ngam,
    set_chemostat_constraints,
)

pytestmark = pytest.mark.slow


def test_undetected_byproducts_are_blocked(ec_model, conditions):
    with ec_model as model:
        constrain_byproducts(model, conditions["CN4"])
        for reaction_id in BYPRODUCT_RXNS.values():
            assert model.reactions.get_by_id(reaction_id).bounds == (0.0, 0.0)


def test_detected_byproducts_keep_ten_percent_headroom(ec_model, conditions):
    with ec_model as model:
        constrain_byproducts(model, conditions["CN75"])
        ethanol = model.reactions.get_by_id(BYPRODUCT_RXNS["Ethanol"])
        assert ethanol.upper_bound == pytest.approx(1.1 * 2.282585358)
        assert ethanol.lower_bound == 0.0


def test_uptake_is_capped_five_percent_above_the_measurement(ec_model, conditions):
    with ec_model as model:
        constrain_uptake(model, conditions["CN22"])
        glucose = model.reactions.get_by_id(C_SOURCE)
        assert glucose.lower_bound == pytest.approx(-1.05 * 1.477902724)
        assert glucose.upper_bound == 0.0


def test_maintenance_ships_pinned_and_is_released(ec_model):
    assert ec_model.reactions.get_by_id(NGAM_RXN).bounds == (0.7, 0.7)
    with ec_model as model:
        free_ngam(model)
        assert model.reactions.get_by_id(NGAM_RXN).bounds == (0.0, 1000.0)


def test_chemostat_fixes_growth_and_pins_uptake(ec_model, conditions):
    condition = conditions["CN4"]
    with ec_model as model:
        free_ngam(model)
        constrain_byproducts(model, condition)
        uptake = set_chemostat_constraints(model, condition, glucose=condition.glucose)

        biomass = model.reactions.get_by_id(BIO_RXN)
        assert biomass.lower_bound == pytest.approx(0.99 * condition.d_rate)

        glucose = model.reactions.get_by_id(C_SOURCE)
        assert glucose.upper_bound == pytest.approx(-uptake)
        assert glucose.lower_bound == pytest.approx(-1.01 * uptake)
        assert uptake >= 0.99 * condition.glucose * (1 - 1e-9)

        assert model.objective.direction == "min"
        assert POOL_RXN in str(model.objective.expression)


def test_the_measured_uptake_barely_covers_the_dilution_rate(ec_model, conditions):
    """Left to minimise uptake, the model needs slightly more glucose than
    was measured -- 1.158 against 1.116 mmol/gDW/h for CN4. The five
    percent of headroom the pipeline allows on the measurement is what
    makes the condition feasible at all, so it is not a formality."""
    condition = conditions["CN4"]
    with ec_model as model:
        free_ngam(model)
        constrain_byproducts(model, condition)
        uptake = set_chemostat_constraints(model, condition)
        assert condition.glucose < uptake <= 1.05 * condition.glucose
