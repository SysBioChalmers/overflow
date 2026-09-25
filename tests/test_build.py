"""End-to-end build of one condition.

Marked ``integration``: it solves several hundred LPs and takes minutes.
Run it with ``pytest -m integration``.
"""
import numpy as np
import pytest

from overflow.build import build_condition, summary_table
from overflow.config import BIO_RXN, C_SOURCE, POOL_RXN

pytestmark = [pytest.mark.integration, pytest.mark.slow]


@pytest.fixture(scope="module")
def built(conditions):
    return build_condition(
        conditions["CN4"], ngam_steps=20, uptake_flex=1.08, verbose=False
    )


def test_the_model_grows_at_the_dilution_rate(built, conditions):
    target = conditions["CN4"].d_rate
    assert built.growth >= 0.99 * target


def test_uptake_stays_within_the_measurement_plus_headroom(built, conditions):
    measured = conditions["CN4"].glucose
    assert 0.99 * measured <= built.glucose <= 1.08 * measured


def test_undetected_byproducts_carry_no_flux(built):
    for reaction_id in ("r_1808", "r_1634", "r_1761", "r_1793"):
        assert abs(float(built.fluxes[reaction_id])) < 1e-9


def test_protein_usage_stays_within_the_budget(built):
    assert float(built.fluxes[POOL_RXN]) <= built.flexibilization.pool_after + 1e-6


def test_every_measured_enzyme_keeps_a_concentration(built):
    concentrations = np.asarray(built.model.ec.concs, dtype=float)
    assert int((~np.isnan(concentrations)).sum()) == built.n_measured
    assert np.nanmin(concentrations) > 0


def test_flexibilization_only_ever_raises(built):
    table = built.flexibilization.table()
    if len(table):
        assert (table["modified_mg_gDW"] >= table["previous_mg_gDW"]).all()


def test_the_summary_reports_the_condition(built):
    table = summary_table([built])
    assert list(table["condition"]) == ["CN4"]
    assert table.loc[0, "growth"] >= 0.099
