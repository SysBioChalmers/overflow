import numpy as np
import pytest

from overflow.config import NGAM_RXN
from overflow.ngam import fit_ngam, relative_residual


def test_a_perfect_prediction_has_no_residual():
    assert relative_residual(np.array([1.0, 2.0]), np.array([1.0, 2.0])) == 0.0


def test_the_residual_is_relative_not_absolute():
    """Being wrong by one unit on a rate of two must count for more than
    being wrong by one unit on a rate of a hundred."""
    small = relative_residual(np.array([3.0]), np.array([2.0]))
    large = relative_residual(np.array([101.0]), np.array([100.0]))
    assert small > large


def test_residuals_over_several_rates_add_in_quadrature():
    residual = relative_residual(np.array([1.1, 2.2]), np.array([1.0, 2.0]))
    assert residual == pytest.approx(np.sqrt(0.1**2 + 0.1**2))


@pytest.mark.slow
def test_fitting_leaves_maintenance_as_a_requirement(ec_model, conditions):
    """The fitted value is a lower bound: maintenance is a demand the
    cell must meet, not a fixed expenditure."""
    condition = conditions["CN4"]
    with ec_model as model:
        from overflow.constraints import (
            constrain_byproducts,
            free_ngam,
            set_chemostat_constraints,
        )

        free_ngam(model)
        constrain_byproducts(model, condition)
        set_chemostat_constraints(model, condition)
        fit = fit_ngam(model, condition, bounds=(0.0, 5.0), steps=6)

        reaction = model.reactions.get_by_id(NGAM_RXN)
        assert reaction.lower_bound == pytest.approx(fit.value)
        assert reaction.upper_bound == 1000.0
        assert 0.0 <= fit.value <= 5.0
        assert np.isfinite(fit.error)
