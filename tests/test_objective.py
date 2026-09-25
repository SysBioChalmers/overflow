import pytest

from overflow.build import OBJECTIVES, build_condition, metabolic_reactions
from overflow.config import POOL_RXN
from overflow.constraints import block_reactions
from routes import two_route_model
from tiny import tiny_ec_model


def test_metabolic_reactions_leave_out_enzyme_usage_and_the_pool():
    model = tiny_ec_model()
    ids = {r.id for r in metabolic_reactions(model)}
    assert POOL_RXN in {r.id for r in model.reactions}
    assert POOL_RXN not in ids
    assert not any(i.startswith("usage_prot_") for i in ids)
    assert any(r.id.startswith("usage_prot_") for r in model.reactions)
    assert len(ids) > 0


def test_an_unknown_objective_is_refused(conditions):
    with pytest.raises(ValueError, match="objective must be one of"):
        build_condition(conditions["CN4"], objective="speed")


@pytest.mark.parametrize("objective", [o for o in OBJECTIVES if o != "protein"])
def test_flux_objectives_need_rate_fitting(conditions, objective):
    with pytest.raises(ValueError, match="needs fit_rates"):
        build_condition(conditions["CN4"], objective=objective)


def test_parsimony_over_a_subset_ignores_the_others_but_returns_every_flux():
    from overflow.ngam import parsimonious_solution

    model = two_route_model()
    metabolic = [r for r in model.reactions if not r.id.startswith("usage_prot_")]

    everything = parsimonious_solution(model)
    assert everything.fluxes["step1"] == pytest.approx(1.0)
    assert everything.fluxes["direct"] == pytest.approx(0.0)

    subset = parsimonious_solution(model, metabolic)
    assert subset.fluxes["direct"] == pytest.approx(1.0)
    assert subset.fluxes["step1"] == pytest.approx(0.0)
    assert subset.fluxes["usage_prot_E"] == pytest.approx(10.0)
    assert set(subset.fluxes.index) == {r.id for r in model.reactions}


def test_a_blocked_reaction_stays_shut_on_a_small_model():
    model = two_route_model()
    block_reactions(model, ["direct"])
    assert model.reactions.direct.bounds == (0.0, 0.0)
    assert model.slim_optimize() == pytest.approx(1.0)
    assert model.optimize().fluxes["step1"] == pytest.approx(1.0)


def test_blocking_an_unknown_reaction_is_refused_on_a_small_model():
    with pytest.raises(KeyError, match="cannot block r_missing"):
        block_reactions(two_route_model(), ["r_missing"])
