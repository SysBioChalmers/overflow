"""Abundance reconciliation treats an enzyme complex as one unit.

    supply -> S -+- R_E1 (complex E1a + E1b) -+-> P -> biomass
                 \\- R_E2 (enzyme E2, closed) -/

A unit of flux on R_E1 costs 10 mg of each subunit, so the numbers below
can be worked out by hand.
"""
import cobra
import numpy as np
import pytest
from geckopy import EcModel
from raven_toolbox.io.ec_data import EcData
from scipy import sparse

from overflow.complexes import complex_families
from overflow.flexibilize import (
    apply_concentrations,
    minimum_usage_pass,
    relax_to_growth,
    relax_to_measured_rates,
)
from tiny import protein_coefficient, tiny_adapter

ENZYMES = ["E1a", "E1b", "E2"]
COST = protein_coefficient(kcat=1.0, mw=36_000.0)  # 10 mg per unit flux


def complex_model(pool: float = 1000.0) -> EcModel:
    model = cobra.Model("tiny_complex")
    model.compartments = {"c": "cytoplasm"}
    mets = {n: cobra.Metabolite(n, compartment="c") for n in ("S", "P", "biomass", "prot_pool", *(f"prot_{e}" for e in ENZYMES))}
    model.add_metabolites(list(mets.values()))

    def reaction(rid, stoichiometry, lb=0.0, ub=1000.0):
        r = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        r.add_metabolites({mets[m]: c for m, c in stoichiometry.items()})
        return r

    model.add_reactions([
        reaction("EX_S", {"S": 1}, ub=10.0),
        reaction("R_E1", {"S": -1, "P": 1, "prot_E1a": -COST, "prot_E1b": -COST}),
        reaction("R_E2", {"S": -1, "P": 1, "prot_E2": -COST}, ub=0.0),
        reaction("BIO", {"P": -1, "biomass": 1}),
        reaction("EX_BIO", {"biomass": -1}),
        reaction("prot_pool_exchange", {"prot_pool": 1}, ub=pool),
        *(reaction(f"usage_prot_{e}", {"prot_pool": -1, f"prot_{e}": 1}) for e in ENZYMES),
    ])
    model.objective = "BIO"
    ec_model = EcModel.from_cobra(model, adapter=tiny_adapter())
    ec_model.ec = EcData(
        gecko_light=False, rxns=["R_E1", "R_E2"], kcat=np.array([1.0, 1.0]),
        source=["test"] * 2, notes=[""] * 2, eccodes=[""] * 2,
        genes=[f"gene_{e}" for e in ENZYMES], enzymes=list(ENZYMES), mw=np.full(3, 36_000.0),
        sequence=["M"] * 3, concs=np.full(3, np.nan),
        rxn_enz_mat=sparse.csr_matrix([[1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    ec_model.ec.validate()
    return ec_model


def concentrations(**values):
    return np.array([values.get(e, np.nan) for e in ENZYMES], dtype=float)


def test_the_two_subunits_are_one_complex():
    families = complex_families(complex_model())
    assert sorted(families.members("E1a")) == ["E1a", "E1b"]


def test_the_minimum_usage_pass_raises_the_whole_complex_by_one_factor():
    model = complex_model()
    model.ec.concs = concentrations(E1a=10.0, E1b=40.0)
    result = minimum_usage_pass(model, target_growth=5.0, bio_rxn="BIO")
    # 5 units need 50 mg of each; E1a needs 1.01 * 50 / 10 = 5.05-fold
    raised = dict(zip(ENZYMES, result.concentrations))
    assert raised["E1a"] == pytest.approx(10.0 * 5.05)
    assert raised["E1b"] == pytest.approx(40.0 * 5.05)
    assert sorted(result.raised) == ["E1a", "E1b"]


def test_the_rate_fit_relaxes_the_whole_complex_by_one_factor():
    model = complex_model()
    apply_concentrations(model, concentrations(E1a=10.0, E1b=40.0))
    result = relax_to_measured_rates(model, target_growth=5.0, bio_rxn="BIO")
    assert result.feasible
    folds = {p: r / m for p, m, r in zip(result.proteins, result.measured, result.relaxed)}
    assert folds["E1a"] == pytest.approx(5.0)
    assert folds["E1b"] == pytest.approx(5.0)


def test_a_subunit_measured_at_zero_does_not_inflate_its_complex_in_the_rate_fit():
    model = complex_model()
    apply_concentrations(model, concentrations(E1a=0.0, E1b=40.0))
    result = relax_to_measured_rates(model, target_growth=5.0, bio_rxn="BIO")
    assert result.feasible
    relaxed = dict(zip(result.proteins, result.relaxed))
    assert relaxed.get("E1b", 40.0) == pytest.approx(50.0, rel=1e-6)


def test_absolute_weighting_does_not_make_a_zero_measurement_free():
    model = complex_model()
    apply_concentrations(model, concentrations(E1a=0.0, E1b=40.0))
    result = relax_to_measured_rates(model, target_growth=5.0, bio_rxn="BIO", weight="absolute")
    assert result.feasible
    assert result.added_mg == pytest.approx(50.0 + 10.0, rel=1e-6)


def test_the_global_relaxation_raises_the_complex_together_and_leaves_the_pool():
    model = complex_model()
    apply_concentrations(model, concentrations(E1a=10.0, E1b=10.0))
    result = relax_to_growth(model, target_growth=5.0, bio_rxn="BIO", pool_rxn="prot_pool_exchange")
    assert result.reached_target
    assert result.pool_increases == 0
    caps = [model.reactions.get_by_id(f"usage_prot_{e}").upper_bound for e in ("E1a", "E1b")]
    assert caps == pytest.approx([50.0, 50.0])


def test_the_global_relaxation_grows_the_pool_only_when_it_has_to():
    model = complex_model(pool=60.0)  # 5 units need 100 mg in all
    apply_concentrations(model, concentrations(E1a=50.0, E1b=50.0))
    result = relax_to_growth(model, target_growth=5.0, bio_rxn="BIO", pool_rxn="prot_pool_exchange")
    assert result.reached_target
    assert result.pool_increases == 1
    assert model.reactions.prot_pool_exchange.upper_bound == pytest.approx(100.0)
    assert result.released == []
