"""Respiratory-chain kcats derived from the sampled fluxes and the measured complexes."""
import numpy as np
import pandas as pd
import pytest

from overflow import load_model
from overflow.config import RESULTS_DIR
from overflow.respiratory_kcats import (
    COMPLEXES,
    apparent_turnovers,
    apply_respiratory_kcats,
    derived_kcats,
    variants,
)

pytestmark = pytest.mark.slow


def test_the_derived_kcat_is_the_highest_turnover_and_names_its_condition():
    turnovers = pd.DataFrame({
        "condition": ["a", "b", "a", "b"],
        "complex": ["complex III", "complex III", "complex IV", "complex IV"],
        "reaction": ["r_0439", "r_0439", "r_0438", "r_0438"],
        "apparent_turnover_per_s": [100.0, 300.0, 50.0, 20.0],
        "ecmodel_kcat_per_s": [290.0, 290.0, 1376.0, 1376.0],
    })
    kcats = derived_kcats(turnovers).set_index("complex")
    assert kcats.loc["complex III", "kcat_per_s"] == 300.0
    assert kcats.loc["complex III", "condition"] == "b"
    assert kcats.loc["complex IV", "condition"] == "a"


def test_every_complex_has_a_finite_turnover_in_every_condition(ec_model, conditions):
    sampled = pd.read_csv(RESULTS_DIR / "randomSampling" / "allFluxes.tsv", sep="\t").set_index("rxnID")
    turnovers = apparent_turnovers(ec_model, conditions, sampled)
    assert len(turnovers) == len(conditions) * len(COMPLEXES)
    assert np.isfinite(turnovers["apparent_turnover_per_s"]).all()
    assert (turnovers["apparent_turnover_per_s"] > 0).all()
    assert np.allclose(turnovers.groupby("complex")["relative_to_highest"].max(), 1.0)


def test_applying_changes_only_the_respiratory_complexes():
    model, reference = load_model(), load_model()
    kcats = pd.DataFrame({"reaction": ["r_0438"], "kcat_per_s": [2 * float(model.ec.kcat[list(model.ec.rxns).index("r_0438_EXP_1")])]})
    updated = apply_respiratory_kcats(model, kcats)
    assert set(updated) == set(variants(model, "r_0438"))
    for rid in updated:
        new = {m.id: c for m, c in model.reactions.get_by_id(rid).metabolites.items() if m.id.startswith("prot_")}
        old = {m.id: c for m, c in reference.reactions.get_by_id(rid).metabolites.items() if m.id.startswith("prot_")}
        assert new.keys() == old.keys()
        assert all(np.isclose(new[k], old[k] / 2) for k in new)
    untouched = reference.reactions.get_by_id("r_0439_EXP_1").metabolites
    assert {m.id: c for m, c in model.reactions.get_by_id("r_0439_EXP_1").metabolites.items()} == {m.id: c for m, c in untouched.items()}
