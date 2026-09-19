from dataclasses import dataclass
from typing import Any

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from overflow.complexes import complex_subunits, fix_complexes
from overflow.config import OXPHOS_RXNS


@dataclass
class _Ec:
    rxns: list
    enzymes: list
    rxn_enz_mat: Any


class _Model:
    def __init__(self, rxns, enzymes, matrix):
        self.ec = _Ec(rxns, enzymes, csr_matrix(np.array(matrix, dtype=float)))


def test_subunits_are_read_with_their_stoichiometry():
    model = _Model(["r_1_EXP_1"], ["P1", "P2", "P3"], [[1, 2, 0]])
    (reaction_id, subunits), = complex_subunits(model, "r_1")
    assert reaction_id == "r_1_EXP_1"
    assert subunits == [("P1", 1.0), ("P2", 2.0)]


def test_a_base_reaction_matches_all_its_isozyme_variants():
    model = _Model(["r_1_EXP_1", "r_1_EXP_2", "r_10_EXP_1"],
                   ["P1", "P2"], [[1, 0], [0, 1], [1, 1]])
    assert [v[0] for v in complex_subunits(model, "r_1")] == ["r_1_EXP_1", "r_1_EXP_2"]


def test_subunits_are_levelled_to_one_abundance_per_complex():
    """Two subunits at 1:1 measured 2 and 4 both become 3."""
    model = _Model(["r_1"], ["P1", "P2"], [[1, 1]])
    result = fix_complexes(model, ["r_1"], ["P1", "P2"], np.array([2.0, 4.0]))
    assert dict(zip(result.uniprot_ids, result.abundances)) == {"P1": 3.0, "P2": 3.0}


def test_stoichiometry_scales_the_shared_abundance():
    """A subunit present twice per complex gets twice the abundance."""
    model = _Model(["r_1"], ["P1", "P2"], [[1, 2]])
    result = fix_complexes(model, ["r_1"], ["P1", "P2"], np.array([3.0, 6.0]))
    values = dict(zip(result.uniprot_ids, result.abundances))
    assert values["P1"] == pytest.approx(3.0)
    assert values["P2"] == pytest.approx(6.0)


def test_an_unmeasured_subunit_is_added_at_its_share():
    model = _Model(["r_1"], ["P1", "P2"], [[1, 1]])
    result = fix_complexes(model, ["r_1"], ["P1"], np.array([5.0]))
    values = dict(zip(result.uniprot_ids, result.abundances))
    assert result.added == ["P2"]
    assert values["P2"] == pytest.approx(5.0)


def test_a_complex_with_no_measurement_is_left_alone():
    model = _Model(["r_1"], ["P1", "P2"], [[1, 1]])
    result = fix_complexes(model, ["r_1"], ["P3"], np.array([1.0]))
    assert result.uniprot_ids == ["P3"]
    assert result.added == []


def test_zero_measurements_do_not_drag_the_complex_down():
    """A subunit measured at zero is treated as unseen, not as absent."""
    model = _Model(["r_1"], ["P1", "P2"], [[1, 1]])
    result = fix_complexes(model, ["r_1"], ["P1", "P2"], np.array([4.0, 0.0]))
    values = dict(zip(result.uniprot_ids, result.abundances))
    assert values["P1"] == pytest.approx(4.0)
    assert values["P2"] == pytest.approx(4.0)


@pytest.mark.slow
def test_respiratory_complexes_are_found_in_the_real_model(ec_model):
    total = sum(len(complex_subunits(ec_model, base)) for base in OXPHOS_RXNS)
    assert total > 0
    for base in OXPHOS_RXNS:
        for _, subunits in complex_subunits(ec_model, base):
            assert all(s > 0 for _, s in subunits)
