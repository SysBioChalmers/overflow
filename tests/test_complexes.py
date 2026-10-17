from dataclasses import dataclass
from typing import Any

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from overflow.complexes import complex_families, complex_subunits, level_complexes


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


def _families(rows):
    """``rows``: {reaction: {enzyme: copies}} -> (model, families)."""
    enzymes = sorted({e for r in rows.values() for e in r})
    matrix = [[r.get(e, 0) for e in enzymes] for r in rows.values()]
    model = _Model(list(rows), enzymes, matrix)
    return model, complex_families(model)


def _levelled(model, ids, values):
    result = level_complexes(model, ids, np.array(values, dtype=float))
    return dict(zip(result.uniprot_ids, result.abundances)), result


def test_subunits_are_levelled_to_the_median_per_copy():
    # per copy: P1 = 1, P2 = 4 / 2 = 2, P3 = 6; median 2
    model, _ = _families({"r_1": {"P1": 1, "P2": 2, "P3": 1}})
    values, result = _levelled(model, ["P1", "P2", "P3"], [1.0, 4.0, 6.0])
    assert values == pytest.approx({"P1": 2.0, "P2": 4.0, "P3": 2.0})
    assert sorted(result.adjusted) == ["P1", "P3"]


def test_an_unmeasured_subunit_is_added_at_its_share():
    model, _ = _families({"r_1": {"P1": 1, "P2": 1}})
    values, result = _levelled(model, ["P1"], [5.0])
    assert result.added == ["P2"]
    assert values["P2"] == pytest.approx(5.0)


def test_a_complex_with_no_measurement_is_left_alone():
    model, _ = _families({"r_1": {"P1": 1, "P2": 1}, "r_2": {"P3": 1}})
    values, result = _levelled(model, ["P3"], [1.0])
    assert values == {"P3": 1.0}
    assert result.added == []


def test_zero_measurements_count_as_unseen():
    model, _ = _families({"r_1": {"P1": 1, "P2": 1}})
    values, _ = _levelled(model, ["P1", "P2"], [4.0, 0.0])
    assert values == pytest.approx({"P1": 4.0, "P2": 4.0})


def test_isozyme_variants_are_one_complex():
    _, families = _families({"r_1_EXP_1": {"A": 1, "B": 1, "C": 1}, "r_1_EXP_2": {"A": 1, "B": 1, "D": 1}})
    assert sorted(families.members("A")) == ["A", "B", "C", "D"]


def test_one_shared_protein_does_not_chain_complexes_together():
    _, families = _families({"r_1": {"A": 1, "B": 1}, "r_2": {"A": 1, "C": 1}, "r_3": {"A": 1, "D": 1, "E": 1, "F": 1}})
    assert families.shared == {"A"}
    assert families.members("B") == ["B"]


def test_a_carrier_in_several_complexes_keeps_its_measurement():
    rows = {"r_1_EXP_1": {"B": 1, "C": 1, "D": 1, "S": 1}, "r_1_EXP_2": {"B": 1, "C": 1, "D": 1, "H": 1},
            "r_2": {"E": 1, "S": 1}, "r_3": {"F": 1, "S": 1}}
    model, families = _families(rows)
    assert "S" in families.shared
    values, _ = _levelled(model, ["B", "C", "D", "S"], [1.0, 1.0, 1.0, 50.0])
    assert values["S"] == 50.0
    assert values["B"] == 1.0


def test_a_subunit_in_two_complexes_belongs_to_the_one_with_more_reactions():
    rows = {"r_1": {"A": 1, "B": 1, "C": 1, "D": 1}, "r_2": {"A": 1, "B": 1, "C": 1, "D": 1},
            "r_3": {"A": 1, "B": 1, "E": 1, "F": 1, "G": 1, "H": 1}}
    _, families = _families(rows)
    assert sorted(families.members("A")) == ["A", "B", "C", "D"]
    assert sorted(families.members("E")) == ["E", "F", "G", "H"]


def test_an_enzyme_that_also_works_alone_is_shared():
    _, families = _families({"r_1": {"A": 1, "B": 1}, "r_2": {"A": 1}})
    assert "A" in families.shared


@pytest.mark.slow
def test_the_respiratory_complexes_are_found_in_the_real_model(ec_model):
    families = complex_families(ec_model)
    complex_iv = families.members("P00401")  # COX1
    complex_iii = families.members("P07256")  # COR1
    assert "P00410" in complex_iv and "P04037" in complex_iv  # COX2, COX4
    assert "P07257" in complex_iii and "P00401" not in complex_iii  # QCR2
    assert {"P00044", "P00045", "P09624"} <= families.shared  # cytochrome c (CYC1, CYC7), LPD1
