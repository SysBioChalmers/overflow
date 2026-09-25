"""Filter, unit conversion and the derived quantities.

The filter is checked against the abundances the MATLAB implementation
wrote into ``legacy_matlab/results/modelGeneration/``, which are the
only committed record of its output.
"""
import numpy as np
import pandas as pd
import pytest

from overflow.config import CONDITION_ORDER, ROOT, load_conditions
from overflow.proteomics import (
    condition_prot_data,
    drop_proteins,
    f_factor,
    filter_prot_data,
    mean_abundances,
    molecular_masses,
    read_proteomics,
    replicate_matrix,
    rescaled_p_tot,
    to_mass,
)

LEGACY = ROOT / "legacy_matlab" / "results" / "modelGeneration"


@pytest.fixture(scope="module")
def table():
    return read_proteomics()


@pytest.fixture(scope="module")
def masses():
    return molecular_masses()


# --- filter semantics -------------------------------------------------

def test_abundance_is_the_upper_edge_of_the_replicate_interval():
    data = np.array([[1.0, 2.0, 3.0]])
    result = filter_prot_data(["P1"], data)
    expected = 2.0 + 1.96 * np.std([1.0, 2.0, 3.0], ddof=1)
    assert result.abundances[0] == pytest.approx(expected)


def test_a_missing_replicate_drops_the_protein():
    data = np.array([[1.0, 2.0, np.nan]])
    assert filter_prot_data(["P1"], data).n_kept == 0


def test_a_protein_absent_from_most_replicates_is_dropped():
    data = np.array([[0.0, 0.0, 3.0]])
    assert filter_prot_data(["P1"], data).n_kept == 0


def test_a_protein_that_varies_more_than_its_median_is_dropped():
    data = np.array([[0.01, 0.01, 10.0]])
    assert filter_prot_data(["P1"], data).n_kept == 0


def test_median_filter_is_the_stricter_of_the_two():
    data = np.array([[1.0, 1.0, 3.2]])
    assert filter_prot_data(["P1"], data, median_filter=True).n_kept == 0
    assert filter_prot_data(["P1"], data, median_filter=False).n_kept == 1


def test_empty_identifiers_are_skipped():
    data = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    assert filter_prot_data(["", "P2"], data).uniprot_ids == ["P2"]


# --- parity with the MATLAB implementation ----------------------------

def _legacy_modified(condition):
    table = pd.read_csv(LEGACY / f"modifiedEnzymes_{condition}.txt", sep="\t")
    return {
        row["protein_IDs"].replace("prot_", ""): float(row["previous_values"])
        for _, row in table.iterrows()
    }


@pytest.mark.parametrize("condition", CONDITION_ORDER)
def test_filter_matches_the_abundances_matlab_carried_forward(table, condition):
    """Every abundance the MATLAB run reported as a starting value either
    equals this filter's output, or exceeds it because the minimum-usage
    pass raised it before flexibilization. It is never lower."""
    ids, data = replicate_matrix(table, condition)
    kept = filter_prot_data(ids, data)
    ours = dict(zip(kept.uniprot_ids, kept.abundances))

    exact = 0
    for protein, previous in _legacy_modified(condition).items():
        if protein not in ours:
            continue
        if previous == pytest.approx(ours[protein], rel=1e-9):
            exact += 1
        else:
            assert previous > ours[protein], (
                f"{protein}: MATLAB carried {previous:.6e} forward, below this "
                f"filter's {ours[protein]:.6e}"
            )
    assert exact >= 1


def test_cn4_filter_output_matches_matlab_exactly(table):
    """CN4 flexibilized nothing before integration, so every one of its
    modified enzymes pins the filter output exactly."""
    ids, data = replicate_matrix(table, "CN4")
    kept = filter_prot_data(ids, data)
    ours = dict(zip(kept.uniprot_ids, kept.abundances))
    legacy = _legacy_modified("CN4")
    assert len(legacy) == 7
    for protein, previous in legacy.items():
        assert ours[protein] == pytest.approx(previous, rel=1e-12)


# --- reading ----------------------------------------------------------

def test_duplicated_protein_is_collapsed(table):
    assert table["Protein.IDs"].duplicated().sum() == 0
    assert (table["Protein.IDs"] == "P61830").sum() == 1
    assert len(table) == 2928


def test_disagreeing_duplicates_are_rejected(tmp_path):
    path = tmp_path / "prot.txt"
    path.write_text(
        "Protein.IDs\tGene\tCN4_1_abs\n" "P1\tYA\t1e-6\n" "P1\tYB\t2e-6\n"
    )
    with pytest.raises(ValueError, match="differing"):
        read_proteomics(path)


def test_hgr_has_four_replicates(table):
    _, data = replicate_matrix(table, "hGR")
    assert data.shape == (2928, 4)


# --- unit conversion --------------------------------------------------

def test_daltons_convert_mmol_to_milligrams(masses):
    result = to_mass(["P00330"], np.array([1e-6]), masses)
    assert result[0] == pytest.approx(1e-6 * masses["P00330"])


def test_unknown_protein_is_an_error_not_a_zero():
    with pytest.raises(KeyError, match="molecular mass"):
        to_mass(["NOT_A_PROTEIN"], np.array([1.0]), {"P1": 1.0})


def test_every_measured_protein_has_a_mass(table, masses):
    missing = [p for p in table["Protein.IDs"] if p not in masses]
    assert missing == []


@pytest.mark.parametrize("condition", CONDITION_ORDER)
def test_measured_proteome_mass_is_close_to_the_measured_protein_content(
    table, masses, condition
):
    """A unit slip anywhere in the conversion shows up here as orders of
    magnitude, not percent."""
    ids, means = mean_abundances(table, condition)
    total_mg = float(np.nansum(to_mass(ids, means, masses)))
    p_tot_mg = load_conditions()[condition].p_tot * 1000
    assert 0.5 * p_tot_mg < total_mg < 1.6 * p_tot_mg


# --- derived quantities -----------------------------------------------

def test_f_factor_is_a_fraction_of_the_measured_mass():
    enzymes = ["P1", "P2"]
    ids = ["P1", "P2", "P3"]
    mass = np.array([1.0, 2.0, 7.0])
    assert f_factor(enzymes, ids, mass) == pytest.approx(0.3)


def test_rescaling_is_a_molar_ratio(table):
    condition = load_conditions()["CN4"]
    ids, means = mean_abundances(table, condition)
    _, data = replicate_matrix(table, condition)
    kept = filter_prot_data(ids, data)
    expected = condition.p_tot * (kept.abundances.sum() / np.nansum(means))
    assert rescaled_p_tot(condition, means, kept) == pytest.approx(expected)


@pytest.mark.slow
@pytest.mark.parametrize("condition", CONDITION_ORDER)
def test_condition_prot_data_is_positive_milligrams(ec_model, table, masses, condition):
    result = condition_prot_data(
        load_conditions()[condition], ec_model, table=table, masses=masses,
        fix_complex_subunits=False,
    )
    assert (result.prot_data.abundances > 0).all()
    assert len(result.prot_data.uniprot_ids) == result.filtered.n_kept
    assert 0.3 < result.f_factor < 0.6
    assert 0.1 < result.p_tot < 1.0
    in_model = set(result.prot_data.uniprot_ids) & set(ec_model.ec.enzymes)
    assert len(in_model) > 500


@pytest.mark.slow
def test_f_factor_agrees_with_geckopy(ec_model, table, masses):
    """This module computes the f-factor without needing a model; it must
    not drift from geckopy's."""
    from geckopy import calculate_f_factor
    from geckopy.databases import ProtData

    ids, means = mean_abundances(table, "CN4")
    mass = to_mass(ids, means, masses)
    ours = f_factor(ec_model.ec.enzymes, ids, mass)
    theirs = calculate_f_factor(ec_model, ProtData(uniprot_ids=ids, abundances=mass))
    assert ours == pytest.approx(theirs, rel=1e-12)


def test_dropping_proteins_keeps_the_others_and_their_values():
    ids, values = drop_proteins(["A", "B", "C"], np.array([1.0, 2.0, 3.0]), ["B", "Z"])
    assert ids == ["A", "C"]
    assert list(values) == [1.0, 3.0]


def test_nothing_is_dropped_by_default():
    ids, values = drop_proteins(["A", "B"], np.array([1.0, 2.0]), ())
    assert ids == ["A", "B"] and list(values) == [1.0, 2.0]


@pytest.mark.slow
def test_an_unmeasured_enzyme_is_left_out_of_the_abundances(ec_model, table, masses):
    condition = load_conditions()["CN38"]
    kwargs = dict(table=table, masses=masses, fix_complex_subunits=False)
    with_rki1 = condition_prot_data(condition, ec_model, **kwargs)
    without = condition_prot_data(condition, ec_model, unmeasured=["Q12189"], **kwargs)
    assert "Q12189" in with_rki1.prot_data.uniprot_ids
    assert "Q12189" not in without.prot_data.uniprot_ids
    assert len(without.prot_data.uniprot_ids) == len(with_rki1.prot_data.uniprot_ids) - 1
    assert without.f_factor == with_rki1.f_factor and without.p_tot == with_rki1.p_tot
