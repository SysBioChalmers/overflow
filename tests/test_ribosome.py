"""Putting the ribosome in the way of protein synthesis.

Every number here is worked out by hand from the tiny model: one unit of
protein costs two amino acids, so a ribosome elongating at 10.5 aa/s
turns over 5.25 times a second in units of protein, and a 30 kDa subunit
costs 30000 x 2 / (10.5 x 3600) = 1.5873 mg per mmol of protein made.
"""
import numpy as np
import pandas as pd
import pytest

from overflow.ribosome import (
    AA_PER_SECOND,
    AMINO_ACID_MET,
    MIN_MEAN_ABUNDANCE,
    RibosomeSubunits,
    add_ribosome,
    amino_acid_demand,
    candidate_means,
    core_subunits,
    read_ribosome,
    subunit_coefficient,
    translation_kcat,
)
from tiny import TINY_AA_PER_PROTEIN, tiny_translation_model

SUBUNIT_MASS = 30_000.0
EXPECTED_COEFFICIENT = SUBUNIT_MASS * TINY_AA_PER_PROTEIN / (AA_PER_SECOND * 3600.0)


@pytest.fixture
def model():
    return tiny_translation_model(pool=200.0, supply=10.0)


@pytest.fixture
def subunits():
    ids = ["R1", "R2"]
    return RibosomeSubunits(
        uniprot_ids=ids,
        genes={u: f"gene_{u}" for u in ids},
        masses={u: SUBUNIT_MASS for u in ids},
        mean_abundance={u: 1e-4 for u in ids},
        n_candidates=2,
    )


def _translation(model):
    return next(
        r for r in model.reactions
        if r.id == "translation" or r.id.startswith("translation_")
    )


# --- the arithmetic ---------------------------------------------------

def test_amino_acid_demand_is_read_from_the_pseudoreaction(model):
    assert amino_acid_demand(model, "r_protein") == pytest.approx(TINY_AA_PER_PROTEIN)


def test_the_ribosome_turns_over_once_per_protein_it_finishes(model):
    """10.5 amino acids a second, two per protein, so 5.25 proteins a second."""
    assert translation_kcat(model, protein_rxn="r_protein") == pytest.approx(5.25)


def test_a_protein_that_takes_more_amino_acids_costs_more_ribosome():
    cheap = subunit_coefficient(SUBUNIT_MASS, demand=2.0)
    dear = subunit_coefficient(SUBUNIT_MASS, demand=4.0)
    assert dear == pytest.approx(2 * cheap)


# --- the rewiring -----------------------------------------------------

def test_protein_now_comes_out_of_translation(model, subunits):
    add_ribosome(model, subunits, protein_rxn="r_protein")

    pseudoreaction = model.reactions.get_by_id("r_protein")
    protein = model.metabolites.get_by_id("protein")
    assert protein not in pseudoreaction.metabolites
    assert model.metabolites.get_by_id(AMINO_ACID_MET) in pseudoreaction.metabolites

    translation = _translation(model)
    assert translation.metabolites[protein] == pytest.approx(1.0)
    assert translation.metabolites[model.metabolites.get_by_id(AMINO_ACID_MET)] == pytest.approx(-1.0)


def test_each_subunit_becomes_an_enzyme_drawing_on_the_pool(model, subunits):
    add_ribosome(model, subunits, protein_rxn="r_protein")

    for unit in subunits.uniprot_ids:
        assert unit in model.ec.enzymes
        index = model.ec.enzymes.index(unit)
        assert model.ec.mw[index] == pytest.approx(SUBUNIT_MASS)
        usage = model.reactions.get_by_id(f"usage_prot_{unit}")
        assert any(m.id == "prot_pool" for m in usage.metabolites)


def test_the_subunit_cost_matches_the_elongation_rate(model, subunits):
    """The round trip that matters: a kcat written through geckopy has to
    come back out as the protein cost the elongation rate implies. This
    is where a confusion between daltons and kDa would show up, as a
    thousandfold error in the cost of translation."""
    add_ribosome(model, subunits, protein_rxn="r_protein")

    translation = _translation(model)
    for unit in subunits.uniprot_ids:
        coefficient = translation.metabolites[model.metabolites.get_by_id(f"prot_{unit}")]
        assert coefficient == pytest.approx(-EXPECTED_COEFFICIENT)
    assert EXPECTED_COEFFICIENT == pytest.approx(1.5873, abs=1e-4)


# --- what it does to the model ---------------------------------------

def test_growth_is_unchanged_when_ribosomes_are_affordable(model, subunits):
    before = model.optimize().objective_value
    add_ribosome(model, subunits, protein_rxn="r_protein")
    after = model.optimize().objective_value
    assert after == pytest.approx(before)


def test_the_ribosome_is_paid_for_out_of_the_protein_pool(model, subunits):
    before = model.optimize().fluxes["prot_pool_exchange"]
    add_ribosome(model, subunits, protein_rxn="r_protein")
    after = model.optimize().fluxes["prot_pool_exchange"]

    growth = model.optimize().objective_value
    expected = before + len(subunits.uniprot_ids) * growth * EXPECTED_COEFFICIENT
    assert after == pytest.approx(expected)


def test_without_ribosomes_there_is_no_protein_and_no_growth(model, subunits):
    """The point of the rewiring: protein cannot bypass translation."""
    add_ribosome(model, subunits, protein_rxn="r_protein")
    for unit in subunits.uniprot_ids:
        model.reactions.get_by_id(f"usage_prot_{unit}").upper_bound = 0.0

    solution = model.optimize()
    assert (solution.objective_value or 0.0) == pytest.approx(0.0, abs=1e-9)


def test_a_capped_subunit_limits_growth_in_proportion(model, subunits):
    """Half the ribosome, half the protein: growth 5 becomes 2.5."""
    add_ribosome(model, subunits, protein_rxn="r_protein")
    full = model.optimize().objective_value
    for unit in subunits.uniprot_ids:
        model.reactions.get_by_id(f"usage_prot_{unit}").upper_bound = (
            0.5 * full * EXPECTED_COEFFICIENT
        )

    assert model.optimize().objective_value == pytest.approx(0.5 * full)


def test_the_result_reports_what_it_added(model, subunits):
    result = add_ribosome(model, subunits, protein_rxn="r_protein")
    assert set(result.enzymes_added) == {"R1", "R2"}
    assert result.amino_acid_demand == pytest.approx(TINY_AA_PER_PROTEIN)
    assert result.kcat == pytest.approx(5.25)
    assert result.reactions_added


# --- refusals ---------------------------------------------------------

def test_adding_the_ribosome_twice_is_refused(model, subunits):
    add_ribosome(model, subunits, protein_rxn="r_protein")
    with pytest.raises(ValueError, match="already"):
        add_ribosome(model, subunits, protein_rxn="r_protein")


def test_an_empty_subunit_set_is_refused(model):
    empty = RibosomeSubunits(uniprot_ids=[], genes={}, masses={}, mean_abundance={})
    with pytest.raises(ValueError, match="no ribosomal subunits"):
        add_ribosome(model, empty, protein_rxn="r_protein")


def test_a_pseudoreaction_without_protein_is_refused(model, subunits):
    model.reactions.get_by_id("r_protein").subtract_metabolites(
        {model.metabolites.get_by_id("protein"): 1.0}
    )
    with pytest.raises(ValueError, match="protein"):
        add_ribosome(model, subunits, protein_rxn="r_protein")


# --- selecting the core ----------------------------------------------

def _ribosome_table(ids):
    return pd.DataFrame(
        {
            "uniprot": ids,
            "gene_name": [f"N{i}" for i in ids],
            "gene": [f"Y{i}" for i in ids],
            "mass": [30000.0] * len(ids),
            "sequence": ["M"] * len(ids),
        }
    )


def _proteomics(values):
    return pd.DataFrame(
        {
            "Protein.IDs": list(values),
            "Gene": [f"Y{i}" for i in values],
            "CN4_1_abs": [v[0] for v in values.values()],
            "CN4_2_abs": [v[1] for v in values.values()],
        }
    )


def test_a_subunit_below_the_threshold_is_left_out():
    table = _ribosome_table(["A", "B"])
    data = _proteomics({"A": (1e-4, 1e-4), "B": (1e-6, 1e-6)})
    selected = core_subunits(table, data, threshold=MIN_MEAN_ABUNDANCE)
    assert selected.uniprot_ids == ["A"]


def test_a_subunit_that_was_never_measured_is_left_out():
    table = _ribosome_table(["A", "B"])
    data = _proteomics({"A": (1e-4, 1e-4)})
    assert core_subunits(table, data).uniprot_ids == ["A"]


def test_a_subunit_measured_only_as_missing_is_left_out():
    table = _ribosome_table(["A"])
    data = _proteomics({"A": (np.nan, np.nan)})
    assert core_subunits(table, data).uniprot_ids == []


def test_selection_averages_over_every_replicate_of_every_condition():
    """A subunit seen strongly in one condition and not at all in another
    still belongs to the ribosome; the average decides, not one column."""
    table = _ribosome_table(["A"])
    data = _proteomics({"A": (2e-5, 2e-6)})
    assert core_subunits(table, data).uniprot_ids == ["A"]
    assert core_subunits(table, data).mean_abundance["A"] == pytest.approx(1.1e-5)


def test_candidate_means_keep_the_subunits_below_the_threshold():
    table = _ribosome_table(["A", "B", "C"])
    data = _proteomics({"A": (1e-4, 3e-4), "B": (1e-7, 1e-7)})
    means = candidate_means(table, data)
    assert list(means.index) == ["A", "B"]
    assert means["A"] == pytest.approx(2e-4)
    assert means["B"] == pytest.approx(1e-7)


def test_the_density_integrates_to_one_and_peaks_where_the_data_do():
    from overflow.plots import subunit_density

    log_values = np.r_[np.full(30, -4.0), np.full(10, -6.0)]
    grid, density = subunit_density(log_values, bandwidth=0.1)
    area = float(np.sum(0.5 * (density[1:] + density[:-1]) * np.diff(grid)))
    assert area == pytest.approx(1.0, abs=0.02)
    assert grid[np.argmax(density)] == pytest.approx(-4.0, abs=0.1)
    assert grid.min() == pytest.approx(-6.3) and grid.max() == pytest.approx(-3.7)


def test_the_figure_is_the_published_matlab_plot(tmp_path):
    from overflow.plots import subunit_abundance_figure

    means = pd.Series({f"P{i}": 10 ** x for i, x in enumerate(np.linspace(-7, -4, 40))})
    path = tmp_path / "riboSubunits.pdf"
    figure = subunit_abundance_figure(means, path)
    assert path.stat().st_size > 1000

    axis = figure.axes[0]
    assert axis.get_title() == "Distribution of average ribosomal subunit abundances"
    assert axis.title.get_fontweight() == "bold"
    assert axis.get_xlabel() == "Subunit abundance (log10(mmol/gDCW))"
    assert axis.get_ylabel() == "Density"
    line = axis.lines[0]
    assert line.get_color() == "#0072BD"
    assert axis.get_ylim()[0] == 0
    box = axis.get_position()
    assert box.width * figure.get_figwidth() * 72 == pytest.approx(341, abs=0.5)
    assert box.height * figure.get_figheight() * 72 == pytest.approx(247, abs=0.5)


# --- against the real data -------------------------------------------

@pytest.mark.slow
def test_the_real_ribosome_table_reads(): 
    table = read_ribosome()
    assert {"uniprot", "gene", "mass", "sequence"} <= set(table.columns)
    assert len(table) == 253
    assert table["mass"].min() > 1000


@pytest.mark.slow
def test_the_core_ribosome_is_selected_from_the_real_data():
    from overflow.proteomics import read_proteomics

    selected = core_subunits(read_ribosome(), read_proteomics())
    assert len(selected) == 48
    assert all(selected.mean_abundance[u] >= MIN_MEAN_ABUNDANCE for u in selected.uniprot_ids)


@pytest.mark.slow
@pytest.mark.parametrize("condition", ["CN4", "CN22", "CN38", "CN75", "hGR"])
def test_the_core_is_exactly_what_the_matlab_run_added(condition):
    """The GECKO 2 ecModel carries 816 enzymes and each finished
    condition model carries 864. Those 48 are the subunits ribosome.m
    selected, and this selection reproduces them one for one."""
    import scipy.io as sio

    from overflow.config import ROOT
    from overflow.proteomics import read_proteomics

    models = ROOT / "legacy_matlab" / "models"
    base = sio.loadmat(
        models / "ecYeastGEM.mat", struct_as_record=False, squeeze_me=True
    )["ecModel"].enzymes
    built = sio.loadmat(
        models / f"ecModel_P_{condition}.mat", struct_as_record=False, squeeze_me=True
    )[f"ecModelP_{condition}"].enzymes
    added = set(map(str, built)) - set(map(str, base))

    selected = core_subunits(read_ribosome(), read_proteomics())
    assert set(selected.uniprot_ids) == added


@pytest.mark.slow
def test_the_real_protein_pseudoreaction_sets_the_elongation_cost(ec_model):
    demand = amino_acid_demand(ec_model)
    assert demand == pytest.approx(4.191011, rel=1e-5)
    assert translation_kcat(ec_model) == pytest.approx(10.5 / demand)


# --- applying the measurements ---------------------------------------

def test_a_measured_subunit_is_capped_at_its_measurement(model, subunits):
    from overflow.ribosome import constrain_subunits

    add_ribosome(model, subunits, protein_rxn="r_protein")
    growth = model.optimize().objective_value
    generous = 2 * growth * EXPECTED_COEFFICIENT

    result = constrain_subunits(model, {"R1": generous, "R2": generous})
    assert result.adjusted == []
    assert model.reactions.get_by_id("usage_prot_R1").upper_bound == generous
    assert model.ec.concs[model.ec.enzymes.index("R1")] == generous


def test_a_subunit_measured_below_what_translation_needs_is_raised(model, subunits):
    """The cell was growing, so a measurement that would stop it is
    taken as an underestimate rather than as a limit."""
    from overflow.ribosome import constrain_subunits

    add_ribosome(model, subunits, protein_rxn="r_protein")
    growth = model.optimize().objective_value
    needed = growth * EXPECTED_COEFFICIENT

    result = constrain_subunits(model, {"R1": 0.1 * needed, "R2": 10 * needed})
    assert result.adjusted == ["R1"]
    assert result.required[0] == pytest.approx(1.01 * needed)
    assert model.optimize().objective_value == pytest.approx(growth, rel=1e-6)


def test_an_unmeasured_subunit_keeps_drawing_on_the_pool(model, subunits):
    from overflow.ribosome import constrain_subunits

    add_ribosome(model, subunits, protein_rxn="r_protein")
    constrain_subunits(model, {"R1": 1.0})
    assert model.reactions.get_by_id("usage_prot_R2").upper_bound == 1000.0


# --- reading the condition's subunit measurements --------------------

def test_subunit_abundances_are_filtered_like_any_other_protein():
    """A subunit belongs to the core on its average across conditions,
    but its cap in one condition still has to survive that condition's
    own quality filter."""
    from overflow.build_ribosome import subunit_abundances
    from overflow.config import load_conditions

    table = pd.DataFrame(
        {
            "Protein.IDs": ["A", "B"],
            "Gene": ["YA", "YB"],
            "CN4_1_abs": [1e-5, 1e-7],
            "CN4_2_abs": [1.1e-5, 9e-5],
            "CN4_3_abs": [1.05e-5, 1e-7],
        }
    )
    selected = RibosomeSubunits(
        uniprot_ids=["A", "B"],
        genes={"A": "YA", "B": "YB"},
        masses={"A": 30000.0, "B": 30000.0},
        mean_abundance={"A": 1e-5, "B": 3e-5},
    )
    result = subunit_abundances(
        load_conditions()["CN4"], selected, table=table,
        masses={"A": 30000.0, "B": 30000.0},
    )
    assert set(result) == {"A"}, "B varies more than its median and should drop out"
    assert result["A"] > 0


def test_subunit_abundances_are_milligrams():
    from overflow.build_ribosome import subunit_abundances
    from overflow.config import load_conditions

    table = pd.DataFrame(
        {
            "Protein.IDs": ["A"],
            "Gene": ["YA"],
            "CN4_1_abs": [1e-5],
            "CN4_2_abs": [1e-5],
            "CN4_3_abs": [1e-5],
        }
    )
    selected = RibosomeSubunits(
        uniprot_ids=["A"], genes={"A": "YA"}, masses={"A": 30000.0},
        mean_abundance={"A": 1e-5},
    )
    result = subunit_abundances(
        load_conditions()["CN4"], selected, table=table, masses={"A": 30000.0}
    )
    assert result["A"] == pytest.approx(1e-5 * 30000.0)


def test_ribosome_cost_is_proportional_to_how_fast_protein_is_made(model, subunits):
    """Translation is charged per unit of protein, so halving the supply
    halves the ribosome. Across the real conditions this shows up as
    11.5 mg/gDW at D=0.1 and 33.3 at D=0.29, a ratio of 2.90 against a
    dilution-rate ratio of 2.9."""
    add_ribosome(model, subunits, protein_rxn="r_protein")
    full = model.optimize()
    full_usage = abs(full.fluxes["usage_prot_R1"])

    model.reactions.get_by_id("EX_S").upper_bound /= 2
    half = model.optimize()
    assert half.objective_value == pytest.approx(full.objective_value / 2)
    assert abs(half.fluxes["usage_prot_R1"]) == pytest.approx(full_usage / 2)


def test_the_candidate_means_agree_with_the_selected_core():
    from overflow.proteomics import read_proteomics

    table = read_ribosome()
    means = candidate_means(table, read_proteomics())
    core = core_subunits(table, read_proteomics())
    assert int((means >= MIN_MEAN_ABUNDANCE).sum()) == len(core.uniprot_ids)
    assert set(means[means >= MIN_MEAN_ABUNDANCE].index) == set(core.uniprot_ids)


def test_a_subunit_missing_from_a_replicate_has_no_average():
    """MATLAB's mean propagates NaN, so the published selection and figure leave
    such a subunit out."""
    table = _ribosome_table(["A", "B"])
    data = _proteomics({"A": (1e-4, np.nan), "B": (1e-4, 1e-4)})
    assert list(candidate_means(table, data).index) == ["B"]
    assert core_subunits(table, data).uniprot_ids == ["B"]
