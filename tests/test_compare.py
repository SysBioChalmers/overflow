"""The comparison against the published results."""
import pandas as pd
import pytest

from overflow.compare import LEGACY, exchange_comparison, report, usage_comparison
from overflow.config import CONDITION_ORDER, ROOT


@pytest.fixture(scope="module")
def results_dir(tmp_path_factory):
    """A results directory holding the legacy outputs, reshaped as this
    pipeline writes them. Lets the comparison be exercised without a run."""
    root = tmp_path_factory.mktemp("results")
    (root / "modelSimulation").mkdir()
    (root / "enzymeUsage").mkdir()

    # GECKO 2 splits reversible reactions, so an uptake lives on the
    # _REV copy as a positive flux. Fold those back onto the base
    # reaction as a negative flux, which is this pipeline's format.
    reversed_uptakes = {"r_1714_REV": "r_1714", "r_1992_REV": "r_1992"}

    for condition in CONDITION_ORDER:
        rows = {}
        for line in (LEGACY / "modelSimulation" / f"allFluxes_{condition}.txt").read_text().splitlines():
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            try:
                flux = float(parts[3])
            except ValueError:
                continue
            reaction = parts[0]
            if reaction in reversed_uptakes:
                rows[reversed_uptakes[reaction]] = {
                    "rxnID": reversed_uptakes[reaction], "rxnName": parts[1],
                    "equation": parts[2], "flux": -flux,
                }
            elif reaction not in rows:
                rows[reaction] = {
                    "rxnID": reaction, "rxnName": parts[1],
                    "equation": parts[2], "flux": flux,
                }
        pd.DataFrame(list(rows.values())).to_csv(
            root / "modelSimulation" / f"allFluxes_{condition}.tsv", sep="\t", index=False
        )

    pd.read_csv(LEGACY / "enzymeUsage" / "capUsage.txt", sep="\t").to_csv(
        root / "enzymeUsage" / "capUsage.tsv", sep="\t", index=False
    )
    return root


def test_the_comparison_covers_every_condition(results_dir):
    table = exchange_comparison(results_dir)
    assert list(table["condition"]) == list(CONDITION_ORDER)


def test_comparing_a_run_against_itself_gives_no_difference(results_dir):
    """Fed the legacy fluxes as though they were new ones, the two
    columns must agree -- otherwise the comparison is reading one of the
    two formats wrongly."""
    table = exchange_comparison(results_dir)
    for quantity in ("glucose", "CO2", "oxygen", "ethanol"):
        assert table[f"{quantity}_gecko2"].to_numpy() == pytest.approx(
            table[f"{quantity}_gecko4"].to_numpy(), abs=1e-6
        )
    assert table["residual_gecko2"].to_numpy() == pytest.approx(
        table["residual_gecko4"].to_numpy()
    )


def test_uptakes_are_compared_as_magnitudes(results_dir):
    """GECKO 2 splits reversible reactions, so its glucose uptake is a
    positive flux on r_1714_REV while this pipeline's is a negative flux
    on r_1714. Comparing them unsigned is the whole point."""
    table = exchange_comparison(results_dir)
    assert (table["glucose_gecko2"] > 0).all()
    assert (table["glucose_gecko4"] > 0).all()


def test_usage_comparison_lines_the_systems_up(results_dir):
    table = usage_comparison(results_dir)
    assert "ETC" in set(table["system"])
    for condition in CONDITION_ORDER:
        assert table[f"{condition}_gecko2"].to_numpy() == pytest.approx(
            table[f"{condition}_gecko4"].to_numpy(), nan_ok=True
        )


def test_the_report_names_what_it_compares(results_dir):
    text = report(results_dir)
    assert "Exchange rates" in text
    assert "Mean residual" in text
    for condition in CONDITION_ORDER:
        assert condition in text


def test_the_sampled_budget_is_compared_when_a_run_is_present(results_dir):
    """The ATP budget lines up only if both files are read with the same
    row names; a mismatch would silently produce an empty table."""
    import shutil

    from overflow.compare import BUDGET_ROWS, budget_comparison

    sampling = results_dir / "randomSampling"
    sampling.mkdir(exist_ok=True)
    legacy = pd.read_csv(LEGACY / "randomSampling" / "selectedFluxes.txt", sep="\t")
    legacy.to_csv(sampling / "selectedFluxes.tsv", sep="\t", index=False)

    table = budget_comparison(results_dir)
    assert set(table["row"]) == set(BUDGET_ROWS)
    for condition in CONDITION_ORDER:
        assert table[f"{condition}_gecko2"].to_numpy() == pytest.approx(
            table[f"{condition}_gecko4"].to_numpy()
        )
    shutil.rmtree(sampling)
