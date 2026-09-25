import math

import pandas as pd
import pytest

from overflow.config import (
    BYPRODUCT_RXNS,
    CONDITION_ORDER,
    FERMENTATION_DATA,
    PROTEOMICS_DATA,
    REPLICATE_PREFIX,
    load_conditions,
)


def test_all_conditions_present_and_ordered(conditions):
    assert tuple(conditions) == CONDITION_ORDER


def test_cn4_matches_the_data_file(conditions):
    cn4 = conditions["CN4"]
    assert cn4.p_tot == pytest.approx(0.5)
    assert cn4.d_rate == pytest.approx(0.1)
    assert cn4.glucose == pytest.approx(1.115503298)
    assert cn4.co2 == pytest.approx(3.087765343)
    assert cn4.oxygen == pytest.approx(2.928264928)


def test_undetected_byproducts_are_zero_and_not_measured(conditions):
    cn4 = conditions["CN4"]
    assert cn4.measured == frozenset()
    assert all(v == 0.0 for v in cn4.byproducts.values())

    cn22 = conditions["CN22"]
    assert cn22.measured == {"Formate"}
    assert cn22.byproducts["Formate"] == pytest.approx(0.048607541)
    assert cn22.byproducts["Ethanol"] == 0.0


def test_byproduct_bounds_block_what_was_not_detected(conditions):
    bounds = conditions["CN22"].byproduct_bounds(flex=1.1)
    assert bounds[BYPRODUCT_RXNS["Formate"]] == pytest.approx(1.1 * 0.048607541)
    assert bounds[BYPRODUCT_RXNS["Ethanol"]] == 0.0


def test_hgr_is_the_fermentative_condition(conditions):
    hgr = conditions["hGR"]
    assert hgr.d_rate == pytest.approx(0.29)
    assert hgr.measured == {"Acetate", "Ethanol"}


def test_replicate_prefixes_match_the_proteomics_columns():
    """hGR's replicate columns are named CN_hGR_*, unlike every other
    condition. A prefix that matches nothing would silently produce an
    empty abundance vector."""
    columns = pd.read_csv(PROTEOMICS_DATA, sep="\t", nrows=0).columns
    counts = {
        name: sum(c.startswith(prefix) for c in columns)
        for name, prefix in REPLICATE_PREFIX.items()
    }
    assert counts == {"CN4": 3, "CN22": 3, "CN38": 3, "CN75": 3, "hGR": 4}


def test_prefixes_do_not_overlap():
    """CN4_ must not also match CN40_-style columns of another condition."""
    columns = pd.read_csv(PROTEOMICS_DATA, sep="\t").columns
    claimed: dict[str, str] = {}
    for name, prefix in REPLICATE_PREFIX.items():
        for column in columns:
            if column.startswith(prefix):
                assert column not in claimed, (
                    f"{column} claimed by both {claimed.get(column)} and {name}"
                )
                claimed[column] = name
    assert len(claimed) == 16


def test_conditions_cover_every_row_of_the_data_file():
    table = pd.read_csv(FERMENTATION_DATA, sep="\t")
    assert set(table["Condition"]) == set(CONDITION_ORDER)
    assert not any(math.isnan(v) for v in table["Ptot"])
