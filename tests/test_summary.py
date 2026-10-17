"""The build summary keeps the rows of conditions built separately."""
import pandas as pd

from overflow.build import update_summary


def rows(*conditions, growth=0.1):
    return pd.DataFrame({"condition": list(conditions), "growth": growth})


def test_conditions_built_one_at_a_time_all_end_up_in_the_summary(tmp_path):
    path = tmp_path / "modelGeneration" / "summary.tsv"
    for name in ("hGR", "CN22", "CN4"):
        update_summary(path, rows(name))
    table = pd.read_csv(path, sep="\t")
    assert list(table["condition"]) == ["CN4", "CN22", "hGR"]


def test_a_rebuilt_condition_replaces_its_row(tmp_path):
    path = tmp_path / "summary.tsv"
    update_summary(path, rows("CN4", "CN22", growth=0.1))
    table = update_summary(path, rows("CN22", growth=0.2))
    assert list(table["condition"]) == ["CN4", "CN22"]
    assert list(table["growth"]) == [0.1, 0.2]


def test_conditions_outside_the_order_follow_it(tmp_path):
    path = tmp_path / "summary.tsv"
    update_summary(path, rows("mu010"))
    table = update_summary(path, rows("CN75"))
    assert list(table["condition"]) == ["CN75", "mu010"]
