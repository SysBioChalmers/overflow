import pytest

from overflow import build_adapter
from overflow.config import BIO_RXN, C_SOURCE


def test_defaults_come_from_the_toml(adapter):
    assert adapter.params.bio_rxn == BIO_RXN
    assert adapter.params.c_source == C_SOURCE
    assert adapter.params.sigma == pytest.approx(1.0)
    assert adapter.params.conv_gem.is_file()


def test_condition_sets_protein_content_and_reference_growth(conditions):
    hgr = build_adapter(conditions["hGR"])
    assert hgr.params.p_tot == pytest.approx(0.41)
    assert hgr.params.gr_exp == pytest.approx(0.29)


def test_explicit_overrides_win(conditions):
    a = build_adapter(conditions["CN4"], p_tot=0.42, sigma=0.5)
    assert a.params.p_tot == pytest.approx(0.42)
    assert a.params.sigma == pytest.approx(0.5)
    assert a.params.gr_exp == pytest.approx(0.1)


def test_building_an_adapter_does_not_mutate_the_shared_defaults(adapter, conditions):
    before = adapter.params.p_tot
    build_adapter(conditions["CN75"])
    assert adapter.params.p_tot == before
