import pytest

from overflow import build_adapter, load_conditions, load_gem, load_model


@pytest.fixture(scope="session")
def conditions():
    return load_conditions()


@pytest.fixture(scope="session")
def adapter():
    return build_adapter()


@pytest.fixture(scope="session")
def ec_model():
    """The distributed ecModel. Loaded once; do not mutate it in place."""
    return load_model()


@pytest.fixture(scope="session")
def conv_model():
    return load_gem()
