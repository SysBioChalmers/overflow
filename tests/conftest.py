import os

import cobra
import pytest

# Tests run on the open solver by default: they are small enough not to
# need a commercial one, and a licence that is busy or absent should not
# turn into a wall of failures that look like code.
cobra.Configuration().solver = os.environ.get("OVERFLOW_TEST_SOLVER", "glpk")

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
