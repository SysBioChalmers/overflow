import pytest

import overflow.solver as solver


class FakeConfiguration:
    requested = []

    @property
    def solver(self):
        return self.requested[-1] if self.requested else None

    @solver.setter
    def solver(self, name):
        if name == "missing":
            raise ValueError("Solver not available")
        self.requested.append(name)


@pytest.fixture
def fake(monkeypatch):
    FakeConfiguration.requested = []
    monkeypatch.setattr(solver.cobra, "Configuration", FakeConfiguration)
    monkeypatch.delenv("OVERFLOW_SOLVER", raising=False)
    return FakeConfiguration


def test_gurobi_is_the_default(fake):
    assert solver.DEFAULT_SOLVER == "gurobi"
    assert solver.use_solver() == "gurobi"
    assert fake.requested == ["gurobi"]


def test_the_environment_overrides_the_default(fake, monkeypatch):
    monkeypatch.setenv("OVERFLOW_SOLVER", "glpk")
    assert solver.use_solver() == "glpk"


def test_a_named_solver_beats_the_environment(fake, monkeypatch):
    monkeypatch.setenv("OVERFLOW_SOLVER", "glpk")
    assert solver.use_solver("gurobi") == "gurobi"


def test_an_unavailable_solver_is_an_error_not_a_fallback(fake):
    with pytest.raises(RuntimeError, match="cannot use the missing solver"):
        solver.use_solver("missing")
    assert fake.requested == []
