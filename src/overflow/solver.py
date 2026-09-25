"""The LP solver the analysis runs on: Gurobi unless told otherwise."""
import os
from typing import Optional

import cobra

DEFAULT_SOLVER = "gurobi"
HELP = "LP solver to use, e.g. gurobi or glpk (default: gurobi, or $OVERFLOW_SOLVER)"


def use_solver(name: Optional[str] = None) -> str:
    """Make ``name`` cobra's solver, and return it.

    Without a name it is ``$OVERFLOW_SOLVER``, else Gurobi. Asking for a solver that
    is not available is an error rather than a quiet switch to another one.
    """
    name = name or os.environ.get("OVERFLOW_SOLVER") or DEFAULT_SOLVER
    try:
        cobra.Configuration().solver = name
    except Exception as error:
        raise RuntimeError(
            f"cannot use the {name} solver ({error}). Gurobi is the default: install "
            "gurobipy with a licence, or choose another solver with --solver or "
            "OVERFLOW_SOLVER."
        ) from error
    return name
