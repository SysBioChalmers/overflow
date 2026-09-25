"""Fitting non-growth associated maintenance to the measured gas rates.

The maintenance requirement is the one parameter left free once growth
and uptake are fixed: it is what the model spends carbon on without
producing biomass, and it shows up in the CO2 and oxygen rates.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import numpy as np
from cobra.exceptions import OptimizationError
from cobra.util import solver as solver_util
from optlang.symbolics import Zero

from overflow.config import C_SOURCE, CO2_RXN, NGAM_RXN, O2_RXN, Condition

if TYPE_CHECKING:  # pragma: no cover - typing only
    from geckopy import EcModel


@dataclass
class NgamFit:
    value: float
    error: float
    at_upper_bound: bool
    scanned: np.ndarray
    errors: np.ndarray


def relative_residual(predicted: np.ndarray, measured: np.ndarray) -> float:
    """Root sum of squared relative deviations.

    Relative, so that a rate which is large in absolute terms does not
    dominate the fit simply by being large.
    """
    predicted = np.asarray(predicted, dtype=float)
    measured = np.asarray(measured, dtype=float)
    return float(np.sqrt(np.sum(((predicted - measured) / measured) ** 2)))


def parsimonious_solution(model: "EcModel", reactions: Optional[list] = None):
    """Minimise total flux at the current optimum, over ``reactions`` or all of them.

    ``cobra.flux_analysis.pfba`` minimises every reaction whatever it is
    asked to return; here only ``reactions`` enter the objective, and the
    solution still carries every reaction's flux.
    """
    reactions = model.reactions if reactions is None else reactions
    with model:
        solver_util.fix_objective_as_constraint(model, fraction=1.0)
        model.objective = model.problem.Objective(
            Zero, direction="min", sloppy=True, name="_parsimony_objective"
        )
        model.objective.set_linear_coefficients(
            {
                variable: 1.0
                for reaction in reactions
                for variable in (reaction.forward_variable, reaction.reverse_variable)
            }
        )
        return model.optimize()


def fit_ngam(
    model: "EcModel",
    condition: Condition,
    bounds: tuple[float, float] = (0.0, 5.0),
    steps: int = 100,
    ngam_rxn: str = NGAM_RXN,
    targets: tuple[str, str, str] = (C_SOURCE, CO2_RXN, O2_RXN),
    parsimony_reactions: Optional[list] = None,
) -> NgamFit:
    """Scan the maintenance requirement for the best fit to the gas rates.

    The residual is relative, over glucose uptake, CO2 production and
    oxygen uptake together, so none of the three dominates by magnitude.
    The model keeps whatever objective it was given; each point is
    evaluated parsimoniously so that the gas rates are determined.

    The fitted value becomes the reaction's lower bound, as maintenance
    is a requirement rather than a fixed expense.
    """
    reaction = model.reactions.get_by_id(ngam_rxn)
    reaction.bounds = (bounds[0], 1000.0)

    measured = np.array([condition.glucose, condition.co2, condition.oxygen], float)

    scanned = np.array([bounds[0] + (bounds[1] - bounds[0]) * i / steps for i in range(steps)])
    errors = np.full(steps, np.inf)

    for index, value in enumerate(scanned):
        reaction.lower_bound = value
        try:
            solution = parsimonious_solution(model, parsimony_reactions)
        except OptimizationError:
            continue
        if solution.status != "optimal":
            continue
        predicted = np.abs(np.array([solution.fluxes[r] for r in targets], float))
        errors[index] = relative_residual(predicted, measured)

    best = int(np.argmin(errors))
    at_upper = best == steps - 1
    value = bounds[0] if at_upper else float(scanned[best])
    reaction.lower_bound = value
    return NgamFit(
        value=value,
        error=float(errors[best]),
        at_upper_bound=at_upper,
        scanned=scanned,
        errors=errors,
    )
