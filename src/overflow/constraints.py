"""Applying the measured rates to a model.

The constraints follow the original analysis: the carbon source is
allowed a little headroom over its measurement, a byproduct that was
detected may exceed its measurement by 10% and one that was not is
blocked, and the chemostat is set by fixing growth and then minimising
first uptake and then protein.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import cobra

from overflow.config import (
    BIO_RXN,
    C_SOURCE,
    CO2_RXN,
    NGAM_RXN,
    O2_RXN,
    POOL_RXN,
    Condition,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from geckopy import EcModel


class InfeasibleCondition(RuntimeError):
    """The measured rates leave the model with no feasible solution."""


def constrain_byproducts(model: cobra.Model, condition: Condition, flex: float = 1.1) -> None:
    """Bound the overflow metabolite exchanges to the measurement.

    A byproduct that was not detected is blocked, which is what keeps the
    model from paying for it in the conditions where it was absent.
    """
    for reaction_id, upper in condition.byproduct_bounds(flex).items():
        model.reactions.get_by_id(reaction_id).bounds = (0.0, upper)


def constrain_uptake(model: cobra.Model, condition: Condition, flex: float = 1.05) -> None:
    """Cap glucose uptake at ``flex`` times the measured rate."""
    model.reactions.get_by_id(C_SOURCE).bounds = (-flex * condition.glucose, 0.0)


def block_reactions(model: cobra.Model, reaction_ids) -> None:
    """Set the bounds of each named reaction to zero.

    An unknown ID is an error rather than a reaction that quietly stays open.
    """
    for reaction_id in reaction_ids:
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            raise KeyError(f"cannot block {reaction_id}: the model has no such reaction") from None
        reaction.bounds = (0.0, 0.0)


def free_ngam(model: cobra.Model, lower: float = 0.0) -> None:
    """Release the maintenance reaction so it can be fitted.

    The distributed model pins it, which would otherwise fix maintenance
    at a value belonging to another condition.
    """
    model.reactions.get_by_id(NGAM_RXN).bounds = (lower, 1000.0)


def set_chemostat_constraints(
    model: "EcModel",
    condition: Condition,
    flex: float = 0.01,
    glucose: Optional[float] = None,
    minimise_protein: bool = True,
) -> float:
    """Fix the chemostat and leave the model minimising protein usage.

    Growth is held at ``(1 - flex)`` times the dilution rate. Uptake is
    minimised first, then pinned to the value found, so that the protein
    minimisation that follows cannot buy protein with extra carbon.
    Passing ``glucose`` additionally requires uptake to reach
    ``(1 - flex)`` times that measurement.

    Returns the glucose uptake the model settled on, as a positive rate.
    """
    biomass = model.reactions.get_by_id(BIO_RXN)
    biomass.lower_bound = (1 - flex) * condition.d_rate

    carbon = model.reactions.get_by_id(C_SOURCE)
    carbon.bounds = (-1000.0, 0.0)
    if glucose is not None:
        carbon.upper_bound = -(1 - flex) * glucose

    # Uptake is a negative flux, so the least uptake is the largest flux.
    model.objective = C_SOURCE
    model.objective.direction = "max"
    solution = model.optimize()
    if solution.status != "optimal":
        raise InfeasibleCondition(
            f"{condition.name}: no solution at growth "
            f">= {(1 - flex) * condition.d_rate:.4f} /h"
        )

    uptake = -float(solution.fluxes[C_SOURCE])
    if minimise_protein:
        carbon.bounds = (-(1 + flex) * uptake, -uptake)
        # Written as maximising the negative rather than minimising,
        # because a model file records objective coefficients but not
        # the direction: read back, this still minimises protein, where
        # a plain minimisation would come back maximising it.
        model.objective = {model.reactions.get_by_id(POOL_RXN): -1.0}
        model.objective.direction = "max"
    return uptake


def constrain_measured_rates(
    model: cobra.Model,
    condition: Condition,
    tolerance: float = 0.05,
    minimum_secretion: float = 0.01,
) -> dict[str, tuple[float, float]]:
    """Hold every measured rate inside a band around its measurement.

    The pipeline otherwise constrains only uptake, growth and the
    byproducts, leaving the model free to dispose of carbon however is
    cheapest in protein. Banding the gas rates as well makes the
    measurements something the model has to account for rather than
    something it is merely allowed to approach.
    """
    from overflow.config import BYPRODUCT_RXNS

    bounds: dict[str, tuple[float, float]] = {}
    bounds[C_SOURCE] = (-(1 + tolerance) * condition.glucose,
                        -(1 - tolerance) * condition.glucose)
    bounds[CO2_RXN] = ((1 - tolerance) * condition.co2,
                       (1 + tolerance) * condition.co2)
    bounds[O2_RXN] = (-(1 + tolerance) * condition.oxygen,
                      -(1 - tolerance) * condition.oxygen)
    for name, reaction_id in BYPRODUCT_RXNS.items():
        value = condition.byproducts[name]
        if value == 0:
            bounds[reaction_id] = (0.0, minimum_secretion)
        else:
            bounds[reaction_id] = ((1 - tolerance) * value, (1 + tolerance) * value)

    for reaction_id, (low, high) in bounds.items():
        model.reactions.get_by_id(reaction_id).bounds = (low, high)
    return bounds
