"""Reconciling measured enzyme abundances with what the model needs.

Two passes, in the order the original analysis applies them:

1. A minimum-usage pass on a model constrained only by the protein pool.
   An enzyme the model cannot run at the measured dilution rate for less
   than its measured abundance has that abundance raised to what the
   model needs. This catches enzymes whose kcat is too low for the
   measurement, before any of them are allowed to limit growth.
2. Iterative flexibilization of whatever still limits growth. The most
   limiting enzyme is released outright, one at a time, until the model
   reaches the dilution rate; if nothing single is limiting, the protein
   pool itself is grown. The released enzymes are then re-tightened to
   the usage they take in a minimum-protein solution, so releasing an
   enzyme costs only what the model actually spends on it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

import numpy as np
import pandas as pd
from geckopy import constrain_enz_concs, get_conc_control_coeffs

from overflow.config import BIO_RXN, POOL_RXN

if TYPE_CHECKING:  # pragma: no cover - typing only
    from geckopy import EcModel

USAGE_PREFIX = "usage_prot_"
#: Upper bound that marks an enzyme as unconstrained.
FREE = 1000.0


@dataclass
class MinimumUsageResult:
    """Outcome of the minimum-usage pass."""

    concentrations: np.ndarray
    raised: list[str] = field(default_factory=list)
    previous: list[float] = field(default_factory=list)
    minimum: list[float] = field(default_factory=list)

    def table(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "protein": self.raised,
                "measured_mg_gDW": self.previous,
                "required_mg_gDW": self.minimum,
                "fold_change": [
                    m / p if p else np.inf for p, m in zip(self.previous, self.minimum)
                ],
            }
        )


@dataclass
class FlexibilizationResult:
    """Outcome of the iterative flexibilization."""

    released: list[str] = field(default_factory=list)
    previous: list[float] = field(default_factory=list)
    modified: list[float] = field(default_factory=list)
    pool_increases: int = 0
    pool_before: float = 0.0
    pool_after: float = 0.0
    growth: float = 0.0
    reached_target: bool = False

    def table(self) -> pd.DataFrame:
        previous = np.asarray(self.previous, dtype=float)
        modified = np.asarray(self.modified, dtype=float)
        return pd.DataFrame(
            {
                "protein": self.released,
                "previous_mg_gDW": previous,
                "modified_mg_gDW": modified,
                "flex_mass_mg_gDW": modified - previous,
                "fold_change": np.divide(
                    modified, previous, out=np.full_like(modified, np.inf),
                    where=previous > 0,
                ),
            }
        )


def measured_enzymes(model: "EcModel") -> list[int]:
    """Indices in ``model.ec.enzymes`` that carry a measured concentration."""
    return [i for i, c in enumerate(model.ec.concs) if not np.isnan(c)]


def minimum_usage_pass(
    model: "EcModel",
    target_growth: float,
    tolerance: float = 1.01,
    bio_rxn: str = BIO_RXN,
) -> MinimumUsageResult:
    """Raise abundances the model cannot run the condition on.

    ``model`` must already be constrained to the condition and limited by
    the protein pool only -- individual enzyme caps would make each
    minimum a function of the others. Its bounds are left untouched;
    only ``ec.concs`` is returned, updated.
    """
    concentrations = np.asarray(model.ec.concs, dtype=float).copy()
    result = MinimumUsageResult(concentrations=concentrations)

    biomass = model.reactions.get_by_id(bio_rxn)
    original_bounds = biomass.bounds
    biomass.lower_bound = target_growth
    try:
        for index in measured_enzymes(model):
            enzyme = model.ec.enzymes[index]
            reaction_id = f"{USAGE_PREFIX}{enzyme}"
            with model:
                model.objective = reaction_id
                model.objective.direction = "min"
                solution = model.optimize()
                usage = (
                    float(solution.fluxes[reaction_id])
                    if solution.status == "optimal"
                    else None
                )
            if usage is not None and usage > concentrations[index]:
                result.raised.append(enzyme)
                result.previous.append(float(concentrations[index]))
                result.minimum.append(tolerance * usage)
                concentrations[index] = tolerance * usage
    finally:
        biomass.bounds = original_bounds

    result.concentrations = concentrations
    return result


def flexibilize_proteins(
    model: "EcModel",
    target_growth: float,
    max_iterations: int = 500,
    pool_step: float = 1.01,
    growth_tolerance: float = 1e-4,
    bound_tolerance: float = 1e-6,
    fold_change: float = 100.0,
    patience: int = 25,
    bio_rxn: str = BIO_RXN,
    pool_rxn: str = POOL_RXN,
) -> FlexibilizationResult:
    """Release limiting enzymes until the model reaches ``target_growth``.

    The model is expected to carry its condition constraints and its
    enzyme caps. Growth is maximised throughout; the caller decides what
    the objective should be afterwards.
    """
    pool = model.reactions.get_by_id(pool_rxn)
    result = FlexibilizationResult(pool_before=pool.upper_bound)

    biomass = model.reactions.get_by_id(bio_rxn)
    biomass.lower_bound = 0.0
    model.objective = bio_rxn
    model.objective.direction = "max"

    original: dict[str, float] = {}
    threshold = target_growth * (1 - growth_tolerance)

    best_growth = -np.inf
    stalled = 0
    for _ in range(max_iterations):
        solution = model.optimize()
        growth = solution.objective_value if solution.status == "optimal" else 0.0
        growth = float(growth or 0.0)
        if growth >= threshold:
            break

        if growth > best_growth + growth_tolerance * target_growth:
            best_growth, stalled = growth, 0
        else:
            stalled += 1
        if stalled > patience:
            break

        enzyme = _most_limiting(model, solution, bound_tolerance, fold_change)
        if enzyme is not None:
            reaction = model.reactions.get_by_id(f"{USAGE_PREFIX}{enzyme}")
            original.setdefault(enzyme, reaction.upper_bound)
            reaction.upper_bound = FREE
            result.released.append(enzyme)
        else:
            pool.upper_bound *= pool_step
            result.pool_increases += 1

    solution = model.optimize()
    result.growth = float(solution.objective_value or 0.0)
    result.reached_target = result.growth >= threshold
    result.pool_after = pool.upper_bound

    _retighten(model, target_growth, original, result, bio_rxn, pool_rxn)
    return result


def binding_caps(model: "EcModel", solution, bound_tolerance: float) -> list[str]:
    """Enzymes sitting at a cap that is not already released.

    Only an enzyme at its own cap can be limiting through that cap. An
    enzyme that merely carries flux also has a positive shadow price
    whenever the protein pool binds -- the price of the pool, not of the
    cap -- so considering those would keep re-releasing enzymes that are
    already free while the real constraint went untouched.
    """
    candidates = []
    for enzyme in model.ec.enzymes:
        reaction = model.reactions.get_by_id(f"{USAGE_PREFIX}{enzyme}")
        if reaction.upper_bound >= FREE:
            continue
        usage = float(solution.fluxes[reaction.id])
        if usage >= (1 - bound_tolerance) * reaction.upper_bound:
            candidates.append(enzyme)
    return candidates


def _most_limiting(
    model: "EcModel",
    solution,
    bound_tolerance: float,
    fold_change: float = 100.0,
) -> Optional[str]:
    """The capped enzyme worth the most growth, if any.

    Shadow prices answer this in one solve, but they go silent on a
    degenerate optimum: when several caps bind at once only one of them
    carries a non-zero dual, and sometimes none does. Falling back to
    probing each candidate cap directly costs one solve per candidate
    but cannot be fooled that way, which is what keeps a degenerate
    model from looking as though nothing at all were limiting it.
    """
    candidates = binding_caps(model, solution, bound_tolerance)
    if not candidates:
        return None

    _, coefficients = get_conc_control_coeffs(model, proteins=candidates)
    if coefficients.size and coefficients.max() > 0:
        return candidates[int(np.argmax(coefficients))]

    baseline = float(solution.objective_value or 0.0)
    best, best_gain = None, 0.0
    for enzyme in candidates:
        reaction = model.reactions.get_by_id(f"{USAGE_PREFIX}{enzyme}")
        previous = reaction.upper_bound
        reaction.upper_bound = min(previous * fold_change, FREE) if previous else FREE
        try:
            probe = model.optimize()
        finally:
            reaction.upper_bound = previous
        if probe.status != "optimal":
            continue
        gain = float(probe.objective_value or 0.0) - baseline
        if gain > best_gain:
            best, best_gain = enzyme, gain
    return best


def _retighten(
    model: "EcModel",
    target_growth: float,
    original: dict[str, float],
    result: FlexibilizationResult,
    bio_rxn: str = BIO_RXN,
    pool_rxn: str = POOL_RXN,
) -> None:
    """Pull the released enzymes back to a minimum-protein solution.

    Released enzymes keep their previous cap where the minimum-protein
    solution does not use them at all, so a release that turned out to be
    unnecessary leaves no trace.
    """
    biomass = model.reactions.get_by_id(bio_rxn)
    biomass.lower_bound = 0.9999 * target_growth
    model.objective = pool_rxn
    model.objective.direction = "min"
    solution = model.optimize()

    result.previous = [original[e] for e in result.released]
    if solution.status != "optimal":
        result.modified = list(result.previous)
        for enzyme, previous in zip(result.released, result.previous):
            model.reactions.get_by_id(f"{USAGE_PREFIX}{enzyme}").upper_bound = previous
        return

    modified: list[float] = []
    for enzyme, previous in zip(result.released, result.previous):
        reaction = model.reactions.get_by_id(f"{USAGE_PREFIX}{enzyme}")
        usage = float(solution.fluxes[reaction.id])
        new_bound = usage if usage > 0 else previous
        reaction.upper_bound = new_bound
        modified.append(new_bound)
    result.modified = modified


def write_back_concentrations(model: "EcModel") -> None:
    """Record the enzyme caps in ``ec.concs`` so they survive saving."""
    concentrations = np.asarray(model.ec.concs, dtype=float).copy()
    for index, enzyme in enumerate(model.ec.enzymes):
        if np.isnan(concentrations[index]):
            continue
        bound = model.reactions.get_by_id(f"{USAGE_PREFIX}{enzyme}").upper_bound
        concentrations[index] = bound
    model.ec.concs = concentrations


def apply_concentrations(model: "EcModel", concentrations: np.ndarray) -> None:
    """Write concentrations into the model and its enzyme caps."""
    model.ec.concs = np.asarray(concentrations, dtype=float)
    constrain_enz_concs(model)


@dataclass
class RateFitResult:
    """Outcome of relaxing abundances until the measured rates fit."""

    proteins: list[str] = field(default_factory=list)
    measured: list[float] = field(default_factory=list)
    relaxed: list[float] = field(default_factory=list)
    added_mg: float = 0.0
    feasible: bool = False
    status: str = ""

    def table(self) -> pd.DataFrame:
        measured = np.asarray(self.measured, dtype=float)
        relaxed = np.asarray(self.relaxed, dtype=float)
        return pd.DataFrame(
            {
                "protein": self.proteins,
                "measured_mg_gDW": measured,
                "relaxed_mg_gDW": relaxed,
                "added_mg_gDW": relaxed - measured,
                "fold_change": np.divide(
                    relaxed, measured, out=np.full_like(relaxed, np.inf),
                    where=measured > 0,
                ),
            }
        )


def relax_to_measured_rates(
    model: "EcModel",
    target_growth: float,
    weight: str = "relative",
    bio_rxn: str = BIO_RXN,
    floor: float = 1e-9,
) -> RateFitResult:
    """Raise measured abundances by the least that fits the measured rates.

    The model is expected to carry its rate constraints already. Each
    measured enzyme's cap becomes ``usage <= measured + slack``, and the
    slacks are minimised, so the answer is the smallest departure from
    the proteomics that accounts for the rates rather than one enzyme's
    cap released outright.

    ``weight`` decides what "smallest" means: ``"relative"`` minimises
    the sum of fold-increases, so an enzyme measured at almost nothing
    is not relaxed freely just because its absolute cost is small;
    ``"absolute"`` minimises the milligrams added.
    """
    from optlang.symbolics import add

    if weight not in ("relative", "absolute"):
        raise ValueError(f"unknown weight {weight!r}")

    biomass = model.reactions.get_by_id(bio_rxn)
    biomass.bounds = (target_growth, biomass.upper_bound or 1000.0)

    # The slack objective is this function's own business; the caller's
    # objective has to survive it, or the model comes back optimising
    # an expression whose variables have just been removed.
    previous_objective = model.objective.expression
    previous_direction = model.objective.direction

    indices = measured_enzymes(model)
    slacks, added = {}, []
    for index in indices:
        enzyme = model.ec.enzymes[index]
        reaction = model.reactions.get_by_id(f"{USAGE_PREFIX}{enzyme}")
        concentration = float(model.ec.concs[index])
        reaction.upper_bound = FREE
        slack = model.problem.Variable(f"relax_{enzyme}", lb=0.0, ub=FREE)
        constraint = model.problem.Constraint(
            reaction.flux_expression - slack, ub=concentration, name=f"cap_{enzyme}"
        )
        slacks[enzyme] = (slack, concentration)
        added.extend([slack, constraint])
    model.add_cons_vars(added)

    try:
        model.objective = model.problem.Objective(
            add([
                (1.0 / max(concentration, floor) if weight == "relative" else 1.0) * slack
                for slack, concentration in slacks.values()
            ]),
            direction="min",
        )
        solution = model.optimize()
        result = RateFitResult(status=solution.status)
        result.feasible = solution.status == "optimal"
        if result.feasible:
            values = model.solver.primal_values
            concentrations = np.asarray(model.ec.concs, dtype=float).copy()
            for enzyme, (slack, concentration) in slacks.items():
                raised = concentration + float(values.get(slack.name, 0.0))
                if raised > concentration * (1 + 1e-9) + 1e-12:
                    result.proteins.append(enzyme)
                    result.measured.append(concentration)
                    result.relaxed.append(raised)
                    result.added_mg += raised - concentration
                concentrations[model.ec.enzymes.index(enzyme)] = raised
            model.ec.concs = concentrations
    finally:
        model.remove_cons_vars(added)
        model.objective = model.problem.Objective(
            previous_objective, direction=previous_direction
        )

    if result.feasible:
        constrain_enz_concs(model)
    return result
