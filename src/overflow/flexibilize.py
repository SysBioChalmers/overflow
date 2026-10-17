"""Reconciling measured enzyme abundances with what the model needs.

Three steps, each treating an enzyme complex as one unit
(:func:`overflow.complexes.complex_families`): when a complex needs more
capacity all its measured subunits are raised by the same factor, so a
complex levelled to one abundance stays levelled.

1. A minimum-usage pass on a model constrained only by the protein pool.
   A complex the model cannot run at the dilution rate for less than its
   measured abundance is raised to what the model needs. This catches
   enzymes whose kcat is too low for the measurement, before any of them
   are allowed to limit growth.
2. The least relaxation that reaches the dilution rate
   (:func:`relax_to_growth`): one LP minimises the sum of the fold
   increases of the measured caps, and of the protein pool at a far
   higher price, so the answer does not depend on the order in which
   limiting caps are found.
3. With measured rates, the least further relaxation that reproduces them
   (:func:`relax_to_measured_rates`), the same LP under the rate
   constraints.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

import numpy as np
import pandas as pd
from geckopy import constrain_enz_concs

from overflow.complexes import complex_families
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
    """Outcome of the relaxation to the dilution rate."""

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
    usages: dict[int, float] = {}
    try:
        for index in measured_enzymes(model):
            reaction_id = f"{USAGE_PREFIX}{model.ec.enzymes[index]}"
            with model:
                model.objective = reaction_id
                model.objective.direction = "min"
                solution = model.optimize()
                if solution.status == "optimal":
                    usages[index] = float(solution.fluxes[reaction_id])
    finally:
        biomass.bounds = original_bounds

    families = complex_families(model)
    position = {e: i for i, e in enumerate(model.ec.enzymes)}
    done: set[int] = set()
    for index in measured_enzymes(model):
        if index in done:
            continue
        group = [position[m] for m in families.members(model.ec.enzymes[index])
                 if m in position and not np.isnan(concentrations[position[m]])]
        done.update(group)
        if not any(usages.get(i, 0.0) > concentrations[i] for i in group):
            continue
        fold = max((tolerance * usages[i] / concentrations[i]
                    for i in group if i in usages and concentrations[i] > 0), default=1.0)
        for i in group:
            raised = concentrations[i] * fold if concentrations[i] > 0 else tolerance * usages.get(i, 0.0)
            if raised > concentrations[i]:
                result.raised.append(model.ec.enzymes[i])
                result.previous.append(float(concentrations[i]))
                result.minimum.append(float(raised))
                concentrations[i] = raised

    result.concentrations = concentrations
    return result


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

    See :func:`_minimal_relaxation`; the model carries its rate constraints.
    """
    return _minimal_relaxation(model, target_growth, weight, bio_rxn, floor)[0]


def relax_to_growth(
    model: "EcModel",
    target_growth: float,
    weight: str = "relative",
    pool_weight: float = 1000.0,
    bio_rxn: str = BIO_RXN,
    pool_rxn: str = POOL_RXN,
    floor: float = 1e-9,
) -> FlexibilizationResult:
    """Raise measured abundances, and if need be the pool, by the least that reaches ``target_growth``.

    One LP chooses every relaxation at once (see :func:`_minimal_relaxation`),
    so the result does not depend on the order in which limiting caps
    would be found. The protein pool can grow by a fold too, at
    ``pool_weight`` times the price of an enzyme's fold, so that it grows
    only where no relaxation of enzymes suffices.
    """
    pool = model.reactions.get_by_id(pool_rxn)
    result = FlexibilizationResult(pool_before=pool.upper_bound)
    fit, pool_fold = _minimal_relaxation(
        model, target_growth, weight, bio_rxn, floor, pool_weight=pool_weight, pool_rxn=pool_rxn,
    )
    result.released, result.previous, result.modified = fit.proteins, fit.measured, fit.relaxed
    result.pool_increases = int(pool_fold > 1e-9)
    result.pool_after = pool.upper_bound
    model.objective = bio_rxn
    model.objective.direction = "max"
    solution = model.optimize()
    result.growth = float(solution.objective_value or 0.0) if solution.status == "optimal" else 0.0
    result.reached_target = fit.feasible and result.growth >= target_growth * (1 - 1e-4)
    return result


def _minimal_relaxation(
    model: "EcModel",
    target_growth: float,
    weight: str = "relative",
    bio_rxn: str = BIO_RXN,
    floor: float = 1e-9,
    pool_weight: Optional[float] = None,
    pool_rxn: str = POOL_RXN,
) -> tuple[RateFitResult, float]:
    """The least relaxation of measured abundances that reaches ``target_growth``.

    Whatever else constrains the model (the measured rates, for the rate
    fit) stays in place. Each measured enzyme's cap becomes
    ``usage <= measured * (1 + x)``, with one ``x`` per complex (shared by
    its measured subunits; a subunit measured at zero is its own unit) or
    per enzyme, and the ``x`` are minimised, so the answer is the smallest
    departure from the proteomics rather than one enzyme's cap released
    outright. With ``pool_weight`` the protein pool may grow by a fold too,
    at that price per fold. Returns the result and the pool's fold.

    ``weight`` decides what "smallest" means: ``"relative"`` minimises
    the sum of fold-increases, so an enzyme measured at almost nothing
    is not relaxed freely just because its absolute cost is small;
    ``"absolute"`` minimises the milligrams added.
    """
    from optlang.symbolics import add

    if weight not in ("relative", "absolute"):
        raise ValueError(f"unknown weight {weight!r}")

    biomass = model.reactions.get_by_id(bio_rxn)
    biomass.bounds = (target_growth, max(biomass.upper_bound or 1000.0, target_growth))

    # The slack objective is this function's own business; the caller's
    # objective has to survive it, or the model comes back optimising
    # an expression whose variables have just been removed.
    previous_objective = model.objective.expression
    previous_direction = model.objective.direction

    indices = measured_enzymes(model)
    families = complex_families(model)
    position = {e: i for i, e in enumerate(model.ec.enzymes)}
    measured = set(indices)
    groups, done = [], set()
    for index in indices:
        if index in done:
            continue
        # a subunit measured at zero cannot share a fold with the others
        if not model.ec.concs[index] > 0:
            group = [index]
        else:
            group = [position[m] for m in families.members(model.ec.enzymes[index])
                     if m in position and position[m] in measured and position[m] not in done
                     and model.ec.concs[position[m]] > 0]
        done.update(group)
        groups.append(group)

    folds, added = [], []
    pool_x = 0.0
    for k, group in enumerate(groups):
        fold = model.problem.Variable(f"relax_{k}", lb=0.0, ub=None)
        added.append(fold)
        for index in group:
            enzyme = model.ec.enzymes[index]
            reaction = model.reactions.get_by_id(f"{USAGE_PREFIX}{enzyme}")
            concentration = float(model.ec.concs[index])
            reaction.upper_bound = FREE
            added.append(model.problem.Constraint(
                reaction.flux_expression - max(concentration, floor) * fold,
                ub=concentration, name=f"cap_{enzyme}",
            ))
        folds.append((fold, group))
    pool_fold = None
    if pool_weight is not None:
        pool = model.reactions.get_by_id(pool_rxn)
        pool_cap = pool.upper_bound
        pool_fold = model.problem.Variable("relax_pool", lb=0.0, ub=None)
        pool.upper_bound = FREE
        added += [pool_fold, model.problem.Constraint(
            pool.flux_expression - pool_cap * pool_fold, ub=pool_cap, name="cap_pool")]
    model.add_cons_vars(added)

    try:
        terms = [
            (1.0 if weight == "relative" else sum(max(float(model.ec.concs[i]), floor) for i in group)) * fold
            for fold, group in folds
        ]
        if pool_fold is not None:
            terms.append(pool_weight * (1.0 if weight == "relative" else pool_cap) * pool_fold)
        model.objective = model.problem.Objective(add(terms), direction="min")
        solution = model.optimize()
        result = RateFitResult(status=solution.status)
        result.feasible = solution.status == "optimal"
        if result.feasible:
            values = model.solver.primal_values
            concentrations = np.asarray(model.ec.concs, dtype=float).copy()
            for fold, group in folds:
                x = float(values.get(fold.name, 0.0))
                for index in group:
                    concentration = float(model.ec.concs[index])
                    raised = concentration + max(concentration, floor) * x
                    if raised > concentration * (1 + 1e-9) + 1e-12:
                        result.proteins.append(model.ec.enzymes[index])
                        result.measured.append(concentration)
                        result.relaxed.append(raised)
                        result.added_mg += raised - concentration
                    concentrations[index] = raised
            model.ec.concs = concentrations
            if pool_fold is not None:
                pool_x = float(values.get(pool_fold.name, 0.0))
    finally:
        model.remove_cons_vars(added)
        if pool_fold is not None:
            model.reactions.get_by_id(pool_rxn).upper_bound = pool_cap * (1 + pool_x) if result.feasible else pool_cap
        model.objective = model.problem.Objective(
            previous_objective, direction=previous_direction
        )

    if result.feasible:
        constrain_enz_concs(model)
    return result, (pool_x if pool_fold is not None and result.feasible else 0.0)
