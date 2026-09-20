"""Random sampling of the conventional model under the measured rates.

No enzyme constraints here: the question is what the measured exchange
rates alone imply about the flux distribution, so the model is the plain
yeast-GEM with its biomass rescaled to the condition and every measured
rate held inside a narrow band around its measurement.

Sampling draws vertices by maximising small random objectives, which is
what the published analysis did. It is not a uniform draw from the
interior; means over such draws are vertex means.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Mapping, Optional

import numpy as np
import pandas as pd
from raven_toolbox.analysis.sampling import random_sampling

from overflow.biomass import (
    GAM_NO_POLYMERIZATION,
    POLYMERIZATION_COST,
    YEAST_BIOMASS,
    apply_protein_content,
)
from overflow.config import BIO_RXN, GROWTH_RXN, NGAM_RXN, Condition

if TYPE_CHECKING:  # pragma: no cover - typing only
    import cobra

#: Exchange reactions held at their measurement, and the sign convention
#: each is measured with: -1 for an uptake, +1 for a secretion.
MEASURED_EXCHANGES = {
    "r_1714": ("Glucose", -1),
    "r_1992": ("Oxygen", -1),
    "r_1672": ("CO2", +1),
    "r_1808": ("Glycerol", +1),
    "r_1634": ("Acetate", +1),
    "r_1761": ("Ethanol", +1),
    "r_1793": ("Formate", +1),
    GROWTH_RXN: ("D", +1),
}

#: Formate dehydrogenase, blocked so that formate leaves the cell rather
#: than being oxidised to CO2 inside it.
FORMATE_DEHYDROGENASE = "r_0445"

FORMATE_EXCHANGE = "r_1793"


def use_fork_start_method() -> None:
    """Have multiprocessing fork rather than use a forkserver.

    Python 3.14 made forkserver the default on Linux. Its helper process
    does not survive in every environment -- on this cluster it dies and
    takes the run with it, as a broken pipe or a reset connection partway
    through sampling. Fork is what cobrapy's process pool was written
    against and works here.
    """
    import multiprocessing

    if "fork" not in multiprocessing.get_all_start_methods():
        return
    if multiprocessing.get_start_method(allow_none=True) == "fork":
        return
    try:
        multiprocessing.set_start_method("fork", force=True)
    except RuntimeError:  # already started a pool; nothing to do
        pass


def measured_rate(condition: Condition, field_name: str) -> float:
    """One measured rate, as a positive magnitude."""
    if field_name == "D":
        return condition.d_rate
    if field_name == "Glucose":
        return condition.glucose
    if field_name == "Oxygen":
        return condition.oxygen
    if field_name == "CO2":
        return condition.co2
    return condition.byproducts.get(field_name, 0.0)


def sampling_bounds(
    condition: Condition,
    tolerance: float = 0.05,
    minimum_secretion: float = 0.01,
    include_formate: bool = True,
) -> dict[str, tuple[float, float]]:
    """Bounds holding every measured rate within ``tolerance`` of itself.

    A byproduct that was not detected gets a small allowance rather than
    being blocked, so that the sampler can say whether the model wants to
    make it.
    """
    bounds: dict[str, tuple[float, float]] = {}
    for reaction_id, (field_name, sign) in MEASURED_EXCHANGES.items():
        if reaction_id == FORMATE_EXCHANGE and not include_formate:
            continue
        value = measured_rate(condition, field_name)
        low = (1 - tolerance) * value
        high = (1 + tolerance) * value
        if sign < 0:
            low, high = -high, -low
        else:
            high = high if high != 0 else minimum_secretion
        bounds[reaction_id] = (low, high)
    return bounds


def prepare_model(
    model: "cobra.Model",
    condition: Condition,
    gam_base: float = GAM_NO_POLYMERIZATION,
) -> dict[str, float]:
    """Set the condition's biomass and open maintenance.

    Biomass is rescaled to the measured protein content and maintenance
    rebuilt from it, which is what this analysis fixes the ATP budget
    against. Unlike the enzyme-constrained models there is no protein
    budget here to make that unaffordable.
    """
    composition = apply_protein_content(
        model, condition.p_tot, gam="polymerization", gam_base=gam_base,
        scale_protein=True,
    )
    model.reactions.get_by_id(FORMATE_DEHYDROGENASE).bounds = (0.0, 0.0)
    model.reactions.get_by_id(NGAM_RXN).bounds = (0.0, 1000.0)
    return composition


def apply_bounds(model: "cobra.Model", bounds: Mapping[str, tuple[float, float]]) -> None:
    for reaction_id, (low, high) in bounds.items():
        model.reactions.get_by_id(reaction_id).bounds = (low, high)


def constrain_maintenance(model: "cobra.Model", fraction: float = 0.95) -> float:
    """Hold maintenance near the most the condition can support.

    Everything else is fixed by the measured rates, so whatever ATP the
    model does not spend on growth it spends here; leaving maintenance
    free would let it dissipate that instead through whichever internal
    cycle happens to be cheapest.
    """
    maintenance = model.reactions.get_by_id(NGAM_RXN)
    maintenance.bounds = (0.0, 1000.0)
    model.objective = NGAM_RXN
    model.objective.direction = "max"
    solution = model.optimize()
    if solution.status != "optimal":
        raise RuntimeError("no feasible solution under the measured rates")
    highest = float(solution.fluxes[NGAM_RXN])
    maintenance.lower_bound = fraction * highest
    return highest


def polymerization_breakdown(model: "cobra.Model", gam_base: float = GAM_NO_POLYMERIZATION) -> dict[str, float]:
    """ATP spent polymerising each component of one gram of biomass."""
    from raven_toolbox.biomass.scale import sum_biomass

    fractions = sum_biomass(model, YEAST_BIOMASS)
    costs = {
        component: cost * fractions.get(component, 0.0)
        for component, cost in POLYMERIZATION_COST.items()
    }
    costs["total"] = sum(costs.values())
    costs["GAEC_noPol"] = gam_base
    costs["GAEC_total"] = costs["total"] + gam_base
    return costs


def loopless_bounds(model: "cobra.Model", processes: Optional[int] = None) -> pd.DataFrame:
    """Flux range of every reaction that does not need a closed cycle.

    Ordinary variability analysis lets a reaction in a thermodynamically
    infeasible cycle run to whatever arbitrary bound the model carries;
    the loop-free range is what the reaction can do while the rest of
    the network stays consistent.
    """
    from cobra.flux_analysis import flux_variability_analysis

    if processes is not None:
        import cobra

        cobra.Configuration().processes = processes
    return flux_variability_analysis(
        model, fraction_of_optimum=0.0, loopless="cycleFreeFlux"
    )


def apply_loopless_bounds(
    model: "cobra.Model", ranges: pd.DataFrame, tolerance: float = 1e-9
) -> int:
    """Tighten every reaction to its loop-free range.

    Every flux distribution free of closed cycles already lies inside
    these bounds, so nothing thermodynamically sensible is excluded --
    what goes is the room a cycle needs to inflate a flux. Returns how
    many reactions were tightened.
    """
    tightened = 0
    for reaction_id, row in ranges.iterrows():
        reaction = model.reactions.get_by_id(reaction_id)
        low = max(reaction.lower_bound, float(row["minimum"]) - tolerance)
        high = min(reaction.upper_bound, float(row["maximum"]) + tolerance)
        if low > high:  # numerical crossing; leave the reaction alone
            continue
        if low > reaction.lower_bound or high < reaction.upper_bound:
            reaction.bounds = (low, high)
            tightened += 1
    return tightened


@dataclass
class SamplingResult:
    """Sampled fluxes for one condition and one formate setting."""

    condition: str
    include_formate: bool
    reactions: list[str] = field(default_factory=list)
    mean: np.ndarray = field(default_factory=lambda: np.empty(0))
    sd: np.ndarray = field(default_factory=lambda: np.empty(0))
    n_samples: int = 0
    max_maintenance: float = 0.0
    loopless_tightened: int = 0

    def to_frame(self) -> pd.DataFrame:
        with np.errstate(divide="ignore", invalid="ignore"):
            relative = np.where(self.mean != 0, self.sd / self.mean * 100, np.nan)
        return pd.DataFrame(
            {"rxnID": self.reactions, "mean": self.mean, "sd": self.sd,
             "sd_percent": relative}
        )

    @property
    def means(self) -> pd.Series:
        return pd.Series(self.mean, index=self.reactions)


def sample_condition(
    model: "cobra.Model",
    condition: Condition,
    n_samples: int = 5000,
    include_formate: bool = True,
    tolerance: float = 0.05,
    seed: Optional[int] = None,
    good_reactions: Optional[list[str]] = None,
    n_proc: Optional[int] = None,
    min_flux: Optional[bool] = None,
    replace_max_bound: bool = False,
    loopless: bool = True,
) -> tuple[SamplingResult, list[str]]:
    """Sample one condition, returning the summary and the good reactions.

    ``good_reactions`` is the loop-free set the sampler may use as random
    objectives; it depends only on the network, so it is computed once
    and handed back for reuse.

    ``replace_max_bound`` opens the arbitrary 1000 bounds to infinity, as
    the published analysis did, so that a reaction in a loop cannot sit
    at 1000 and be mistaken for flux. It is off here: the screening that
    picks the random objectives runs on the finite bounds, so once they
    are opened an objective can turn out unbounded, and this sampler
    raises on that where RAVEN's returned no solution and moved on.

    ``min_flux`` minimises total flux within each draw. The published
    second pass did this to suppress loops; here one awkward draw in
    several thousand raises out of the parsimonious solve and takes the
    whole run with it, so it is off unless asked for.
    """
    apply_bounds(model, sampling_bounds(condition, tolerance, include_formate=include_formate))
    highest = constrain_maintenance(model)

    # Without this the malate dehydrogenases cycle against each other at
    # several hundred mmol/gDW/h, which is the arbitrary 1000 bound
    # showing through rather than anything the cell does.
    tightened = 0
    if loopless:
        tightened = apply_loopless_bounds(model, loopless_bounds(model, n_proc))

    sampled = random_sampling(
        model,
        n_samples,
        method="random_objective",
        n_objectives=2,
        replace_max_bound=replace_max_bound,
        suppress_errors=True,
        good_reactions=good_reactions,
        min_flux=include_formate if min_flux is None else min_flux,
        seed=seed,
        n_proc=n_proc,
    )
    samples = sampled.samples
    result = SamplingResult(
        condition=condition.name,
        include_formate=include_formate,
        reactions=list(samples.columns),
        mean=samples.mean(axis=0).to_numpy(),
        sd=samples.std(axis=0, ddof=1).to_numpy(),
        n_samples=len(samples),
        max_maintenance=highest,
        loopless_tightened=tightened,
    )
    returned = sampled.good_reactions
    return result, list(returned) if returned is not None else (good_reactions or [])
