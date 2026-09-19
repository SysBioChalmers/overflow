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
) -> tuple[SamplingResult, list[str]]:
    """Sample one condition, returning the summary and the good reactions.

    ``good_reactions`` is the loop-free set the sampler may use as random
    objectives; it depends only on the network, so it is computed once
    and handed back for reuse.
    """
    apply_bounds(model, sampling_bounds(condition, tolerance, include_formate=include_formate))
    highest = constrain_maintenance(model)

    sampled = random_sampling(
        model,
        n_samples,
        method="random_objective",
        n_objectives=2,
        replace_max_bound=True,
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
    )
    returned = sampled.good_reactions
    return result, list(returned) if returned is not None else (good_reactions or [])
