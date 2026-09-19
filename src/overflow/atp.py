"""Where the ATP comes from and what turns over the redox cofactors.

The sampled flux distributions are summarised into a budget: how much
ATP respiration makes, how much substrate-level phosphorylation makes,
and how much growth and maintenance spend. Expressed per unit of glucose
and per unit of growth, so conditions with different uptake rates can be
compared.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from overflow.config import BIO_RXN, NGAM_RXN

if TYPE_CHECKING:  # pragma: no cover - typing only
    import cobra

#: ATP made by the respiratory chain and by substrate-level
#: phosphorylation in the TCA cycle.
RESPIRATION = ("r_0226", "r_1022")

#: ATP made in glycolysis, and the steps that spend it. yeast-GEM merged
#: glucokinase (r_4235 in 8.3.4) into hexokinase, so r_0534 now carries
#: both and is counted once.
GLYCOLYSIS_PRODUCING = ("r_0892", "r_0962")
GLYCOLYSIS_CONSUMING = ("r_0886", "r_0534")

GLUCOSE_TRANSPORT = "r_1166"

#: Reactions whose flux the published summary reports directly.
TURNOVER = {
    "rPDH": "r_0961",
    "rIDH": "r_0658",
    "rMDHc": "r_0714",
    "rMDHm": "r_0713",
    "rNDE": "r_0770",
}


def _sum(means: Mapping[str, float], reactions: Sequence[str]) -> float:
    return float(sum(float(means.get(r, 0.0)) for r in reactions))


def atp_budget(
    means: Mapping[str, float],
    gam: float,
    growth_rate: float,
    glucose_rate: Optional[float] = None,
    respiration: Sequence[str] = RESPIRATION,
    glycolysis_producing: Sequence[str] = GLYCOLYSIS_PRODUCING,
    glycolysis_consuming: Sequence[str] = GLYCOLYSIS_CONSUMING,
    glucose_transport: str = GLUCOSE_TRANSPORT,
    turnover: Mapping[str, str] = TURNOVER,
) -> pd.Series:
    """The ATP budget implied by one flux distribution.

    ``gam`` is the growth-associated ATP cost written into the biomass
    reaction. Rates are mmol ATP/gDW/h; yields are per mmol glucose and
    per unit of growth.
    """
    if glucose_rate is None:
        glucose_rate = float(means.get(glucose_transport, 0.0))

    respiration_rate = _sum(means, respiration)
    glycolysis = _sum(means, glycolysis_producing) - _sum(means, glycolysis_consuming)
    growth_cost = gam * float(means.get(BIO_RXN, 0.0))
    maintenance = float(means.get(NGAM_RXN, 0.0))
    metabolism = respiration_rate + glycolysis - growth_cost - maintenance

    rows: dict[str, float] = {"rGlu": glucose_rate}
    for name, rate in (
        ("ETC", respiration_rate),
        ("glycolysis", glycolysis),
        ("GAEC", growth_cost),
        ("NGAM", maintenance),
        ("Metabolism", metabolism),
    ):
        rows[f"{name}_rATP"] = rate
        rows[f"{name}_YATP_glu"] = rate / glucose_rate if glucose_rate else np.nan
        rows[f"{name}_YATP_mu"] = rate / growth_rate if growth_rate else np.nan

    spent = growth_cost + maintenance + metabolism
    rows["GAEC+NGAM+Metabolism_rATP"] = spent
    rows["GAEC+NGAM+Metabolism_YATP_glu"] = spent / glucose_rate if glucose_rate else np.nan
    rows["GAEC+NGAM+Metabolism_YATP_mu"] = spent / growth_rate if growth_rate else np.nan

    for name, reaction in turnover.items():
        rows[name] = float(means.get(reaction, 0.0))

    return pd.Series(rows)


def cofactor_turnover(
    model: "cobra.Model",
    means: Mapping[str, float],
    metabolite_name: str,
    label: str,
) -> pd.Series:
    """Production of one cofactor per compartment.

    Only production counts: consumption is its mirror image, and summing
    both would report zero everywhere.
    """
    metabolites = [m for m in model.metabolites if m.name == metabolite_name]
    rows: dict[str, float] = {}
    total = 0.0
    for metabolite in metabolites:
        produced = 0.0
        for reaction in metabolite.reactions:
            flux = float(means.get(reaction.id, 0.0))
            rate = reaction.metabolites[metabolite] * flux
            if rate > 0:
                produced += rate
        if produced == 0:
            continue
        rows[f"{label}[{metabolite.compartment}]"] = produced
        total += produced
    rows[f"{label}[tot]"] = total
    return pd.Series(rows)


def selected_fluxes(
    model: "cobra.Model",
    means: Mapping[str, float],
    gam: float,
    growth_rate: float,
    polymerization: Optional[Mapping[str, float]] = None,
) -> pd.Series:
    """The published summary of one condition's sampled fluxes."""
    parts = [
        atp_budget(means, gam, growth_rate),
        cofactor_turnover(model, means, "NAD", "rNAD"),
        cofactor_turnover(model, means, "NADP(+)", "rNADP"),
    ]
    if polymerization:
        parts.append(pd.Series({f"GAECpol_{k}": v for k, v in polymerization.items()}))
    return pd.concat(parts)
