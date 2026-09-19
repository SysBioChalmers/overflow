"""Biomass composition handling for yeast-GEM.

raven-toolbox implements the generic ``sumBioMass``/``scaleBioMass``
arithmetic; this module supplies the yeast-GEM layout it needs, plus the
polymerization costs used when expressing maintenance per component.
"""
from __future__ import annotations

from raven_toolbox.biomass.config import BiomassComponent, BiomassConfig

#: Growth-associated ATP cost of polymerising one gram of each component
#: [mmol ATP / g], from Forster et al. 2003, table S8.
POLYMERIZATION_COST = {
    "protein": 37.7,
    "carbohydrate": 12.8,
    "RNA": 26.0,
    "DNA": 26.0,
}

#: Growth-associated maintenance excluding polymerization [mmol ATP/gDW].
GAM_NO_POLYMERIZATION = 34.0

#: Metabolites in the biomass pseudoreaction that carry the GAM cost.
GAM_COFACTORS = ("ATP", "ADP", "H2O", "H+", "phosphate")

YEAST_BIOMASS = BiomassConfig(
    biomass_rxn="r_4041",
    proton_met="s_0794",
    components=(
        BiomassComponent("protein", "protein pseudoreaction", "mw_minus_2h"),
        BiomassComponent("carbohydrate", "carbohydrate pseudoreaction", "mw"),
        BiomassComponent("RNA", "RNA pseudoreaction", "mw_minus_water"),
        BiomassComponent("DNA", "DNA pseudoreaction", "mw_minus_water"),
        BiomassComponent("lipid backbone", "lipid backbone pseudoreaction", "grams"),
        BiomassComponent("lipid chain", "lipid chain pseudoreaction", "grams"),
        BiomassComponent("ion", "ion pseudoreaction", "mw"),
        BiomassComponent("cofactor", "cofactor pseudoreaction", "mw"),
    ),
)
