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

#: Metabolites in the biomass pseudoreaction whose coefficient *is* the
#: GAM. The proton is not among them: yeast-GEM's biomass releases three
#: fewer protons than it hydrolyses ATP, because of the NADPH it also
#: consumes, and that offset has to survive a change of GAM.
GAM_COFACTORS = ("ATP", "ADP", "H2O", "phosphate")

#: The proton, which moves with the GAM but keeps its own offset.
GAM_PROTON = "H+"

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


def polymerization_cost(model) -> float:
    """ATP spent polymerising one gram of biomass [mmol ATP/gDW].

    Each macromolecular component costs its own rate per gram, so the
    cost depends on the composition and has to be recomputed whenever
    the composition is rescaled.
    """
    from raven_toolbox.biomass.scale import sum_biomass

    fractions = sum_biomass(model, YEAST_BIOMASS)
    return sum(
        cost * fractions.get(component, 0.0)
        for component, cost in POLYMERIZATION_COST.items()
    )


def current_gam(model) -> float:
    """The growth-associated ATP cost written into the biomass reaction."""
    reaction = model.reactions.get_by_id(YEAST_BIOMASS.biomass_rxn)
    atp = next(m for m in reaction.metabolites if m.name == "ATP")
    return -reaction.metabolites[atp]


def set_growth_maintenance(model, value: float) -> float:
    """Set the growth-associated ATP cost of the biomass reaction.

    ATP, ADP, water and phosphate take the new value directly. The
    proton is shifted by the same amount instead of being set to it, so
    that the reaction keeps whatever proton offset it was built with;
    setting it outright would quietly add protons to every gram of
    biomass and cost the model growth it should have had.

    Returns the previous value.
    """
    reaction = model.reactions.get_by_id(YEAST_BIOMASS.biomass_rxn)
    previous = current_gam(model)
    delta = value - previous
    if delta == 0:
        return previous

    updates = {}
    for metabolite, coefficient in reaction.metabolites.items():
        if metabolite.name in GAM_COFACTORS:
            updates[metabolite] = (1 if coefficient > 0 else -1) * value - coefficient
        elif metabolite.name == GAM_PROTON and coefficient > 0:
            updates[metabolite] = delta
    reaction.add_metabolites(updates)
    return previous


def apply_protein_content(
    model,
    p_tot: float,
    gam: str = "model",
    gam_base: float = GAM_NO_POLYMERIZATION,
    scale_protein: bool = True,
) -> dict[str, float]:
    """Set the biomass composition and maintenance for one condition.

    With ``scale_protein``, protein is scaled to the measured ``p_tot``
    and carbohydrate takes up the difference, keeping biomass at one
    gram per gram.

    ``gam`` selects the growth-associated maintenance:

    ``"model"``
        keep the value the distributed model carries, which was fitted
        for this model's own stoichiometry.
    ``"polymerization"``
        ``gam_base`` plus the polymerization cost of the resulting
        composition, which is how the MATLAB analysis set it.

    Returns the resulting composition, with the GAM under ``"GAM"``.
    """
    from raven_toolbox.biomass.scale import scale_biomass, sum_biomass

    if gam not in ("model", "polymerization"):
        raise ValueError(f"unknown gam policy {gam!r}")

    kept = current_gam(model)
    if scale_protein:
        scale_biomass(model, YEAST_BIOMASS, "protein", p_tot, balance_out="carbohydrate")

    value = kept if gam == "model" else gam_base + polymerization_cost(model)
    set_growth_maintenance(model, value)
    composition = dict(sum_biomass(model, YEAST_BIOMASS))
    composition["GAM"] = value
    return composition
