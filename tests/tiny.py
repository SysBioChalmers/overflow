r"""A hand-built ecModel, small enough to reason about exactly.

    supply -> S -+- R1 (enzyme E1) -+-> P -> biomass
                 \- R2 (enzyme E2) -/

Both routes make the same product, so which one the model uses is
decided purely by enzyme cost and by the caps placed on E1 and E2.
That makes it a precise instrument for the abundance reconciliation and,
in Phase 4, for adding an enzyme to an existing reaction.
"""
from __future__ import annotations

import cobra
import numpy as np
from geckopy import EcModel, ModelAdapter, ModelParameters
from raven_toolbox.io.ec_data import EcData
from scipy import sparse

#: kcat [1/s] and molecular mass [Da] of the two routes. E2 is the
#: cheaper enzyme per unit flux, so an unconstrained model prefers it.
ENZYMES = {"E1": dict(kcat=1.0, mw=36_000.0), "E2": dict(kcat=2.0, mw=36_000.0)}


def protein_coefficient(kcat: float, mw: float, subunits: float = 1.0) -> float:
    """Protein needed per unit flux [mg/gDW per mmol/gDW/h]."""
    return subunits * mw / (kcat * 3600.0)


def tiny_adapter() -> ModelAdapter:
    """Just enough project configuration for the geckopy helpers."""
    return ModelAdapter(
        ModelParameters(
            path=".", conv_gem="none.yml", org_name="test organism",
            sigma=1.0, p_tot=0.5, f=1.0, gr_exp=1.0,
            c_source="EX_S", bio_rxn="BIO", enzyme_comp="cytoplasm",
        )
    )


def tiny_ec_model(pool: float = 100.0, supply: float = 10.0) -> EcModel:
    """Build the model above, with no enzyme measurements."""
    model = cobra.Model("tiny")
    model.compartments = {"c": "cytoplasm"}

    metabolites = {
        name: cobra.Metabolite(name, compartment="c")
        for name in ("S", "P", "biomass", "prot_pool", "prot_E1", "prot_E2")
    }
    model.add_metabolites(list(metabolites.values()))

    supply_rxn = cobra.Reaction("EX_S", name="supply", lower_bound=0.0, upper_bound=supply)
    supply_rxn.add_metabolites({metabolites["S"]: 1.0})

    routes = []
    for enzyme, properties in ENZYMES.items():
        reaction = cobra.Reaction(f"R_{enzyme}", lower_bound=0.0, upper_bound=1000.0)
        reaction.add_metabolites(
            {
                metabolites["S"]: -1.0,
                metabolites["P"]: 1.0,
                metabolites[f"prot_{enzyme}"]: -protein_coefficient(
                    properties["kcat"], properties["mw"]
                ),
            }
        )
        routes.append(reaction)

    biomass = cobra.Reaction("BIO", name="biomass pseudoreaction",
                             lower_bound=0.0, upper_bound=1000.0)
    biomass.add_metabolites({metabolites["P"]: -1.0, metabolites["biomass"]: 1.0})
    sink = cobra.Reaction("EX_BIO", lower_bound=0.0, upper_bound=1000.0)
    sink.add_metabolites({metabolites["biomass"]: -1.0})

    pool_exchange = cobra.Reaction("prot_pool_exchange", lower_bound=0.0, upper_bound=pool)
    pool_exchange.add_metabolites({metabolites["prot_pool"]: 1.0})

    usage = []
    for enzyme in ENZYMES:
        reaction = cobra.Reaction(f"usage_prot_{enzyme}", lower_bound=0.0, upper_bound=1000.0)
        reaction.add_metabolites(
            {metabolites["prot_pool"]: -1.0, metabolites[f"prot_{enzyme}"]: 1.0}
        )
        usage.append(reaction)

    model.add_reactions([supply_rxn, *routes, biomass, sink, pool_exchange, *usage])
    model.objective = "BIO"

    ec_model = EcModel.from_cobra(model, adapter=tiny_adapter())
    names = list(ENZYMES)
    ec_model.ec = EcData(
        gecko_light=False,
        rxns=[f"R_{e}" for e in names],
        kcat=np.array([ENZYMES[e]["kcat"] for e in names], dtype=float),
        source=["test"] * len(names),
        notes=[""] * len(names),
        eccodes=[""] * len(names),
        genes=[f"gene_{e}" for e in names],
        enzymes=names,
        mw=np.array([ENZYMES[e]["mw"] for e in names], dtype=float),
        sequence=["M"] * len(names),
        concs=np.full(len(names), np.nan),
        rxn_enz_mat=sparse.csr_matrix(np.eye(len(names))),
    )
    ec_model.ec.validate()
    return ec_model

