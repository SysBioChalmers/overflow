"""Paths, reaction identifiers and experimental conditions.

Values that describe the experiment live in ``data/``; this module only
names them and provides typed access. Nothing here imports a model, so
it stays cheap to import from tests and scripts.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

#: Repository root, i.e. the directory holding ``data/`` and ``models/``.
ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = ROOT / "data"
MODELS_DIR = ROOT / "models"
RESULTS_DIR = ROOT / "results"

EC_MODEL = MODELS_DIR / "ecYeastGEM.yml"
CONV_GEM = MODELS_DIR / "yeast-GEM.yml"

FERMENTATION_DATA = DATA_DIR / "fermentationData.txt"
PROTEOMICS_DATA = DATA_DIR / "abs_proteomics.txt"
RIBOSOME_DATA = DATA_DIR / "ribosome.txt"
ANNOTATION_DATA = DATA_DIR / "selectedAnnotation.txt"

# --- reaction identifiers (yeast-GEM) ---------------------------------
BIO_RXN = "r_4041"          # biomass pseudoreaction; its flux is the growth rate
GROWTH_RXN = "r_2111"       # growth exchange
NGAM_RXN = "r_4046"         # non-growth associated maintenance
PROTEIN_RXN = "r_4047"      # protein pseudoreaction
C_SOURCE = "r_1714"         # D-glucose exchange
CO2_RXN = "r_1672"
O2_RXN = "r_1992"
POOL_RXN = "prot_pool_exchange"

#: Byproduct exchange reactions, keyed by their column name in
#: ``fermentationData.txt``.
BYPRODUCT_RXNS = {
    "Glycerol": "r_1808",
    "Acetate": "r_1634",
    "Ethanol": "r_1761",
    "Formate": "r_1793",
}

#: Oxidative phosphorylation reactions whose complexes the original
#: analysis treated as a unit.
OXPHOS_RXNS = ("r_1021", "r_0439", "r_0438", "r_0226")

#: Prefix of the replicate columns in ``abs_proteomics.txt`` per condition.
#: hGR is the odd one out -- its columns are named ``CN_hGR_*``.
REPLICATE_PREFIX = {
    "CN4": "CN4_",
    "CN22": "CN22_",
    "CN38": "CN38_",
    "CN75": "CN75_",
    "hGR": "CN_hGR_",
}

#: Condition order used throughout the analysis and in every output table.
CONDITION_ORDER = ("CN4", "CN22", "CN38", "CN75", "hGR")


@dataclass(frozen=True)
class Condition:
    """One chemostat condition with its measured rates.

    Rates are mmol/gDW/h and positive as measured, i.e. ``glucose`` and
    ``oxygen`` are uptakes and the rest are secretions. ``byproducts``
    holds the four overflow metabolites; a metabolite that was not
    detected is 0.0, which the model constraints read as "blocked".
    """

    name: str
    p_tot: float          # total protein content [g protein / gDW]
    d_rate: float         # dilution rate == growth rate [1/h]
    glucose: float        # glucose uptake
    co2: float            # CO2 production
    oxygen: float         # oxygen uptake
    byproducts: dict[str, float] = field(default_factory=dict)
    measured: frozenset[str] = frozenset()

    @property
    def replicate_prefix(self) -> str:
        return REPLICATE_PREFIX[self.name]

    def byproduct_bounds(self, flex: float = 1.1) -> dict[str, float]:
        """Upper bounds for the byproduct exchanges.

        Mirrors ``DataConstrains.m``: a measured byproduct may exceed its
        measurement by ``flex``; one that was not detected is blocked.
        """
        return {
            BYPRODUCT_RXNS[name]: flex * value
            for name, value in self.byproducts.items()
        }


def load_conditions(path: Path | str = FERMENTATION_DATA) -> dict[str, Condition]:
    """Read ``fermentationData.txt`` into :class:`Condition` objects."""
    table = pd.read_csv(path, sep="\t")
    conditions: dict[str, Condition] = {}
    for _, row in table.iterrows():
        name = str(row["Condition"])
        byproducts, measured = {}, set()
        for column in BYPRODUCT_RXNS:
            value = row[column]
            if value is None or (isinstance(value, float) and math.isnan(value)):
                byproducts[column] = 0.0
            else:
                byproducts[column] = float(value)
                measured.add(column)
        conditions[name] = Condition(
            name=name,
            p_tot=float(row["Ptot"]),
            d_rate=float(row["D"]),
            glucose=float(row["Glucose"]),
            co2=float(row["CO2"]),
            oxygen=float(row["Oxygen"]),
            byproducts=byproducts,
            measured=frozenset(measured),
        )
    return {name: conditions[name] for name in CONDITION_ORDER}
