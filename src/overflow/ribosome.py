"""Putting the ribosome in the way of protein synthesis.

An ecModel charges protein cost to the enzymes that carry metabolic
flux, but not to the machinery that makes those enzymes. Translation is
the largest single protein expense a growing cell has, so leaving it out
understates the cost of growth.

The protein pseudoreaction is split in two: it now produces a pool of
amino acids, and a new ``translation`` reaction turns that pool into
protein while drawing on the ribosomal subunits. Every gram of protein
the cell makes therefore pays for the ribosomes that made it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Optional, Sequence

import cobra
import numpy as np
import pandas as pd
from geckopy import set_kcat_for_reactions
from geckopy.utilities import NewEnzyme, add_new_rxns_to_ec

from overflow.config import PROTEIN_RXN, RIBOSOME_DATA

if TYPE_CHECKING:  # pragma: no cover - typing only
    from geckopy import EcModel

#: Elongation rate of the yeast ribosome [amino acids/s], doi:10.1042/bj1680409.
AA_PER_SECOND = 10.5

#: Subunits averaging less than this across all conditions [mmol/gDW] are
#: treated as accessory rather than part of the core ribosome.
MIN_MEAN_ABUNDANCE = 1e-5

TRANSLATION_RXN = "translation"
AMINO_ACID_MET = "s_aminoAcids"
AMINO_ACID_NAME = "amino acids for protein"
SUBSYSTEM = "sce03010  Ribosome"


@dataclass
class RibosomeSubunits:
    """The subunits taken to make up the core ribosome."""

    uniprot_ids: list[str]
    genes: dict[str, str]
    masses: dict[str, float]
    mean_abundance: dict[str, float]
    n_candidates: int = 0

    def __len__(self) -> int:
        return len(self.uniprot_ids)


@dataclass
class RibosomeResult:
    """What adding the ribosome changed."""

    subunits: list[str] = field(default_factory=list)
    enzymes_added: list[str] = field(default_factory=list)
    reactions_added: list[str] = field(default_factory=list)
    amino_acid_demand: float = 0.0
    kcat: float = 0.0


def read_ribosome(path: Path | str = RIBOSOME_DATA) -> pd.DataFrame:
    """Read the ribosomal subunit table.

    Masses are written with thousands separators, so ``43,758`` is one
    number rather than two; they are daltons, matching ``ec.mw``.
    """
    table = pd.read_csv(path, sep="\t")
    table = table.rename(
        columns={
            "Entry": "uniprot",
            "Gene names  (primary )": "gene_name",
            "Gene names  (ordered locus )": "gene",
            "Mass": "mass",
            "Sequence": "sequence",
        }
    )
    table["mass"] = pd.to_numeric(
        table["mass"].astype(str).str.replace(",", "", regex=False), errors="coerce"
    )
    return table


def _average_abundance(measured: pd.DataFrame, uniprot: str) -> Optional[float]:
    """Mean over every replicate of every condition, or None if never measured."""
    if uniprot not in measured.index:
        return None
    values = pd.to_numeric(measured.loc[uniprot], errors="coerce").to_numpy(float)
    if np.isnan(values).all():
        return None
    mean = float(np.nanmean(values))
    return mean if np.isfinite(mean) else None


def candidate_means(
    ribosome_table: pd.DataFrame, proteomics_table: pd.DataFrame
) -> pd.Series:
    """Average abundance [mmol/gDW] of every candidate subunit that was measured."""
    value_columns = [c for c in proteomics_table.columns if c.endswith("_abs")]
    measured = proteomics_table.set_index("Protein.IDs")[value_columns]
    means = {}
    for uniprot in ribosome_table["uniprot"].astype(str):
        mean = _average_abundance(measured, uniprot)
        if mean is not None:
            means[uniprot] = mean
    return pd.Series(means, dtype=float)


def core_subunits(
    ribosome_table: pd.DataFrame,
    proteomics_table: pd.DataFrame,
    threshold: float = MIN_MEAN_ABUNDANCE,
) -> RibosomeSubunits:
    """Select the subunits abundant enough to count as the core ribosome.

    Abundance is averaged over every replicate of every condition: a
    subunit belongs to the core or it does not, and that should not be
    decided one condition at a time.
    """
    value_columns = [c for c in proteomics_table.columns if c.endswith("_abs")]
    measured = proteomics_table.set_index("Protein.IDs")[value_columns]

    ids, genes, masses, means = [], {}, {}, {}
    for _, row in ribosome_table.iterrows():
        uniprot = str(row["uniprot"])
        if uniprot not in measured.index:
            continue
        mean = _average_abundance(measured, uniprot)
        if mean is None or mean < threshold:
            continue
        ids.append(uniprot)
        genes[uniprot] = str(row["gene"])
        masses[uniprot] = float(row["mass"])
        means[uniprot] = mean

    return RibosomeSubunits(
        uniprot_ids=ids,
        genes=genes,
        masses=masses,
        mean_abundance=means,
        n_candidates=len(ribosome_table),
    )


def amino_acid_demand(model: cobra.Model, protein_rxn: str = PROTEIN_RXN) -> float:
    """Amino acids consumed per unit of protein [mmol/mmol].

    Read from the protein pseudoreaction's substrates, so it follows any
    rescaling of the biomass composition.
    """
    reaction = model.reactions.get_by_id(protein_rxn)
    return float(-sum(c for c in reaction.metabolites.values() if c < 0))


def translation_kcat(
    model: cobra.Model,
    aa_per_second: float = AA_PER_SECOND,
    protein_rxn: str = PROTEIN_RXN,
) -> float:
    """Turnover of the ribosome in units of protein [1/s].

    A ribosome adds ``aa_per_second`` amino acids per second, and one
    unit of protein takes ``amino_acid_demand`` of them, so in units of
    the protein pseudoreaction it turns over that many times slower.
    """
    return aa_per_second / amino_acid_demand(model, protein_rxn)


def add_ribosome(
    model: "EcModel",
    subunits: RibosomeSubunits,
    aa_per_second: float = AA_PER_SECOND,
    protein_rxn: str = PROTEIN_RXN,
) -> RibosomeResult:
    """Route protein synthesis through a ribosome-catalysed reaction.

    The protein pseudoreaction is left making everything it made before
    except protein itself, which now comes out of ``translation``.
    """
    if AMINO_ACID_MET in {m.id for m in model.metabolites} or any(
        r.id == TRANSLATION_RXN or r.id.startswith(TRANSLATION_RXN + "_")
        for r in model.reactions
    ):
        raise ValueError("the ribosome is already in this model")
    if not subunits.uniprot_ids:
        raise ValueError("no ribosomal subunits to add")

    demand = amino_acid_demand(model, protein_rxn)
    kcat = aa_per_second / demand

    reaction = model.reactions.get_by_id(protein_rxn)
    protein = _protein_metabolite(reaction)
    produced = reaction.metabolites[protein]

    pool = cobra.Metabolite(
        AMINO_ACID_MET, name=AMINO_ACID_NAME, compartment=protein.compartment
    )
    model.add_metabolites([pool])
    reaction.add_metabolites({protein: -produced, pool: produced})

    translation = cobra.Reaction(
        TRANSLATION_RXN, name=TRANSLATION_RXN, lower_bound=0.0, upper_bound=1000.0
    )
    translation.add_metabolites({pool: -produced, protein: produced})
    translation.subsystem = SUBSYSTEM
    translation.gene_reaction_rule = " and ".join(
        subunits.genes[u] for u in subunits.uniprot_ids
    )

    added = add_new_rxns_to_ec(
        model,
        [translation],
        [
            NewEnzyme(enzyme=u, gene=subunits.genes[u], mw=subunits.masses[u])
            for u in subunits.uniprot_ids
        ],
    )
    set_kcat_for_reactions(model, [TRANSLATION_RXN], kcat, apply=True)

    return RibosomeResult(
        subunits=list(subunits.uniprot_ids),
        enzymes_added=list(added.enz_added),
        reactions_added=list(added.rxns_added),
        amino_acid_demand=demand,
        kcat=kcat,
    )


def _protein_metabolite(reaction: cobra.Reaction) -> cobra.Metabolite:
    """The protein the pseudoreaction produces."""
    produced = [m for m, c in reaction.metabolites.items() if c > 0 and m.name == "protein"]
    if len(produced) != 1:
        raise ValueError(
            f"{reaction.id} does not produce exactly one metabolite named "
            f"'protein' (found {len(produced)})"
        )
    return produced[0]


def subunit_coefficient(
    mass: float, demand: float, aa_per_second: float = AA_PER_SECOND
) -> float:
    """Ribosomal protein needed per unit of protein made [mg/mmol]."""
    return mass * demand / (aa_per_second * 3600.0)


@dataclass
class SubunitConstraints:
    """How the measured subunit abundances were applied."""

    applied: dict[str, float] = field(default_factory=dict)
    adjusted: list[str] = field(default_factory=list)
    measured: list[float] = field(default_factory=list)
    required: list[float] = field(default_factory=list)

    def table(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "protein": self.adjusted,
                "measured_mg_gDW": self.measured,
                "adjusted_mg_gDW": self.required,
                "fold_change": [
                    r / m if m else np.inf
                    for m, r in zip(self.measured, self.required)
                ],
            }
        )


def constrain_subunits(
    model: "EcModel",
    abundances: dict[str, float],
    subunits: Optional[Iterable[str]] = None,
    tolerance: float = 1.01,
) -> SubunitConstraints:
    """Cap the ribosomal subunits at their measured abundances [mg/gDW].

    What the model needs is read first, with the subunits unconstrained.
    A subunit measured at less than that is raised to what the model
    needs rather than being allowed to stop translation, on the grounds
    that the cell was demonstrably growing.

    A subunit with no measurement keeps drawing on the protein pool.
    """
    if subunits is None:
        subunits = list(abundances)
    result = SubunitConstraints()

    solution = model.optimize()
    for unit in subunits:
        reaction = model.reactions.get_by_id(f"usage_prot_{unit}")
        measured = abundances.get(unit)
        if measured is None:
            continue
        required = float(solution.fluxes[reaction.id]) if solution.status == "optimal" else 0.0
        if required > measured:
            bound = tolerance * required
            result.adjusted.append(unit)
            result.measured.append(measured)
            result.required.append(bound)
        else:
            bound = measured
        reaction.upper_bound = bound
        result.applied[unit] = bound

    concentrations = np.asarray(model.ec.concs, dtype=float).copy()
    for unit, bound in result.applied.items():
        concentrations[model.ec.enzymes.index(unit)] = bound
    model.ec.concs = concentrations
    return result
