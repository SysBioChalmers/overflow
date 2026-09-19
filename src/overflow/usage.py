"""Summarising how much of each enzyme the model uses.

Two quantities per enzyme: how much of it the flux distribution needs
(absolute usage, mg/gDW) and how much of what was available that is
(capacity usage). An enzyme at full capacity is one the condition is
pressing against; one well below it is being carried.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from geckopy import enzyme_usage

from overflow.config import ANNOTATION_DATA, CONDITION_ORDER

if TYPE_CHECKING:  # pragma: no cover - typing only
    from geckopy import EcModel

ID_COLUMNS = ["protID", "geneID", "protName"]


def enzyme_usage_table(model: "EcModel", fluxes: Mapping[str, float]) -> pd.DataFrame:
    """Per-enzyme usage for one condition."""
    result = enzyme_usage(model, fluxes)
    gene_of = dict(zip(model.ec.enzymes, model.ec.genes))

    def name_of(gene: str) -> str:
        try:
            return model.genes.get_by_id(gene).name or gene
        except KeyError:
            return gene

    genes = [gene_of.get(p, "") for p in result.prot_id]
    return pd.DataFrame(
        {
            "protID": list(result.prot_id),
            "geneID": genes,
            "protName": [name_of(g) for g in genes],
            "capUse": np.asarray(result.cap_usage, dtype=float),
            "absUse": np.asarray(result.abs_usage, dtype=float),
            "UB": np.asarray(result.ub, dtype=float),
        }
    )


def combine_usage(tables: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """One wide table over conditions, in the published column layout."""
    conditions = [c for c in CONDITION_ORDER if c in tables]
    base = tables[conditions[0]][ID_COLUMNS].copy()
    for quantity in ("capUse", "absUse", "UB"):
        for condition in conditions:
            table = tables[condition].set_index("protID")[quantity]
            base[f"{quantity}_{condition}"] = base["protID"].map(table)
    return base


def read_annotation(path: Path | str = ANNOTATION_DATA) -> pd.DataFrame:
    """Read the assignment of proteins to systems."""
    table = pd.read_csv(path, sep="\t")
    return table.rename(
        columns={
            "Entry": "protID",
            "Gene names  (ordered locus )": "geneID",
            "Gene names  (primary )": "protName",
            "system": "system",
        }
    )[["protID", "geneID", "protName", "system"]]


def capacity_usage_by_system(
    usage: pd.DataFrame,
    annotation: Optional[pd.DataFrame] = None,
    conditions: Optional[Sequence[str]] = None,
    decimals: int = 3,
) -> pd.DataFrame:
    """Capacity usage as a percentage, for the annotated systems only.

    Enzymes the model never uses in any condition are dropped: a
    capacity usage of zero everywhere says nothing about how the
    condition allocates protein.
    """
    if annotation is None:
        annotation = read_annotation()
    if conditions is None:
        conditions = [
            c for c in CONDITION_ORDER if f"capUse_{c}" in usage.columns
        ]

    columns = [f"capUse_{c}" for c in conditions]
    table = usage[ID_COLUMNS + columns].copy()
    table = table[table[columns].sum(axis=1) != 0]
    table[columns] = table[columns] * 100

    # A protein listed under more than one system keeps the first, as
    # the published summary does.
    systems = annotation.drop_duplicates(subset="protID", keep="first").set_index(
        "protID"
    )["system"]
    table = table[table["protID"].isin(systems.index)]
    table["GOterm"] = table["protID"].map(systems)

    table = table.rename(columns={f"capUse_{c}": c for c in conditions})
    table[list(conditions)] = table[list(conditions)].round(decimals)
    return table.reset_index(drop=True)


def usage_by_system(capacity: pd.DataFrame, conditions: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Median capacity usage per system and condition."""
    if conditions is None:
        conditions = [c for c in CONDITION_ORDER if c in capacity.columns]
    return (
        capacity.groupby("GOterm")[list(conditions)]
        .median()
        .round(2)
        .reset_index()
    )
