"""Reading, filtering and unit conversion of the proteomics data.

Abundances are measured in mmol/gDW. geckopy constrains enzymes in
mg/gDW, so everything is converted on the way in, using the molecular
masses in ``data/uniprot.tsv``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import numpy as np
import pandas as pd

from overflow.config import (
    PROTEOMICS_DATA,
    UNIPROT_DATA,
    Condition,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from geckopy.databases import ProtData

ID_COLUMN = "Protein.IDs"


@dataclass(frozen=True)
class FilteredProteomics:
    """Abundances that survived the quality filter, in mmol/gDW.

    ``abundances`` is the upper edge of the replicate confidence
    interval, ``mean + flex_factor * SD``, not the replicate mean.
    """

    uniprot_ids: list[str]
    abundances: np.ndarray
    n_input: int

    @property
    def n_kept(self) -> int:
        return len(self.uniprot_ids)

    @property
    def fraction_removed(self) -> float:
        return 1.0 - self.n_kept / self.n_input if self.n_input else 0.0


def read_proteomics(path: Path | str = PROTEOMICS_DATA) -> pd.DataFrame:
    """Read the abundance table, collapsing repeated protein IDs.

    A protein listed once per gene of a duplicated gene pair appears on
    several rows. Those rows must agree; disagreement means the table
    cannot be reduced to one abundance per protein.
    """
    table = pd.read_csv(path, sep="\t")
    value_columns = [c for c in table.columns if c.endswith("_abs")]
    table[value_columns] = table[value_columns].apply(pd.to_numeric, errors="coerce")

    duplicated = table[ID_COLUMN].duplicated(keep=False)
    if duplicated.any():
        for protein, group in table[duplicated].groupby(ID_COLUMN):
            values = group[value_columns].to_numpy(dtype=float)
            first = values[0]
            same = np.all((values == first) | (np.isnan(values) & np.isnan(first)))
            if not same:
                raise ValueError(
                    f"{protein} appears on {len(group)} rows with differing "
                    "abundances; cannot reduce to one value per protein."
                )
        table = table.drop_duplicates(subset=ID_COLUMN, keep="first")

    return table.reset_index(drop=True)


def replicate_columns(table: pd.DataFrame, condition: Condition | str) -> list[str]:
    """Columns holding the replicates of one condition."""
    from overflow.config import REPLICATE_PREFIX

    prefix = (
        condition.replicate_prefix
        if isinstance(condition, Condition)
        else REPLICATE_PREFIX[condition]
    )
    columns = [c for c in table.columns if c.startswith(prefix)]
    if not columns:
        raise KeyError(f"no replicate columns start with {prefix!r}")
    return columns


def replicate_matrix(
    table: pd.DataFrame, condition: Condition | str
) -> tuple[list[str], np.ndarray]:
    """Protein IDs and their replicate abundances [mmol/gDW]."""
    columns = replicate_columns(table, condition)
    return (
        table[ID_COLUMN].astype(str).tolist(),
        table[columns].to_numpy(dtype=float),
    )


def filter_prot_data(
    uniprot_ids: list[str],
    data: np.ndarray,
    flex_factor: float = 1.96,
    median_filter: bool = True,
    min_value: float = 0.0,
) -> FilteredProteomics:
    """Keep the measurements that are reproducible across replicates.

    A protein is kept when it is detected in at least two thirds of the
    replicates and its relative standard deviation is below 1. The
    abundance reported for it is ``mean + flex_factor * SD``, so a noisy
    protein is represented by the upper edge of its interval rather than
    by its mean.

    A protein with a missing replicate is dropped: its spread cannot be
    evaluated, so it cannot pass the variability test.

    ``median_filter`` divides the standard deviation by the median
    rather than the mean, which is the more stringent of the two.
    """
    n_replicates = data.shape[1]
    kept_ids: list[str] = []
    kept_values: list[float] = []

    for index, protein in enumerate(uniprot_ids):
        if not protein:
            continue
        row = data[index, :]
        if (row > 0).sum() < (2 / 3) * n_replicates:
            continue
        if np.isnan(row).any():
            continue

        sd = float(row.std(ddof=1)) if n_replicates > 1 else 0.0
        centre = float(np.median(row)) if median_filter else float(row.mean())
        if centre == 0:
            continue
        value = float(row.mean()) + flex_factor * sd
        if sd / centre < 1 and value > min_value:
            kept_ids.append(protein)
            kept_values.append(value)

    return FilteredProteomics(
        uniprot_ids=kept_ids,
        abundances=np.asarray(kept_values, dtype=float),
        n_input=len(uniprot_ids),
    )


def molecular_masses(path: Path | str = UNIPROT_DATA) -> dict[str, float]:
    """UniProt molecular masses [Da], which are also mg/mmol."""
    table = pd.read_csv(path, sep="\t")
    return {
        str(entry): float(mass)
        for entry, mass in zip(table["Entry"], table["Mass"])
        if not pd.isna(mass)
    }


def to_mass(
    uniprot_ids: list[str],
    abundances: np.ndarray,
    masses: Optional[dict[str, float]] = None,
) -> np.ndarray:
    """Convert mmol/gDW to mg/gDW.

    One dalton is one mg/mmol, so the conversion is a multiplication by
    the molecular mass. A protein with no known mass raises rather than
    silently contributing zero.
    """
    if masses is None:
        masses = molecular_masses()
    missing = sorted({p for p in uniprot_ids if p not in masses})
    if missing:
        raise KeyError(f"no molecular mass for {len(missing)} proteins: {missing[:5]}")
    return np.asarray(
        [value * masses[protein] for protein, value in zip(uniprot_ids, abundances)],
        dtype=float,
    )


def mean_abundances(
    table: pd.DataFrame, condition: Condition | str
) -> tuple[list[str], np.ndarray]:
    """Replicate means [mmol/gDW], ignoring missing replicates.

    The f-factor and the protein-content rescaling are both derived from
    the unfiltered means, which cover more proteins than the filtered
    set does.
    """
    ids, data = replicate_matrix(table, condition)
    all_missing = np.isnan(data).all(axis=1)
    means = np.full(data.shape[0], np.nan)
    means[~all_missing] = np.nanmean(data[~all_missing], axis=1)
    return ids, means


def f_factor(
    enzymes: list[str],
    uniprot_ids: list[str],
    abundances_mass: np.ndarray,
) -> float:
    """Mass fraction of the measured proteome that the model accounts for."""
    enzyme_set = set(enzymes)
    in_model = np.array([protein in enzyme_set for protein in uniprot_ids])
    total = np.nansum(abundances_mass)
    if total == 0:
        return 0.0
    return float(np.nansum(abundances_mass[in_model]) / total)


def rescaled_p_tot(
    condition: Condition,
    mean_mmol: np.ndarray,
    filtered: FilteredProteomics,
) -> float:
    """Total protein content scaled by what the filter kept.

    The ratio is taken over molar amounts, as in the original analysis.
    Because the filter reports ``mean + 1.96 SD`` rather than the mean,
    the kept fraction can exceed one and the content can scale upwards.
    """
    total = float(np.nansum(mean_mmol))
    kept = float(np.nansum(filtered.abundances))
    if total == 0:
        return condition.p_tot
    return condition.p_tot * (kept / total)


def condition_prot_data(
    condition: Condition,
    enzymes: list[str],
    table: Optional[pd.DataFrame] = None,
    masses: Optional[dict[str, float]] = None,
    flex_factor: float = 1.96,
) -> tuple["ProtData", float, float, FilteredProteomics]:
    """Everything the model build needs from the proteomics of one condition.

    Returns the filtered abundances as a geckopy ``ProtData`` in mg/gDW,
    the f-factor, the rescaled total protein content, and the filter
    result itself for reporting.
    """
    from geckopy.databases import ProtData

    if table is None:
        table = read_proteomics()
    if masses is None:
        masses = molecular_masses()

    ids, means = mean_abundances(table, condition)
    f = f_factor(enzymes, ids, to_mass(ids, means, masses))

    _, data = replicate_matrix(table, condition)
    filtered = filter_prot_data(ids, data, flex_factor=flex_factor)
    p_tot = rescaled_p_tot(condition, means, filtered)

    prot_data = ProtData(
        uniprot_ids=list(filtered.uniprot_ids),
        abundances=to_mass(filtered.uniprot_ids, filtered.abundances, masses),
    )
    return prot_data, f, p_tot, filtered
