"""Making measured subunit abundances consistent with complex stoichiometry.

Proteomics recovers the subunits of a membrane-bound complex unevenly,
which leaves a complex whose subunits disagree about how much of it
there is. For the respiratory chain the original analysis replaces the
measurements of a complex's subunits with one abundance per complex,
scaled by each subunit's stoichiometry.

Abundances here are molar (mmol/gDW), because the stoichiometry they are
divided by is molar.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Iterable, Sequence

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - typing only
    from geckopy import EcModel


@dataclass
class ComplexFixResult:
    uniprot_ids: list[str]
    abundances: np.ndarray
    added: list[str] = field(default_factory=list)
    adjusted: list[str] = field(default_factory=list)


def complex_subunits(model: "EcModel", base_reaction: str) -> list[tuple[str, np.ndarray]]:
    """Subunit stoichiometries of every isozyme variant of a reaction.

    Returns one entry per enzyme-constrained variant: its id, and an
    array of ``(enzyme, coefficient)`` pairs.
    """
    matrix = model.ec.rxn_enz_mat
    enzymes = list(model.ec.enzymes)
    variants = []
    for row, reaction_id in enumerate(model.ec.rxns):
        if reaction_id != base_reaction and not reaction_id.startswith(base_reaction + "_"):
            continue
        coefficients = np.asarray(matrix[row, :].todense()).ravel()
        present = np.flatnonzero(coefficients)
        if present.size == 0:
            continue
        variants.append(
            (reaction_id, [(enzymes[i], float(coefficients[i])) for i in present])
        )
    return variants


def fix_complexes(
    model: "EcModel",
    base_reactions: Iterable[str],
    uniprot_ids: Sequence[str],
    abundances: np.ndarray,
) -> ComplexFixResult:
    """Rewrite subunit abundances to one abundance per complex.

    Each subunit is set to ``stoichiometry * mean(measured / stoichiometry)``
    over the subunits that were measured. A subunit of the complex that
    was not measured is added to the data at its stoichiometric share, so
    a complex is never limited by a subunit that simply was not seen.

    A complex with no measured subunit at all is left alone.
    """
    values = {protein: float(value) for protein, value in zip(uniprot_ids, abundances)}
    order = list(uniprot_ids)
    result = ComplexFixResult(uniprot_ids=order, abundances=np.asarray(abundances, float))

    for base in base_reactions:
        for _, subunits in complex_subunits(model, base):
            measured = [
                (protein, values[protein] / stoichiometry)
                for protein, stoichiometry in subunits
                if protein in values and values[protein] > 0
            ]
            if not measured:
                continue
            per_complex = float(np.mean([v for _, v in measured]))
            for protein, stoichiometry in subunits:
                new_value = stoichiometry * per_complex
                if protein in values:
                    if values[protein] != new_value:
                        result.adjusted.append(protein)
                    values[protein] = new_value
                else:
                    order.append(protein)
                    values[protein] = new_value
                    result.added.append(protein)

    result.uniprot_ids = order
    result.abundances = np.array([values[p] for p in order], dtype=float)
    return result
