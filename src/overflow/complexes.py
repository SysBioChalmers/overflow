"""Enzyme complexes: grouping the subunits, and one abundance per complex.

Proteomics recovers the subunits of a complex unevenly: subunits present
1:1 can be measured several-fold apart while their replicates agree. Capped
one by one, a complex could then carry no more flux than its least abundant
measurement allows. Every complex is therefore given one abundance, the
median over its measured subunits of abundance / copies, scaled back to each
subunit by its copies, and the build raises a complex's subunits together.

A complex is a set of enzymes that catalyses a reaction together
(``ec.rxn_enz_mat``). Two enzyme sets are one complex when they share at
least two subunits and at least half of the larger set (isozyme variants
that swap a subunit, the same complex on different reactions). A single
shared protein never merges two sets, so electron carriers and shared
subunits (cytochrome c, thioredoxin, the E3 subunit of the 2-oxo-acid
dehydrogenases) do not chain unrelated complexes together. A subunit in
exactly two complexes belongs to the one whose enzyme sets catalyse more
reactions (COX1-3 of complex IV also appear in an assembly reaction). A
subunit in three or more complexes, in two with equal weight, or that also
catalyses a reaction on its own is *shared*: it keeps its own measurement.

Abundances here are molar (mmol/gDW), because the copies they are divided
by are molar.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional, Sequence

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - typing only
    from geckopy import EcModel

logger = logging.getLogger(__name__)


@dataclass
class ComplexFixResult:
    uniprot_ids: list[str]
    abundances: np.ndarray
    added: list[str] = field(default_factory=list)
    adjusted: list[str] = field(default_factory=list)
    shared: list[str] = field(default_factory=list)


@dataclass
class ComplexFamilies:
    """The complexes of an ecModel.

    ``complexes`` holds one ``{enzyme: copies}`` per complex over every
    subunit in any of its enzyme sets; ``owner`` maps each subunit that
    belongs to a complex to its index; ``shared`` are the subunits that
    belong to none.
    """

    complexes: list[dict[str, float]] = field(default_factory=list)
    owner: dict[str, int] = field(default_factory=dict)
    shared: set[str] = field(default_factory=set)

    def owned(self, index: int) -> dict[str, float]:
        """``{enzyme: copies}`` of the subunits complex ``index`` owns."""
        return {e: c for e, c in self.complexes[index].items() if self.owner.get(e) == index}

    def members(self, enzyme: str) -> list[str]:
        """The subunits of the complex ``enzyme`` belongs to, or ``[enzyme]``."""
        index = self.owner.get(enzyme)
        return [enzyme] if index is None else list(self.owned(index))


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


def complex_families(model: "EcModel", min_overlap: float = 0.5, min_shared: int = 2) -> ComplexFamilies:
    """Group the enzymes of multi-subunit reactions into complexes (see the module text)."""
    matrix = model.ec.rxn_enz_mat.tocsr()
    enzymes = list(model.ec.enzymes)
    sets: dict[frozenset, dict[str, float]] = {}
    reactions: dict[frozenset, int] = {}
    alone: set[str] = set()
    for row in range(matrix.shape[0]):
        start, end = matrix.indptr[row], matrix.indptr[row + 1]
        present = {enzymes[j]: float(v) for j, v in zip(matrix.indices[start:end], matrix.data[start:end]) if v}
        if len(present) == 1:
            alone.update(present)
        elif present:
            key = frozenset(present)
            sets.setdefault(key, present)
            reactions[key] = reactions.get(key, 0) + 1

    keys = list(sets)
    parent = list(range(len(keys)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    by_enzyme: dict[str, list[int]] = {}
    for i, key in enumerate(keys):
        for enzyme in key:
            by_enzyme.setdefault(enzyme, []).append(i)
    for indices in by_enzyme.values():
        for x in range(len(indices)):
            for y in range(x + 1, len(indices)):
                i, j = indices[x], indices[y]
                common = len(keys[i] & keys[j])
                if common >= min_shared and common >= min_overlap * max(len(keys[i]), len(keys[j])):
                    parent[find(i)] = find(j)

    grouped: dict[int, dict[str, float]] = {}
    weight: dict[int, dict[str, int]] = {}
    for i, key in enumerate(keys):
        w = weight.setdefault(find(i), {})
        complex_ = grouped.setdefault(find(i), {})
        for enzyme, copies in sets[key].items():
            w[enzyme] = w.get(enzyme, 0) + reactions[key]
            if enzyme in complex_ and complex_[enzyme] != copies:
                logger.warning("%s has %g copies in one variant and %g in another; the first is used",
                               enzyme, complex_[enzyme], copies)
            complex_.setdefault(enzyme, copies)
    roots = list(grouped)
    complexes = [grouped[r] for r in roots]

    found: dict[str, list[tuple[int, int]]] = {}
    for index, root in enumerate(roots):
        for enzyme, w in weight[root].items():
            found.setdefault(enzyme, []).append((w, index))
    owner: dict[str, int] = {}
    for enzyme, hits in found.items():
        if enzyme in alone or len(hits) > 2:
            continue
        hits.sort(reverse=True)
        if len(hits) == 1 or hits[0][0] > hits[1][0]:
            owner[enzyme] = hits[0][1]
    return ComplexFamilies(complexes=complexes, owner=owner, shared=set(found) - set(owner))


def level_complexes(
    model: "EcModel",
    uniprot_ids: Sequence[str],
    abundances: np.ndarray,
    families: Optional[ComplexFamilies] = None,
) -> ComplexFixResult:
    """Give every complex one abundance: the median per copy over its measured subunits.

    Unmeasured subunits a complex owns are added at their share, so a
    complex is never limited by a subunit that simply was not seen. A
    subunit measured at zero counts as unseen. Shared subunits keep their
    measurement.
    """
    families = families or complex_families(model)
    values = {p: float(a) for p, a in zip(uniprot_ids, abundances)}
    order = list(uniprot_ids)
    result = ComplexFixResult(uniprot_ids=order, abundances=np.asarray(abundances, float),
                              shared=sorted(families.shared & set(values)))
    for index in range(len(families.complexes)):
        members = families.owned(index)
        measured = [values[e] / c for e, c in members.items() if values.get(e, 0) > 0]
        if not measured:
            continue
        per_copy = float(np.median(measured))
        for enzyme, copies in members.items():
            level = per_copy * copies
            if values.get(enzyme, 0) > 0:
                if not np.isclose(values[enzyme], level):
                    result.adjusted.append(enzyme)
            else:
                if enzyme not in values:
                    order.append(enzyme)
                result.added.append(enzyme)
            values[enzyme] = level
    result.uniprot_ids = order
    result.abundances = np.array([values[p] for p in order], dtype=float)
    return result
