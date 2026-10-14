"""Respiratory-chain kcats derived from the data.

With the ecModel's kcats the measured respiratory chain cannot carry the
measured respiration in several conditions, so the build raises its
abundances to whatever the fluxes need and its capacity usage reads 100%
by construction. Instead, the kcats of complex III, complex IV and ATP
synthase are derived from the conditions themselves.

In each condition the apparent turnover of a complex is its sampled flux
(``randomSampling/allFluxes.tsv``, the conventional model constrained by
the measured rates, which uses no kcats) over the complex abundance the
build uses: the filtered mean + 1.96 SD of its measured subunits,
levelled to one abundance per complex. The highest apparent turnover
across the conditions becomes the complex's kcat, for every isozyme
variant of the forward reaction. The condition that sets it runs the
complex at full capacity, and the derived kcat is a lower bound.

Writes ``modelGeneration/respiratoryTurnover.tsv`` (every condition) and
``modelGeneration/respiratoryKcats.tsv`` (what ``overflow.build
--respiratory-kcats`` applies).
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Mapping, Optional

import numpy as np
import pandas as pd

from overflow.solver import HELP, use_solver
from overflow.adapter import load_model
from overflow.complexes import complex_subunits
from overflow.config import CONDITION_ORDER, RESULTS_DIR, load_conditions
from overflow.proteomics import condition_prot_data, molecular_masses, read_proteomics

#: The complexes whose kcats are derived, by base reaction.
COMPLEXES = {"complex III": "r_0439", "complex IV": "r_0438", "ATP synthase": "r_0226"}


def variants(model, base: str) -> list[str]:
    """The forward isozyme variants of a reaction in ``model.ec.rxns``."""
    return [r for r in model.ec.rxns if r == base or r.startswith(base + "_EXP_")]


def complex_abundance(model, base: str, molar: Mapping[str, float]) -> float:
    """Abundance of a complex [mmol/gDW]: mean over its measured subunits of abundance / copies."""
    _, subunits = complex_subunits(model, base)[0]
    values = [molar[protein] / copies for protein, copies in subunits if molar.get(protein, 0) > 0]
    return float(np.mean(values)) if values else float("nan")


def apparent_turnovers(model, conditions, sampled: pd.DataFrame, table=None, masses=None) -> pd.DataFrame:
    """Apparent turnover [1/s] of each complex in each condition."""
    if table is None:
        table = read_proteomics()
    if masses is None:
        masses = molecular_masses()
    kcat = dict(zip(model.ec.rxns, np.asarray(model.ec.kcat, dtype=float)))
    rows = []
    for name, condition in conditions.items():
        levelled = condition_prot_data(condition, model, table=table, masses=masses).complex_fix
        molar = dict(zip(levelled.uniprot_ids, levelled.abundances))
        for complex_name, base in COMPLEXES.items():
            flux = abs(float(sampled.loc[base, f"{name}_AVERAGE"]))
            abundance = complex_abundance(model, base, molar)
            rows.append(dict(
                condition=name, complex=complex_name, reaction=base,
                flux_mmol_gDW_h=flux, abundance_nmol_gDW=abundance * 1e6,
                apparent_turnover_per_s=flux / (abundance * 3600),
                ecmodel_kcat_per_s=kcat[variants(model, base)[0]],
            ))
    out = pd.DataFrame(rows)
    out["relative_to_highest"] = out["apparent_turnover_per_s"] / out.groupby("complex")["apparent_turnover_per_s"].transform("max")
    return out


def derived_kcats(turnovers: pd.DataFrame) -> pd.DataFrame:
    """The highest apparent turnover of each complex, and the condition it comes from."""
    best = turnovers.loc[turnovers.groupby("complex")["apparent_turnover_per_s"].idxmax()]
    return best[["complex", "reaction", "condition", "apparent_turnover_per_s", "ecmodel_kcat_per_s"]].rename(
        columns={"apparent_turnover_per_s": "kcat_per_s"}).reset_index(drop=True)


def apply_respiratory_kcats(model, kcats: pd.DataFrame) -> list[str]:
    """Set each complex's kcat on all its forward variants and rewrite their enzyme coefficients."""
    from geckopy.ec_model.pipeline.apply_kcat import apply_kcat_constraints

    index = {r: i for i, r in enumerate(model.ec.rxns)}
    updated = []
    for base, value in zip(kcats["reaction"], kcats["kcat_per_s"]):
        for reaction in variants(model, base):
            model.ec.kcat[index[reaction]] = float(value)
            updated.append(reaction)
    apply_kcat_constraints(model, update_rxns=updated)
    return updated


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--solver", help=HELP)
    args = parser.parse_args(argv)
    use_solver(args.solver)

    sampled = pd.read_csv(args.results_dir / "randomSampling" / "allFluxes.tsv", sep="\t").set_index("rxnID")
    conditions = load_conditions()
    turnovers = apparent_turnovers(load_model(), {c: conditions[c] for c in CONDITION_ORDER}, sampled)
    kcats = derived_kcats(turnovers)
    generation = args.results_dir / "modelGeneration"
    generation.mkdir(parents=True, exist_ok=True)
    turnovers.round(6).to_csv(generation / "respiratoryTurnover.tsv", sep="\t", index=False)
    kcats.round(6).to_csv(generation / "respiratoryKcats.tsv", sep="\t", index=False)
    print(turnovers.round(3).to_string(index=False))
    print()
    print(kcats.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
