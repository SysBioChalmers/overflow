"""Adding the ribosome to each condition model.

Run after ``overflow.build``, on the models it wrote:
``python -m overflow.build_ribosome``.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from geckopy import load_ec_model, save_ec_model

from overflow.adapter import build_adapter
from overflow.config import (
    BIO_RXN,
    CONDITION_ORDER,
    MODELS_DIR,
    POOL_RXN,
    RESULTS_DIR,
    Condition,
    load_conditions,
)
from overflow.proteomics import (
    filter_prot_data,
    molecular_masses,
    read_proteomics,
    replicate_matrix,
    to_mass,
)
from overflow.ribosome import (
    RibosomeResult,
    RibosomeSubunits,
    SubunitConstraints,
    add_ribosome,
    constrain_subunits,
    core_subunits,
    read_ribosome,
)


@dataclass
class RibosomeBuild:
    condition: str
    model: object = None
    ribosome: Optional[RibosomeResult] = None
    constraints: Optional[SubunitConstraints] = None
    growth_before: float = 0.0
    growth_after: float = 0.0
    pool_before: float = 0.0
    pool_after: float = 0.0
    ribosome_mass: float = 0.0


def subunit_abundances(
    condition: Condition,
    subunits: RibosomeSubunits,
    table: Optional[pd.DataFrame] = None,
    masses: Optional[dict[str, float]] = None,
) -> dict[str, float]:
    """Measured subunit abundances for one condition [mg/gDW].

    Filtered on the same terms as every other protein, so a subunit that
    was measured unreliably in this condition carries no cap here even
    though its average across conditions put it in the core.
    """
    if table is None:
        table = read_proteomics()
    if masses is None:
        masses = molecular_masses()

    ids, replicates = replicate_matrix(table, condition)
    filtered = filter_prot_data(ids, replicates)
    wanted = set(subunits.uniprot_ids)
    keep = [(p, v) for p, v in zip(filtered.uniprot_ids, filtered.abundances) if p in wanted]
    if not keep:
        return {}
    proteins = [p for p, _ in keep]
    values = np.array([v for _, v in keep], dtype=float)
    return dict(zip(proteins, to_mass(proteins, values, masses)))


def add_ribosome_to_condition(
    condition: Condition,
    models_dir: Path = MODELS_DIR,
    subunits: Optional[RibosomeSubunits] = None,
    table: Optional[pd.DataFrame] = None,
    masses: Optional[dict[str, float]] = None,
    verbose: bool = True,
) -> RibosomeBuild:
    """Load a condition model, put the ribosome in it and constrain it."""
    if subunits is None:
        subunits = core_subunits(read_ribosome(), table if table is not None else read_proteomics())

    path = models_dir / f"ecModel_P_{condition.name}.yml"
    model = load_ec_model(str(path), adapter=build_adapter(condition))
    result = RibosomeBuild(condition=condition.name, model=model)

    solution = model.optimize()
    result.growth_before = float(solution.fluxes[BIO_RXN])
    result.pool_before = float(solution.fluxes[POOL_RXN])

    result.ribosome = add_ribosome(model, subunits)
    result.constraints = constrain_subunits(
        model, subunit_abundances(condition, subunits, table, masses), subunits.uniprot_ids
    )

    solution = model.optimize()
    result.growth_after = float(solution.fluxes[BIO_RXN])
    result.pool_after = float(solution.fluxes[POOL_RXN])
    result.ribosome_mass = sum(
        abs(float(solution.fluxes[f"usage_prot_{u}"])) for u in subunits.uniprot_ids
    )

    if verbose:
        print(
            f"[{condition.name}] {len(subunits)} subunits, kcat "
            f"{result.ribosome.kcat:.4f} /s ({result.ribosome.amino_acid_demand:.4f} "
            f"amino acids per protein)\n"
            f"  {len(result.constraints.adjusted)} subunits raised to what "
            f"translation needs\n"
            f"  growth {result.growth_before:.5f} -> {result.growth_after:.5f}; "
            f"protein pool {result.pool_before:.1f} -> {result.pool_after:.1f} mg/gDW, "
            f"of which ribosome {result.ribosome_mass:.1f}",
            flush=True,
        )
    return result


def summary_table(results: list[RibosomeBuild]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "condition": r.condition,
                "subunits": len(r.ribosome.subunits),
                "translation_kcat": round(r.ribosome.kcat, 4),
                "aa_per_protein": round(r.ribosome.amino_acid_demand, 4),
                "subunits_raised": len(r.constraints.adjusted),
                "growth_before": round(r.growth_before, 5),
                "growth_after": round(r.growth_after, 5),
                "pool_before": round(r.pool_before, 1),
                "pool_after": round(r.pool_after, 1),
                "ribosome_mg_gDW": round(r.ribosome_mass, 1),
                "ribosome_fraction": round(r.ribosome_mass / r.pool_after, 4)
                if r.pool_after
                else float("nan"),
            }
            for r in results
        ]
    )


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("conditions", nargs="*", default=list(CONDITION_ORDER))
    parser.add_argument("--models-dir", type=Path, default=MODELS_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument(
        "--solver",
        help="LP solver to use, e.g. gurobi or glpk (default: cobrapy's)",
    )
    args = parser.parse_args(argv)

    if args.solver:
        import cobra

        cobra.Configuration().solver = args.solver

    conditions = load_conditions()
    table = read_proteomics()
    masses = molecular_masses()
    subunits = core_subunits(read_ribosome(), table)
    print(f"core ribosome: {len(subunits)} of {subunits.n_candidates} subunits", flush=True)

    generation = args.results_dir / "modelGeneration"
    generation.mkdir(parents=True, exist_ok=True)

    results = []
    for name in args.conditions or CONDITION_ORDER:
        result = add_ribosome_to_condition(
            conditions[name], models_dir=args.models_dir,
            subunits=subunits, table=table, masses=masses,
        )
        save_ec_model(result.model, args.models_dir / f"ecModel_P_{name}_ribosome.yml")
        result.constraints.table().to_csv(
            generation / f"modifiedRibosomeSubunits_{name}.tsv", sep="\t", index=False
        )
        results.append(result)

    summary = summary_table(results)
    summary.to_csv(generation / "ribosomeSummary.tsv", sep="\t", index=False)
    print()
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
