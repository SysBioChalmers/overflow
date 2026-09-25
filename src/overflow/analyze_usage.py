"""Summarising enzyme usage across the condition models.

Run after the models are built: ``python -m overflow.analyze_usage``.
By default it reads the models with the ribosome in them.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import pandas as pd
from geckopy import load_ec_model

from overflow.solver import HELP, use_solver
from overflow.adapter import build_adapter
from overflow.config import (
    BIO_RXN,
    CONDITION_ORDER,
    MODELS_DIR,
    POOL_RXN,
    RESULTS_DIR,
    load_conditions,
)
from overflow.plots import (
    SELECTED_SYSTEMS,
    SUPPLEMENT_SYSTEMS,
    capacity_usage_figure,
)
from overflow.usage import (
    capacity_usage_by_system,
    combine_usage,
    enzyme_usage_table,
    read_annotation,
    usage_by_system,
)


def analyse(
    conditions: Optional[list[str]] = None,
    models_dir: Path = MODELS_DIR,
    suffix: str = "_ribosome",
    verbose: bool = True,
) -> dict[str, pd.DataFrame]:
    """Per-enzyme usage for each condition, from its model's own optimum."""
    all_conditions = load_conditions()
    tables: dict[str, pd.DataFrame] = {}

    for name in conditions or CONDITION_ORDER:
        condition = all_conditions[name]
        path = models_dir / f"ecModel_P_{name}{suffix}.yml"
        model = load_ec_model(str(path), adapter=build_adapter(condition))
        solution = model.optimize()
        if solution.status != "optimal":
            raise RuntimeError(f"{name}: {path.name} has no solution ({solution.status})")
        tables[name] = enzyme_usage_table(model, solution.fluxes)
        if verbose:
            used = tables[name]["absUse"] > 0
            print(
                f"[{name}] growth {solution.fluxes[BIO_RXN]:.5f}, protein "
                f"{solution.fluxes[POOL_RXN]:.1f} mg/gDW, "
                f"{int(used.sum())} of {len(used)} enzymes carrying usage",
                flush=True,
            )
    return tables


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("conditions", nargs="*", default=list(CONDITION_ORDER))
    parser.add_argument("--models-dir", type=Path, default=MODELS_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--suffix", default="_ribosome",
                        help="model file suffix, '' for the models without a ribosome")
    parser.add_argument("--solver", help=HELP)
    args = parser.parse_args(argv)

    use_solver(args.solver)

    tables = analyse(args.conditions, args.models_dir, args.suffix)
    usage = combine_usage(tables)
    capacity = capacity_usage_by_system(usage, read_annotation())

    out = args.results_dir / "enzymeUsage"
    out.mkdir(parents=True, exist_ok=True)
    usage.to_csv(out / "enzymeUsages.tsv", sep="\t", index=False)
    capacity.to_csv(out / "capUsage.tsv", sep="\t", index=False)
    medians = usage_by_system(capacity)
    medians.to_csv(out / "systemMedians.tsv", sep="\t", index=False)

    capacity_usage_figure(capacity, SELECTED_SYSTEMS, out / "selectedSystemUsage.pdf")
    capacity_usage_figure(capacity, SUPPLEMENT_SYSTEMS, out / "supplementSystemUsage.pdf")

    print()
    print("median capacity usage per system (%)")
    print(medians.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
