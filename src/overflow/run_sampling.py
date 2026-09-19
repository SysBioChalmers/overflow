"""Random sampling of every condition.

Two passes per condition: one with formate excluded from the measured
set, which asks the model which byproducts it would make if it were not
told, and one with every measured rate applied.

``python -m overflow.run_sampling``.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Optional

import pandas as pd

from overflow.adapter import load_gem
from overflow.atp import selected_fluxes
from overflow.biomass import GAM_NO_POLYMERIZATION, current_gam
from overflow.config import (
    BIO_RXN,
    CONDITION_ORDER,
    RESULTS_DIR,
    load_conditions,
)
from overflow.sampling import (
    polymerization_breakdown,
    prepare_model,
    sample_condition,
)

#: Mean flux above which a byproduct counts as predicted [mmol/gDW/h].
DETECTION = 0.001


def run(
    conditions: Optional[list[str]] = None,
    n_samples: int = 5000,
    seed: Optional[int] = 0,
    n_proc: Optional[int] = None,
    verbose: bool = True,
) -> dict[str, dict]:
    """Sample each condition twice and summarise."""
    all_conditions = load_conditions()
    results: dict[str, dict] = {}
    good_free: Optional[list[str]] = None
    good_full: Optional[list[str]] = None

    for name in conditions or CONDITION_ORDER:
        condition = all_conditions[name]
        started = time.time()

        model = load_gem()
        prepare_model(model, condition, gam_base=GAM_NO_POLYMERIZATION)
        gam = current_gam(model)
        polymerization = polymerization_breakdown(model)

        free, good_free = sample_condition(
            model.copy(), condition, n_samples=n_samples, include_formate=False,
            seed=seed, good_reactions=good_free, n_proc=n_proc,
        )
        full, good_full = sample_condition(
            model, condition, n_samples=n_samples, include_formate=True,
            seed=seed, good_reactions=good_full, n_proc=n_proc,
        )

        summary = selected_fluxes(
            model, full.means, gam=gam, growth_rate=condition.d_rate,
            polymerization=polymerization,
        )
        results[name] = {
            "free": free, "full": full, "summary": summary,
            "gam": gam, "model": model,
        }
        if verbose:
            print(
                f"[{name}] {full.n_samples} samples, GAM {gam:.2f}, "
                f"maintenance up to {full.max_maintenance:.2f}; "
                f"growth {full.means.get(BIO_RXN, float('nan')):.5f} "
                f"[{time.time() - started:.0f}s]",
                flush=True,
            )
    return results


def alternative_exchanges(results: dict[str, dict], model, threshold: float = DETECTION) -> pd.DataFrame:
    """Byproducts the model secretes when it is not told to.

    Taken from the pass with formate left out of the measured set.
    """
    exchanges = [r for r in model.reactions if len(r.metabolites) == 1]
    rows = []
    for reaction in exchanges:
        means = {
            name: float(result["free"].means.get(reaction.id, 0.0))
            for name, result in results.items()
        }
        if max(means.values()) <= threshold:
            continue
        rows.append({"rxnID": reaction.id, "rxnName": reaction.name, **means})
    return pd.DataFrame(rows).sort_values("rxnID")


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("conditions", nargs="*", default=list(CONDITION_ORDER))
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--samples", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--procs", type=int, default=None)
    parser.add_argument("--solver")
    args = parser.parse_args(argv)

    if args.solver:
        import cobra

        cobra.Configuration().solver = args.solver

    results = run(args.conditions, n_samples=args.samples, seed=args.seed, n_proc=args.procs)

    out = args.results_dir / "randomSampling"
    out.mkdir(parents=True, exist_ok=True)

    any_model = next(iter(results.values()))["model"]
    alternative_exchanges(results, any_model).to_csv(
        out / "altExchangeFlux.tsv", sep="\t", index=False
    )

    wide = pd.DataFrame({"rxnID": next(iter(results.values()))["full"].reactions})
    wide["rxnName"] = [any_model.reactions.get_by_id(r).name for r in wide["rxnID"]]
    for name, result in results.items():
        frame = result["full"].to_frame().set_index("rxnID")
        wide[f"{name}_AVERAGE"] = wide["rxnID"].map(frame["mean"])
        wide[f"{name}_STDEV"] = wide["rxnID"].map(frame["sd"])
        wide[f"{name}_STDEV_percent"] = wide["rxnID"].map(frame["sd_percent"])
    wide.to_csv(out / "allFluxes.tsv", sep="\t", index=False)

    summary = pd.DataFrame({name: r["summary"] for name, r in results.items()})
    summary.index.name = "Row"
    summary.to_csv(out / "selectedFluxes.tsv", sep="\t")

    print()
    print(summary.head(24).round(4).to_string())
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
