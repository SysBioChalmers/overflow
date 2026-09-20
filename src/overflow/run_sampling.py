"""Random sampling of every condition.

Two passes per condition: one with formate excluded from the measured
set, which asks the model which byproducts it would make if it were not
told, and one with every measured rate applied.

``python -m overflow.run_sampling``.
"""
from __future__ import annotations

import argparse
import json
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
    use_fork_start_method,
)

#: Mean flux above which a byproduct counts as predicted [mmol/gDW/h].
DETECTION = 0.001


def load_good_reactions(cache: Optional[Path], name: str) -> Optional[list[str]]:
    """The cached loop-free reaction set, if one was kept."""
    if cache is None:
        return None
    path = cache / f"goodReactions_{name}.json"
    if not path.is_file():
        return None
    return json.loads(path.read_text())


def save_good_reactions(cache: Optional[Path], name: str, reactions: list[str]) -> None:
    if cache is None or not reactions:
        return
    cache.mkdir(parents=True, exist_ok=True)
    (cache / f"goodReactions_{name}.json").write_text(json.dumps(sorted(reactions)))


def run(
    conditions: Optional[list[str]] = None,
    n_samples: int = 5000,
    seed: Optional[int] = 0,
    n_proc: Optional[int] = None,
    cache: Optional[Path] = None,
    replace_max_bound: bool = False,
    min_flux: bool = False,
    loopless: bool = True,
    verbose: bool = True,
) -> dict[str, dict]:
    """Sample each condition twice and summarise.

    Which reactions are free of loops depends on the network and not on
    the condition, so that screening is done once and reused; ``cache``
    keeps it across runs, where it is otherwise the dominant cost.
    """
    all_conditions = load_conditions()
    results: dict[str, dict] = {}
    good_free = load_good_reactions(cache, "free")
    good_full = load_good_reactions(cache, "full")
    if verbose and (good_free or good_full):
        print(
            f"reusing cached loop-free sets: "
            f"{len(good_free or [])} / {len(good_full or [])} reactions",
            flush=True,
        )

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
            replace_max_bound=replace_max_bound, min_flux=min_flux,
            loopless=loopless,
        )
        full, good_full = sample_condition(
            model, condition, n_samples=n_samples, include_formate=True,
            seed=seed, good_reactions=good_full, n_proc=n_proc,
            replace_max_bound=replace_max_bound, min_flux=min_flux,
            loopless=loopless,
        )

        summary = selected_fluxes(
            model, full.means, gam=gam, growth_rate=condition.d_rate,
            polymerization=polymerization,
        )
        save_good_reactions(cache, "free", good_free)
        save_good_reactions(cache, "full", good_full)
        results[name] = {
            "free": free, "full": full, "summary": summary,
            "gam": gam, "model": model,
        }
        if verbose:
            print(
                f"[{name}] {full.n_samples} samples, GAM {gam:.2f}, "
                f"maintenance up to {full.max_maintenance:.2f}; "
                f"{full.loopless_tightened} reactions held to their loop-free range; "
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
    parser.add_argument(
        "--cache", type=Path, default=None,
        help="directory to keep the loop-free reaction screening in",
    )
    parser.add_argument(
        "--replace-max-bound", action="store_true",
        help="open the arbitrary 1000 bounds to infinity, as the published "
             "analysis did; the sampler raises if a random objective then "
             "turns out unbounded",
    )
    parser.add_argument(
        "--min-flux", action="store_true",
        help="minimise total flux within each draw, as the published second "
             "pass did; one numerically awkward draw then aborts the run",
    )
    parser.add_argument(
        "--keep-loops", action="store_true",
        help="do not tighten reactions to their loop-free range first",
    )
    parser.add_argument("--solver")
    args = parser.parse_args(argv)

    if args.procs:
        use_fork_start_method()

    if args.solver or args.procs:
        import cobra

        if args.solver:
            cobra.Configuration().solver = args.solver
        if args.procs:
            # The loop-free screening is a flux variability analysis, which
            # reads its process count from here rather than from the
            # sampler's argument. Left unset it runs on one core and
            # dominates the run.
            cobra.Configuration().processes = args.procs

    results = run(
        args.conditions, n_samples=args.samples, seed=args.seed,
        n_proc=args.procs, cache=args.cache,
        replace_max_bound=args.replace_max_bound, min_flux=args.min_flux,
        loopless=not args.keep_loops,
    )

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
