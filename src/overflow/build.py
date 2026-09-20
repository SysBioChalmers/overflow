"""Building the condition-specific proteome-constrained ecModels.

One model per chemostat condition: biomass rescaled to the measured
protein content, exchange rates constrained to the measurements, enzyme
abundances integrated and reconciled with what the model needs, and
maintenance fitted to the measured gas rates.

Run as ``python -m overflow.build``.
"""
from __future__ import annotations

import argparse
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from geckopy import fill_enz_concs, save_ec_model, set_prot_pool_size

from overflow.adapter import build_adapter, load_model
from overflow.biomass import apply_protein_content
from overflow.config import (
    BIO_RXN,
    C_SOURCE,
    CO2_RXN,
    CONDITION_ORDER,
    MODELS_DIR,
    NGAM_RXN,
    O2_RXN,
    POOL_RXN,
    RESULTS_DIR,
    Condition,
    load_conditions,
)
from overflow.constraints import (
    InfeasibleCondition,
    constrain_byproducts,
    constrain_measured_rates,
    constrain_uptake,
    free_ngam,
    set_chemostat_constraints,
)
from overflow.flexibilize import (
    FlexibilizationResult,
    MinimumUsageResult,
    RateFitResult,
    apply_concentrations,
    flexibilize_proteins,
    minimum_usage_pass,
    relax_to_measured_rates,
    write_back_concentrations,
)
from overflow.ngam import NgamFit, fit_ngam
from overflow.proteomics import (
    condition_prot_data,
    molecular_masses,
    read_proteomics,
)


@dataclass
class BuildResult:
    """What one condition's build produced, for reporting and testing."""

    condition: str
    d_rate: float = 0.0
    model: object = None
    f_factor: float = 0.0
    p_tot_measured: float = 0.0
    p_tot_rescaled: float = 0.0
    pool_size: float = 0.0
    composition: dict = field(default_factory=dict)
    n_measured: int = 0
    n_kept: int = 0
    minimum_usage: Optional[MinimumUsageResult] = None
    flexibilization: Optional[FlexibilizationResult] = None
    rate_fit: Optional[RateFitResult] = None
    ngam: Optional[NgamFit] = None
    fluxes: Optional[pd.Series] = None
    growth: float = 0.0
    glucose: float = 0.0
    seconds: float = 0.0

    @property
    def measured_mass(self) -> float:
        concentrations = np.asarray(self.model.ec.concs, dtype=float)
        return float(np.nansum(concentrations))


def build_condition(
    condition: Condition,
    table: Optional[pd.DataFrame] = None,
    masses: Optional[dict] = None,
    fix_complex_subunits: bool = True,
    ngam_steps: int = 100,
    gam: str = "model",
    scale_protein: bool = False,
    fit_rates: bool = False,
    rate_tolerance: float = 0.05,
    verbose: bool = True,
) -> BuildResult:
    """Build one proteome-constrained ecModel.

    ``scale_protein`` rescales the biomass equation to the measured
    protein content. It is off by default: the model has no carbon to
    spare at the measured glucose uptake, so raising biomass protein
    puts the higher-protein conditions below their dilution rate. The
    measured protein content still sets the protein pool either way.

    ``fit_rates`` additionally holds CO2, oxygen and the byproducts
    within ``rate_tolerance`` of their measurements and raises the
    measured abundances by the least that makes that feasible. Without
    it the gas rates are free, and the model disposes of surplus carbon
    through whatever exit is cheapest in protein rather than respiring
    it.
    """
    started = time.time()
    if table is None:
        table = read_proteomics()
    if masses is None:
        masses = molecular_masses()

    result = BuildResult(
        condition=condition.name,
        d_rate=condition.d_rate,
        p_tot_measured=condition.p_tot,
    )

    model = load_model(build_adapter(condition))
    free_ngam(model)

    # Biomass is rescaled to the protein content as measured, not to the
    # value the proteomics coverage implies.
    result.composition = apply_protein_content(
        model, condition.p_tot, gam=gam, scale_protein=scale_protein
    )

    proteomics = condition_prot_data(
        condition, model, table=table, masses=masses,
        fix_complex_subunits=fix_complex_subunits,
    )
    result.f_factor = proteomics.f_factor
    result.n_kept = proteomics.filtered.n_kept
    result.p_tot_rescaled = proteomics.p_tot

    fill_enz_concs(model, proteomics.prot_data)
    result.n_measured = int((~np.isnan(model.ec.concs)).sum())
    result.pool_size = set_prot_pool_size(
        model, p_tot=result.p_tot_rescaled, f=result.f_factor, sigma=1.0
    )

    constrain_byproducts(model, condition)
    constrain_uptake(model, condition)

    # The model carries no individual enzyme caps yet, so each minimum is
    # what the enzyme needs on its own rather than what the others leave it.
    result.minimum_usage = minimum_usage_pass(model, condition.d_rate)
    apply_concentrations(model, result.minimum_usage.concentrations)

    result.flexibilization = flexibilize_proteins(model, condition.d_rate)
    if fit_rates and result.flexibilization.reached_target:
        constrain_measured_rates(model, condition, tolerance=rate_tolerance)
        result.rate_fit = relax_to_measured_rates(model, condition.d_rate)
        if not result.rate_fit.feasible:
            raise InfeasibleCondition(
                f"{condition.name}: no amount of enzyme relaxation reproduces the "
                f"measured rates within {rate_tolerance:.0%} "
                f"({result.rate_fit.status}). The limit is not the proteome."
            )
    if not result.flexibilization.reached_target:
        raise InfeasibleCondition(
            f"{condition.name}: the model reaches "
            f"{result.flexibilization.growth:.5f} /h after releasing "
            f"{len(result.flexibilization.released)} enzyme caps and growing the "
            f"protein pool {result.flexibilization.pool_increases} times, short of "
            f"the dilution rate {condition.d_rate}. Protein is not what limits it: "
            "check the carbon budget, the biomass composition and the maintenance "
            "cost before flexibilizing further."
        )

    if fit_rates:
        # The measured rates are already the constraint; re-minimising
        # uptake would only pull the model back off them.
        model.objective = {model.reactions.get_by_id(POOL_RXN): -1.0}
        model.objective.direction = "max"
        result.ngam = fit_ngam(model, condition, steps=ngam_steps)
        result.glucose = -float(model.optimize().fluxes[C_SOURCE])
    else:
        set_chemostat_constraints(model, condition)
        result.ngam = fit_ngam(model, condition, steps=ngam_steps)
        result.glucose = set_chemostat_constraints(
            model, condition, glucose=condition.glucose
        )

    solution = model.optimize()
    if solution.status != "optimal":
        raise InfeasibleCondition(
            f"{condition.name}: the finished model has no solution "
            f"({solution.status}) under its own chemostat constraints"
        )
    result.fluxes = solution.fluxes
    result.growth = float(solution.fluxes[BIO_RXN])

    write_back_concentrations(model)
    result.model = model
    result.seconds = time.time() - started

    if verbose:
        report(result)
    return result


def report(result: BuildResult) -> None:
    flex = result.flexibilization
    print(
        f"[{result.condition}] f={result.f_factor:.4f} "
        f"Ptot {result.p_tot_measured:.3f}->{result.p_tot_rescaled:.4f} "
        f"GAM={result.composition.get('GAM', float('nan')):.1f} "
        f"pool={result.pool_size:.1f} mg/gDW\n"
        f"  {result.n_kept} proteins kept, {result.n_measured} are model enzymes, "
        f"{len(result.minimum_usage.raised)} raised to their minimum usage\n"
        f"  flexibilized {len(flex.released)} enzymes, pool grown "
        f"{flex.pool_increases}x ({flex.pool_before:.1f}->{flex.pool_after:.1f})\n"
        + (f"  relaxed {len(result.rate_fit.proteins)} abundances to fit the "
           f"measured rates, adding {result.rate_fit.added_mg:.2f} mg/gDW\n"
           if result.rate_fit else "")
        + f"  NGAM={result.ngam.value:.3f} (error {result.ngam.error:.4f})"
        f"{' AT SCAN EDGE' if result.ngam.at_upper_bound else ''}\n"
        f"  growth={result.growth:.5f} (D={result.d_rate}) "
        f"glucose={result.glucose:.4f}  [{result.seconds:.0f}s]"
    )


def write_outputs(result: BuildResult, models_dir: Path, results_dir: Path) -> None:
    """Write the model and the tables describing how it was built."""
    generation = results_dir / "modelGeneration"
    simulation = results_dir / "modelSimulation"
    generation.mkdir(parents=True, exist_ok=True)
    simulation.mkdir(parents=True, exist_ok=True)

    save_ec_model(result.model, models_dir / f"ecModel_P_{result.condition}.yml")

    result.minimum_usage.table().to_csv(
        generation / f"minimumUsage_{result.condition}.tsv", sep="\t", index=False
    )
    result.flexibilization.table().to_csv(
        generation / f"modifiedEnzymes_{result.condition}.tsv", sep="\t", index=False
    )
    if result.rate_fit is not None:
        result.rate_fit.table().to_csv(
            generation / f"relaxedEnzymes_{result.condition}.tsv", sep="\t", index=False
        )

    model = result.model
    fluxes = result.fluxes
    pd.DataFrame(
        {
            "rxnID": [r.id for r in model.reactions],
            "rxnName": [r.name for r in model.reactions],
            "equation": [r.build_reaction_string(use_metabolite_names=True)
                         for r in model.reactions],
            "flux": [float(fluxes[r.id]) for r in model.reactions],
        }
    ).to_csv(simulation / f"allFluxes_{result.condition}.tsv", sep="\t", index=False)


def summary_table(results: list[BuildResult]) -> pd.DataFrame:
    rows = []
    for result in results:
        model = result.model
        rows.append(
            {
                "condition": result.condition,
                "f": round(result.f_factor, 4),
                "Ptot_measured": result.p_tot_measured,
                "Ptot_rescaled": round(result.p_tot_rescaled, 4),
                "GAM": round(result.composition.get("GAM", float("nan")), 2),
                "pool_mg_gDW": round(result.pool_size, 1),
                "proteins_kept": result.n_kept,
                "model_enzymes_measured": result.n_measured,
                "raised_to_minimum": len(result.minimum_usage.raised),
                "flexibilized": len(result.flexibilization.released),
                "pool_increases": result.flexibilization.pool_increases,
                "rates_relaxed": len(result.rate_fit.proteins) if result.rate_fit else 0,
                "rates_added_mg": round(result.rate_fit.added_mg, 2) if result.rate_fit else 0.0,
                "NGAM": round(result.ngam.value, 3),
                "NGAM_error": round(result.ngam.error, 4),
                "growth": round(result.growth, 5),
                "glucose": round(result.glucose, 4),
                "CO2": round(float(result.fluxes[CO2_RXN]), 4),
                "O2": round(-float(result.fluxes[O2_RXN]), 4),
                "pool_used_mg_gDW": round(float(result.fluxes[POOL_RXN]), 1),
            }
        )
    return pd.DataFrame(rows)


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("conditions", nargs="*", default=list(CONDITION_ORDER))
    parser.add_argument("--models-dir", type=Path, default=MODELS_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--ngam-steps", type=int, default=100)
    parser.add_argument(
        "--solver",
        help="LP solver to use, e.g. gurobi or glpk (default: cobrapy's)",
    )
    parser.add_argument("--no-complex-fix", action="store_true")
    parser.add_argument("--gam", choices=("model", "polymerization"), default="model")
    parser.add_argument(
        "--fit-rates", action="store_true",
        help="hold CO2, oxygen and the byproducts at their measurements and "
             "raise abundances by the least that makes that feasible",
    )
    parser.add_argument("--rate-tolerance", type=float, default=0.05)
    parser.add_argument(
        "--scale-protein",
        action="store_true",
        help=(
            "rescale biomass protein to the measured content. Costs enough "
            "carbon that the higher-protein conditions no longer reach their "
            "dilution rate within the measured uptake."
        ),
    )
    args = parser.parse_args(argv)

    if args.solver:
        import cobra

        cobra.Configuration().solver = args.solver

    conditions = load_conditions()
    table = read_proteomics()
    masses = molecular_masses()

    results = []
    for name in args.conditions or CONDITION_ORDER:
        result = build_condition(
            conditions[name],
            table=table,
            masses=masses,
            fix_complex_subunits=not args.no_complex_fix,
            ngam_steps=args.ngam_steps,
            gam=args.gam,
            scale_protein=args.scale_protein,
            fit_rates=args.fit_rates,
            rate_tolerance=args.rate_tolerance,
        )
        write_outputs(result, args.models_dir, args.results_dir)
        results.append(result)

    summary = summary_table(results)
    args.results_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.results_dir / "modelGeneration" / "summary.tsv",
                   sep="\t", index=False)
    print()
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
