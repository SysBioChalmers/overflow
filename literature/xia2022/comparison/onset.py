"""Glucose uptake at which each condition's proteome forces carbon overflow.

Each built ecModel is held at its dilution rate. Caps the pipeline raised
because the measurement could not support growth (raised to minimum,
flexibilized) are lifted, since they carry no capacity information. Two
variants for the caps the rate fit raised: kept as built ("built"), or put
back to the measured abundance ("measured"). Maintenance ATP is free, so
surplus glucose can be respired without overflow until a measured enzyme cap
stops it; the onset is the largest uptake with no forced overflow.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import cobra

cobra.Configuration().solver = "gurobi"
cobra.Configuration().processes = 1

import overflow.config as config
from geckopy import load_ec_model
from overflow import build_adapter
from overflow.config import BIO_RXN, C_SOURCE, CO2_RXN, NGAM_RXN, O2_RXN

HERE = Path(__file__).resolve().parent
OUT = HERE / "onset"
REPO = HERE.parents[2]
XIA = HERE.parent
DATASETS = {
    "Kumar": (REPO, ("CN4", "CN22", "CN38", "CN75", "hGR")),  # this repository's data and results
    "Xia": (XIA, ("mu010", "mu028", "mu033")),
    "XiaCoverage": (XIA / "coverage", ("mu010", "mu028", "mu033")),
}
OVERFLOW = {"r_1761": 2, "r_1634": 2, "r_1808": 3, "r_1793": 1}  # ethanol, acetate, glycerol, formate (C)
OPEN = {CO2_RXN, "r_2100", "r_1832", "r_2111", *OVERFLOW}
THRESHOLD = 0.05  # mmol C/gDW/h of forced overflow counted as onset

uniprot = pd.read_csv(REPO / "data" / "uniprot.tsv", sep="\t")
GENE = dict(zip(uniprot["Entry"], uniprot["Gene Names (ordered locus)"]))
annotation = pd.read_csv(REPO / "data" / "selectedAnnotation.txt", sep="\t")
NAME = dict(zip(annotation.iloc[:, 0], annotation.iloc[:, 2]))
SYSTEM = dict(zip(annotation.iloc[:, 0], annotation.iloc[:, 3]))


def label(protein):
    return f"{NAME.get(protein, GENE.get(protein, protein))} ({SYSTEM.get(protein, 'other')})"


def conditions(root, names):
    config.CONDITION_ORDER = names
    return config.load_conditions(root / "data" / "fermentationData.txt")


def prepare(model, cond, generation, variant):
    for r in model.reactions:
        if (len(r.metabolites) == 1 and r.id.startswith("r_") and r.id not in OPEN
                and r.id not in (C_SOURCE, O2_RXN) and r.upper_bound > 0):
            r.upper_bound = 0.001
    for rid in OVERFLOW:
        model.reactions.get_by_id(rid).bounds = (0.0, 100.0)
    model.reactions.get_by_id(CO2_RXN).bounds = (0.0, 1000.0)
    model.reactions.get_by_id(O2_RXN).bounds = (-1000.0, 0.0)
    model.reactions.get_by_id(NGAM_RXN).bounds = (0.0, 1000.0)
    model.reactions.get_by_id(BIO_RXN).bounds = (0.99 * cond.d_rate, cond.d_rate)

    read = lambda kind: pd.read_csv(generation / f"{kind}_{cond.name}.tsv", sep="\t")
    no_information = set(read("minimumUsage")["protein"]) | set(read("modifiedEnzymes")["protein"])
    for p in no_information:
        rid = f"usage_prot_{p}"
        if rid in model.reactions:
            model.reactions.get_by_id(rid).upper_bound = 1000.0
    restored = 0
    if variant == "measured":
        relaxed_path = generation / f"relaxedEnzymes_{cond.name}.tsv"
        if relaxed_path.exists():
            for p, measured in pd.read_csv(relaxed_path, sep="\t")[["protein", "measured_mg_gDW"]].itertuples(index=False):
                rid = f"usage_prot_{p}"
                if p not in no_information and rid in model.reactions:
                    model.reactions.get_by_id(rid).upper_bound = float(measured)
                    restored += 1
    return len(no_information), restored


def overflow_expression(model):
    return sum(w * model.reactions.get_by_id(r).flux_expression for r, w in OVERFLOW.items())


def limiting(model, uptake):
    """Measured caps that bind when overflow is minimised just past the onset."""
    glc = model.reactions.get_by_id(C_SOURCE)
    with model:
        glc.bounds = (-uptake, -uptake)
        model.objective = model.problem.Objective(overflow_expression(model), direction="min")
        sol = model.optimize()
        if sol.status != "optimal":
            return []
        hits = []
        for r in model.reactions:
            if r.id.startswith("usage_prot_") and r.upper_bound < 1000:
                if sol.fluxes[r.id] >= r.upper_bound - 1e-9 and abs(sol.reduced_costs[r.id]) > 1e-6:
                    hits.append((abs(sol.reduced_costs[r.id]), r.id.removeprefix("usage_prot_")))
    return [p for _, p in sorted(hits, reverse=True)[:3]]


def run(dataset, root, cond, variant):
    model = load_ec_model(str(root / "models" / f"ecModel_P_{cond.name}.yml"), adapter=build_adapter(cond))
    lifted, restored = prepare(model, cond, root / "results" / "modelGeneration", variant)
    glc = model.reactions.get_by_id(C_SOURCE)

    window = []
    for direction in ("max", "min"):  # uptake is negative flux
        with model:
            glc.bounds = (-50.0, 0.0)
            model.objective = {glc: 1.0}
            model.objective.direction = direction
            window.append(-model.slim_optimize(error_value=float("nan")))
    umin, umax = window
    row = dict(dataset=dataset, variant=variant, condition=cond.name, D=cond.d_rate,
               glucose_measured=cond.glucose, ethanol_measured=cond.byproducts.get("Ethanol", 0.0),
               caps_lifted=lifted, caps_restored=restored, uptake_min=umin, uptake_max=umax)
    curve = []
    if not np.isfinite(umin):
        row.update(onset=float("nan"), limiting="infeasible")
        return row, curve

    with model:
        glc.bounds = (-50.0, 0.0)
        model.add_cons_vars(model.problem.Constraint(overflow_expression(model), ub=THRESHOLD, name="no_overflow"))
        model.objective = {glc: 1.0}
        model.objective.direction = "min"
        onset = -model.slim_optimize(error_value=float("nan"))
    hits = limiting(model, min(onset * 1.02, umax)) if np.isfinite(onset) and onset < umax - 1e-6 else []
    row.update(onset=onset, limiting="; ".join(label(p) for p in hits) or "none (window end)")

    for u in np.linspace(umin, min(umax, 3.0 * max(onset, cond.glucose)), 30):
        with model:
            glc.bounds = (-u, -u)
            model.objective = model.problem.Objective(overflow_expression(model), direction="min")
            sol = model.optimize()
            if sol.status == "optimal":
                curve.append(dict(dataset=dataset, variant=variant, condition=cond.name, glucose=u,
                                  forced_overflow_C=sol.objective_value,
                                  ethanol=sol.fluxes["r_1761"], CO2=sol.fluxes[CO2_RXN], O2=-sol.fluxes[O2_RXN]))
    return row, curve


if __name__ == "__main__":
    dataset, name, variant = sys.argv[1:4]
    root, names = DATASETS[dataset]
    cond = conditions(root, names)[name]
    row, curve = run(dataset, root, cond, variant)
    OUT.mkdir(exist_ok=True)
    pd.DataFrame([row]).to_csv(OUT / f"onset_{dataset}_{name}_{variant}.tsv", sep="\t", index=False)
    pd.DataFrame(curve).to_csv(OUT / f"curve_{dataset}_{name}_{variant}.tsv", sep="\t", index=False)
    print(pd.DataFrame([row]).T.to_string())
