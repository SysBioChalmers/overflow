"""Comparing this pipeline's results with the published ones.

The MATLAB implementation's outputs are committed under
``legacy_matlab/results``. Where a quantity means the same thing in both,
it is put side by side here rather than described in prose.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from overflow.config import CONDITION_ORDER, RESULTS_DIR, ROOT, load_conditions
from overflow.ngam import relative_residual
from overflow.usage import usage_by_system

LEGACY = ROOT / "legacy_matlab" / "results"

#: Exchange reactions as the two model formats name them. GECKO 2 split
#: every reversible reaction, so an uptake is a separate reaction there.
EXCHANGES = {
    "glucose": ("r_1714_REV", "r_1714"),
    "CO2": ("r_1672", "r_1672"),
    "oxygen": ("r_1992_REV", "r_1992"),
    "ethanol": ("r_1761", "r_1761"),
}


def _legacy_fluxes(condition: str) -> dict[str, float]:
    fluxes: dict[str, float] = {}
    path = LEGACY / "modelSimulation" / f"allFluxes_{condition}.txt"
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            fluxes[parts[0]] = float(parts[3])
        except ValueError:
            continue
    return fluxes


def _new_fluxes(condition: str, results_dir: Path) -> dict[str, float]:
    table = pd.read_csv(
        results_dir / "modelSimulation" / f"allFluxes_{condition}.tsv", sep="\t"
    )
    return dict(zip(table["rxnID"], table["flux"]))


def exchange_comparison(results_dir: Path = RESULTS_DIR) -> pd.DataFrame:
    """Predicted exchange rates and how far each is from the measurement."""
    conditions = load_conditions()
    rows = []
    for name in CONDITION_ORDER:
        condition = conditions[name]
        legacy, new = _legacy_fluxes(name), _new_fluxes(name, results_dir)
        measured = np.array([condition.glucose, condition.co2, condition.oxygen])
        legacy_rates = np.array([abs(legacy.get(EXCHANGES[k][0], 0.0))
                                 for k in ("glucose", "CO2", "oxygen")])
        new_rates = np.array([abs(new.get(EXCHANGES[k][1], 0.0))
                              for k in ("glucose", "CO2", "oxygen")])
        rows.append(
            {
                "condition": name,
                "glucose_measured": condition.glucose,
                "glucose_gecko2": legacy_rates[0], "glucose_gecko4": new_rates[0],
                "CO2_measured": condition.co2,
                "CO2_gecko2": legacy_rates[1], "CO2_gecko4": new_rates[1],
                "oxygen_measured": condition.oxygen,
                "oxygen_gecko2": legacy_rates[2], "oxygen_gecko4": new_rates[2],
                "ethanol_measured": condition.byproducts["Ethanol"],
                "ethanol_gecko2": legacy.get("r_1761", 0.0),
                "ethanol_gecko4": new.get("r_1761", 0.0),
                "residual_gecko2": relative_residual(legacy_rates, measured),
                "residual_gecko4": relative_residual(new_rates, measured),
            }
        )
    return pd.DataFrame(rows).round(4)


def usage_comparison(results_dir: Path = RESULTS_DIR) -> pd.DataFrame:
    """Median capacity usage per system, side by side."""
    legacy = usage_by_system(
        pd.read_csv(LEGACY / "enzymeUsage" / "capUsage.txt", sep="\t")
    ).set_index("GOterm")
    new = usage_by_system(
        pd.read_csv(results_dir / "enzymeUsage" / "capUsage.tsv", sep="\t")
    ).set_index("GOterm")

    rows = []
    for system in sorted(set(legacy.index) | set(new.index)):
        row: dict[str, object] = {"system": system}
        for condition in CONDITION_ORDER:
            row[f"{condition}_gecko2"] = (
                legacy.loc[system, condition] if system in legacy.index else np.nan
            )
            row[f"{condition}_gecko4"] = (
                new.loc[system, condition] if system in new.index else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


#: Rows of the sampling summary worth putting side by side.
BUDGET_ROWS = (
    "rGlu", "ETC_rATP", "glycolysis_rATP", "GAEC_rATP", "NGAM_rATP",
    "GAEC+NGAM+Metabolism_rATP", "rPDH", "rIDH", "rMDHc", "rMDHm", "rNDE",
)


def budget_comparison(results_dir: Path = RESULTS_DIR) -> pd.DataFrame:
    """The sampled ATP and redox budget, side by side."""
    legacy = pd.read_csv(
        LEGACY / "randomSampling" / "selectedFluxes.txt", sep="\t"
    ).set_index("Row")
    new = pd.read_csv(
        results_dir / "randomSampling" / "selectedFluxes.tsv", sep="\t"
    ).set_index("Row")

    rows = []
    for name in BUDGET_ROWS:
        if name not in legacy.index or name not in new.index:
            continue
        row: dict[str, object] = {"row": name}
        for condition in CONDITION_ORDER:
            row[f"{condition}_gecko2"] = float(legacy.loc[name, condition])
            row[f"{condition}_gecko4"] = float(new.loc[name, condition])
        rows.append(row)
    return pd.DataFrame(rows)


def report(results_dir: Path = RESULTS_DIR) -> str:
    """A written comparison of the two implementations."""
    exchanges = exchange_comparison(results_dir)
    lines = [
        "# This pipeline against the published one",
        "",
        "Generated by `python -m overflow.compare`. GECKO 2 values are read from",
        "`legacy_matlab/results/`, which is what the MATLAB implementation wrote.",
        "",
        "## Exchange rates",
        "",
        "Measured / GECKO 2 / GECKO 4, mmol/gDW/h. The residual is the relative",
        "deviation over glucose, CO2 and oxygen together.",
        "",
        "| condition | glucose | CO2 | oxygen | ethanol | residual G2 | residual G4 |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in exchanges.itertuples():
        lines.append(
            f"| {row.condition} "
            f"| {row.glucose_measured:.3f} / {row.glucose_gecko2:.3f} / {row.glucose_gecko4:.3f} "
            f"| {row.CO2_measured:.3f} / {row.CO2_gecko2:.3f} / {row.CO2_gecko4:.3f} "
            f"| {row.oxygen_measured:.3f} / {row.oxygen_gecko2:.3f} / {row.oxygen_gecko4:.3f} "
            f"| {row.ethanol_measured:.3f} / {row.ethanol_gecko2:.3f} / {row.ethanol_gecko4:.3f} "
            f"| {row.residual_gecko2:.3f} | {row.residual_gecko4:.3f} |"
        )
    lines += [
        "",
        f"Mean residual: GECKO 2 {exchanges.residual_gecko2.mean():.3f}, "
        f"GECKO 4 {exchanges.residual_gecko4.mean():.3f}.",
        "",
    ]

    try:
        usage = usage_comparison(results_dir)
    except FileNotFoundError:
        return "\n".join(lines)

    lines += [
        "## Median capacity usage per system (%)",
        "",
        "| system | " + " | ".join(CONDITION_ORDER) + " |",
        "|---" * (len(CONDITION_ORDER) + 1) + "|",
    ]
    for row in usage.itertuples(index=False):
        values = []
        for condition in CONDITION_ORDER:
            old = getattr(row, f"{condition}_gecko2")
            new = getattr(row, f"{condition}_gecko4")
            values.append(f"{old:.1f} / {new:.1f}")
        lines.append(f"| {row.system} | " + " | ".join(values) + " |")
    lines.append("")

    try:
        budget = budget_comparison(results_dir)
    except FileNotFoundError:
        return "\n".join(lines)

    lines += [
        "## Sampled ATP and redox budget",
        "",
        "mmol/gDW/h, GECKO 2 / GECKO 4.",
        "",
        "| | " + " | ".join(CONDITION_ORDER) + " |",
        "|---" * (len(CONDITION_ORDER) + 1) + "|",
    ]
    for row in budget.itertuples(index=False):
        values = [
            f"{getattr(row, f'{c}_gecko2'):.3g} / {getattr(row, f'{c}_gecko4'):.3g}"
            for c in CONDITION_ORDER
        ]
        lines.append(f"| {row.row} | " + " | ".join(values) + " |")
    lines.append("")
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    text = report(args.results_dir)
    destination = args.output or (args.results_dir / "COMPARISON.md")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
