"""Run one overflow CLI step on the Xia et al. (2022) chemostat series.

Points overflow.config at this folder's data/, models/ and results/ and at its
condition names before the step's module is imported, then calls its main().

usage: python literature/xia2022/scripts/run_step.py <module> [args...]

XIA_DATA names another data folder of the same layout, such as the rescaled
one of the coverage check.
"""
import importlib
import os
import sys
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("XIA_DATA", RUN / "data")).resolve()
NAMES = ("mu010", "mu028", "mu033")

import overflow.config as config  # noqa: E402

# the package of this checkout, not another installed copy
assert config.ROOT == RUN.parents[1], f"overflow imported from {config.ROOT}"

config.CONDITION_ORDER = NAMES
config.REPLICATE_PREFIX.clear()
config.REPLICATE_PREFIX.update({name: f"{name}_" for name in NAMES})
config.FERMENTATION_DATA = DATA / "fermentationData.txt"
config.PROTEOMICS_DATA = DATA / "abs_proteomics.txt"
# so a step run without --models-dir/--results-dir cannot write into the repo
config.MODELS_DIR = RUN / "models"
config.RESULTS_DIR = RUN / "results"
config.load_conditions.__defaults__ = (config.FERMENTATION_DATA,)

import overflow.proteomics as proteomics  # noqa: E402

proteomics.read_proteomics.__defaults__ = (config.PROTEOMICS_DATA,)

#: The core ribosome: the 48 subunits of this study, with the paralog this
#: dataset quantifies at the seven positions where it measures the other one
#: (data/ribosome_core.tsv). Re-selecting the core by abundance on this
#: dataset gives 80 subunits, as the paralog groups it quantifies cross the
#: threshold, and leaves mu033 without the protein to grow.
CORE = DATA / "ribosome_core.tsv"


def core_ids() -> list[str]:
    import pandas as pd

    return pd.read_csv(CORE, sep="\t")["uniprot"].tolist()


def fixed_core(ribosome_table, proteomics_table, threshold=None):
    """Stands in for overflow.ribosome.core_subunits."""
    from overflow.ribosome import RibosomeSubunits, _average_abundance

    core = core_ids()
    rows = ribosome_table.set_index(ribosome_table["uniprot"].astype(str))
    value_columns = [c for c in proteomics_table.columns if c.endswith("_abs")]
    measured = proteomics_table.set_index("Protein.IDs")[value_columns]
    means = {u: _average_abundance(measured, u) for u in core if u in measured.index}
    return RibosomeSubunits(
        uniprot_ids=core,
        genes={u: str(rows.loc[u, "gene"]) for u in core},
        masses={u: float(rows.loc[u, "mass"]) for u in core},
        mean_abundance={u: m for u, m in means.items() if m is not None},
        n_candidates=len(ribosome_table),
    )


def mark_core(argv: list[str]) -> None:
    """The step marks the core by the abundance threshold; mark the one used."""
    import argparse

    import pandas as pd

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--results-dir", type=Path, default=config.RESULTS_DIR)
    results_dir = parser.parse_known_args(argv)[0].results_dir
    path = results_dir / "modelGeneration" / "averageRiboSubunitAbundance.tsv"
    table = pd.read_csv(path, sep="\t")
    table["core"] = table["uniprot"].isin(core_ids())
    table.to_csv(path, sep="\t", index=False)


# A forkserver worker re-imports this file; it must not run the step again.
if __name__ == "__main__":
    step = sys.argv[1]
    # the step must bind the patched names, so it cannot have been imported yet
    assert f"overflow.{step}" not in sys.modules
    module = importlib.import_module(f"overflow.{step}")
    assert getattr(module, "CONDITION_ORDER", NAMES) == NAMES
    if step == "build_ribosome":
        module.core_subunits = fixed_core
    status = module.main(sys.argv[2:])
    if step == "build_ribosome":
        mark_core(sys.argv[2:])
    raise SystemExit(status)
