[![Zenodo DOI](https://zenodo.org/badge/288161152.svg)](https://zenodo.org/badge/latestdoi/288161152)

# overflow

Proteome-constrained enzyme-constrained model analysis of overflow metabolism
in *Saccharomyces cerevisiae*, across four carbon/nitrogen ratios and one
high growth rate condition.

The analysis is built on the GECKO 4 `ecYeastGEM` and runs in Python through
[geckopy](https://github.com/SysBioChalmers/geckopy) and
[raven-toolbox](https://github.com/SysBioChalmers/RAVEN). The MATLAB
implementation that produced the published results, which used GECKO 2, is
kept in [`legacy_matlab/`](legacy_matlab/).

## Installation

```bash
pip install -e ".[dev]"
```

Python 3.11 or later. A linear programming solver is required; the analysis
scripts are developed against Gurobi, and the test suite runs on the GLPK
solver that ships with cobrapy.

## Layout

| Path | Contents |
|---|---|
| `data/` | Measured proteomics, fermentation rates, ribosome subunits, annotation |
| `models/` | The GECKO 4 ecModel and the conventional GEM it derives from, see [`models/PROVENANCE.md`](models/PROVENANCE.md) |
| `src/overflow/` | The analysis package |
| `tests/` | Test suite |
| `results/` | Output tables and figures |
| `legacy_matlab/` | The GECKO 2 MATLAB implementation and its results |

## Experimental data

`data/fermentationData.txt` holds, per condition, the total protein content,
the dilution rate and the measured exchange rates. `data/abs_proteomics.txt`
holds absolute protein abundances in mmol/gDW: three biological replicates
per condition, four for hGR.

```python
from overflow import load_conditions

conditions = load_conditions()
conditions["CN4"].d_rate          # 0.1
conditions["CN4"].byproduct_bounds()   # undetected byproducts are blocked
```

## Model

```python
from overflow import build_adapter, load_model

model = load_model(build_adapter(conditions["CN4"]))
```

`build_adapter` reads `model_adapter.toml` and layers the condition's measured
protein content and dilution rate on top.

## Building the condition models

```bash
python -m overflow.build                 # every condition
python -m overflow.build CN4 --solver glpk
```

Each condition gets `models/ecModel_P_<cond>.yml` and, in `results/`, the
enzymes whose measured abundance the model could not run on, the enzymes whose
caps had to be released, and the full flux distribution.

Options change what the models are asked to account for:

- `--fit-rates` holds CO2, oxygen and the byproducts within a tolerance of their
  measurements (`--rate-tolerance`, default 5%) and raises the measured abundances
  by the least that makes that feasible. Without it those rates are free, and the
  model disposes of surplus carbon through whichever exit is cheapest in protein
  rather than respiring it.
- `--uptake-flex` caps glucose uptake at that multiple of the measured rate
  (default 1.05). CN4 needs 1.08: at its dilution rate the model needs 6.5% more
  glucose than was measured, even with protein unlimited.
- `--objective` chooses what the finished model optimises: the smallest protein
  pool (`protein`, the default), or the smallest total flux at the dilution rate
  over every reaction (`flux`) or over the metabolic reactions only
  (`metabolic-flux`). The flux objectives need `--fit-rates`.
- `--scale-protein` rescales biomass to the measured protein content. Off by
  default: the model has no carbon to spare at the measured glucose uptake, and
  CN4 then falls short of its dilution rate.

## Adding the ribosome

```bash
python -m overflow.build_ribosome
```

Writes `models/ecModel_P_<cond>_ribosome.yml`. The protein pseudoreaction is
split so that protein comes out of a `translation` reaction catalysed by the
ribosomal subunits, at 10.5 amino acids per second per ribosome. The core is
the 48 subunits whose average abundance across all conditions reaches
1e-5 mmol/gDW.

## Summarising enzyme usage

```bash
python -m overflow.analyze_usage
```

Writes, per enzyme and condition, how much of it the model uses and what
fraction of what was available that is, plus the capacity usage of the
annotated systems and the two figures over them.

## Tests

```bash
pytest                  # everything
pytest -m "not slow"    # skip the tests that load a genome-scale model
```
