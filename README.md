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

## Tests

```bash
pytest                  # everything
pytest -m "not slow"    # skip the tests that load a genome-scale model
```
