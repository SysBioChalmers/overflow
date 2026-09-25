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

Python 3.11 or later. Gurobi (`gurobipy` with a licence) is the default solver for
the analysis and the tests. Another solver can be chosen with `--solver`, or with
the `OVERFLOW_SOLVER` environment variable, which also applies to the tests; an
unavailable solver is an error, not a silent switch.

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
python -m overflow.build CN4 --solver glpk   # instead of the default, Gurobi
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

## Random sampling

```bash
python -m overflow.run_sampling --procs 16
```

Samples the conventional model under the measured rates, with and without
formate in the measured set, and writes the sampled means and standard
deviations, the byproducts the model secretes when it is not told to, and the
ATP and redox budget.

Sampling draws vertices by maximising small random objectives and then minimising
total flux at each, as the published analysis did; means over such draws are vertex
means rather than an average over the interior of the flux space. `--no-min-flux`
skips the second step, which lets the sampler wander into high-flux routes: the
pentose phosphate pathway then carries several times the published flux.

## Comparing against the published results

```bash
python -m overflow.compare
```

Writes `results/COMPARISON.md`: predicted exchange rates, median capacity usage
per system and the sampled ATP and redox budget, this pipeline beside the MATLAB
one, read from `legacy_matlab/results/`.

## Running the whole analysis

The analysis is run with `--fit-rates`, implemented in `overflow.build`
(`relax_to_measured_rates`): it holds CO2, oxygen and the byproducts at their
measurements and raises the measured enzyme abundances by the least that makes
that feasible. Without it those rates are free, and on the earlier tutorial model
the enzyme-constrained solution disposed of surplus carbon through unmeasured
exits instead of respiring it, so the flux distribution did not describe the
measured physiology. `--uptake-flex 1.08` gives CN4 the glucose it needs.

```bash
python -m overflow.build --fit-rates --rate-tolerance 0.08 --uptake-flex 1.08
python -m overflow.build_ribosome
python -m overflow.analyze_usage
python -m overflow.run_sampling --procs 12
python -m overflow.compare
```

Each step reads the models the previous one wrote, so `--models-dir` and
`--results-dir` keep a run self-contained.

## Tests

```bash
pytest                       # everything except the full condition builds
pytest -m "not slow"         # skip the tests that load a genome-scale model
pytest -m integration        # build a condition end to end; minutes
```

The suite runs on the same solver as the analysis, Gurobi by default. Set
`OVERFLOW_SOLVER=glpk` to run it without a licence, as CI does.
