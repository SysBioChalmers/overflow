# Glucose-limited chemostat series (Xia et al. 2022)

The analysis of this repository run on three conditions of the glucose-limited chemostat series of Xia et al. (2022), *Nat Commun* 13:2819, [doi:10.1038/s41467-022-30513-2](https://doi.org/10.1038/s41467-022-30513-2), and compared with the results of this study in `results/`.

| Condition | D [1/h] | Ptot [g/gDW] | Ethanol [mmol/gDW/h] |
|---|---|---|---|
| mu010 | 0.10225 | 0.371 | – |
| mu028 | 0.28425 | 0.391 | – |
| mu033 | 0.33375 | 0.419 | 2.19 |

How the data were converted is in [`data/README.md`](data/README.md).

## Running

`scripts/run_step.py <step> [args]` points `overflow.config` at this folder's `data/`, condition names, `models/` and `results/`, then runs the step's `main`. For `build_ribosome` it supplies the core ribosome of `data/ribosome_core.tsv`. From the repository root, the recipe of this study:

```bash
X=literature/xia2022
O="--models-dir $X/models --results-dir $X/results"
python $X/scripts/run_step.py run_sampling --procs 12 --results-dir $X/results
python $X/scripts/run_step.py respiratory_kcats --results-dir $X/results
python $X/scripts/run_step.py build --fit-rates --rate-tolerance 0.12 --uptake-flex 1.08 --block r_2129 --respiratory-kcats $X/results/modelGeneration/respiratoryKcats.tsv $O
python $X/scripts/run_step.py build_ribosome $O
python $X/scripts/run_step.py analyze_usage $O
```

The kcats of complex III, complex IV and ATP synthase are derived from this series, as they are from this study's data (see the repository README). `results/` holds the tables and figures of this run. The models (`models/`) are not committed.

### Departures from the recipe of this study

- `--rate-tolerance 0.12` instead of 0.08. At 0.08 mu010 has no solution: holding glucose, CO2 and O2 together needs at least 11.6%, even with protein unlimited. As in this study, the protein objective puts predicted glucose at the top of its band: 1.12× measured here, 1.08× there.
- GAM is the model's (55.3), as in this study. The series' respiratory rates imply a lower GAM (about 32 with NGAM 0.64), which was not used.
- The core ribosome is fixed rather than re-selected; see [`data/README.md`](data/README.md#core-ribosome).
- `overflow.compare` is not run: it compares against the MATLAB results of this study's conditions only.
- `--block r_2129` applies to the ecModels. Sampling uses the conventional GEM, where the proton leak stays open, as in this study.

The series attributes all of Ptot to measured proteins, while in this study the measured proteins add up to 72–83% of it. For the respiratory chain this difference in calibration cancels, because each dataset's kcats are derived from its own abundances.

## Comparison with this study

`comparison/` puts both datasets side by side. It needs the ecModels of both runs, so run this study's recipe first.

```bash
for c in CN4 CN22 CN38 CN75 hGR; do for v in built measured; do python $X/comparison/respiratory_limit.py ThisStudy $c $v; done; done
for c in mu010 mu028 mu033;      do for v in built measured; do python $X/comparison/respiratory_limit.py Xia $c $v; done; done
python $X/comparison/figures.py
```

`respiratory_limit.py` finds a condition's respiratory limit: the highest glucose uptake its measured proteome can respire completely at the condition's growth rate. Above it, the model has to secrete overflow carbon. The ecModel is held at the dilution rate with maintenance ATP free, so surplus glucose can be respired until a measured enzyme cap stops it. Caps the pipeline raised so that the model could grow (minimum usage, relaxation to the dilution rate) carry no capacity information and are lifted, for the whole complex. Unlisted exits are closed except glycine export, which levelled models need at some growth rates. Caps the rate fit raised are kept as built (`built`), or put back to the measured abundance (`measured`). The limit is the largest uptake with no more than 0.05 mmol C/gDW/h of forced ethanol, acetate, glycerol and formate. Each run writes its curve and limit to `comparison/respiratory_limit/` (not committed); `figures.py` collects them into Table 5 and figure 5. A condition is predicted to overflow when its measured glucose uptake exceeds its limit; the onset of overflow itself is the measured physiology.

| Output (`comparison/`, figures as PDF) | Contents |
|---|---|
| `figures/fig1_capacity_usage` | Capacity usage of glycolysis, TCA cycle, ETC and ribosome, both datasets (`figS1`: the supplementary systems). |
| `figures/fig2_usage_vs_uptake` | Capacity usage per system (enzyme mass used over the measured capacity of its enzymes) against measured glucose uptake. |
| `figures/fig3_proteome_reserve` | Enzyme mass used in the simulation against the measured capacity, per system. |
| `figures/fig4_atp_budget` | Sampled ATP supply and demand, and the respiratory share of ATP. |
| `figures/fig5_respiratory_limit` | Least overflow the proteome allows against glucose uptake, with the respiratory limits (caps as built) and the measured byproduct carbon. |
| `figures/fig6_apparent_turnover` | Apparent turnover of complex III, complex IV and ATP synthase relative to the highest in each dataset: sampled flux over measured complex, which needs no kcat. |
| `tables/tables.md` | Tables 1–6 (also as `table<n>.tsv`): conditions and fit, capacity usage per system, enzyme mass used and measured capacity, ATP and redox budget, respiratory limit, and the derived respiratory-chain kcats. |
