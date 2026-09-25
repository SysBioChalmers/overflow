# MATLAB implementation

The GECKO 2 MATLAB pipeline that produced the published results, together
with the models and result tables it generated. It is kept for reference and
is not maintained; `code/prepareEnvironment.m` pins GECKO v2.0.0 and expects
RAVEN 2.4.0.

| | |
|---|---|
| `code/` | The five analysis scripts plus the four patched GECKO functions in `code/customGECKO/` |
| `models/` | GECKO 2 ecModels (`.mat`) and the yeast-GEM releases they were built from |
| `results/` | Tables and figures produced by those scripts |

The Python pipeline in `src/overflow/` reproduces this analysis on a GECKO 4
model. Where a step could be checked against these committed outputs rather
than reimplemented on trust, the test suite does so.

## What does not carry over

GECKO 2 ecModels have no `ec` structure: enzymes appear as `prot_<id>`
metabolites written directly into reaction stoichiometry, with one
`prot_<id>_exchange` reaction per measured enzyme and a `prot_pool` covering
only the unmeasured remainder. geckopy works through `model.ec` and
`usage_prot_<id>` reactions in which every enzyme draws from one shared pool.
The models here therefore cannot be driven from the Python pipeline, and the
four patched functions in `code/customGECKO/` have no direct counterpart —
their behaviour is reimplemented against GECKO 4 primitives instead.
