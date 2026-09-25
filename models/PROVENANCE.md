# Model provenance

## ecYeastGEM.yml

The full enzyme-constrained yeast model, in GECKO 4 / RAVEN 3 YAML format.

| | |
|---|---|
| Source | [SysBioChalmers/ecModels](https://github.com/SysBioChalmers/ecModels), `ecYeastGEM/models/ecYeastGEM.yml` |
| Branch | `gecko4` |
| Commit | `880288747ad98220764ce07f8b7a396048e0537d` |
| Blob | `d684eb01ddaf49da9ac9f5891aae372e16b34542` |
| Base model | yeast-GEM 9.1.1 |
| Built with | geckopy 4.0.0b1 |

1144 enzymes over 4850 enzyme-constrained reactions. The kcats were tuned
against measured growth rates and exchange fluxes with CMA-ES, see the
`ecYeastGEM` README in the source repository.

The model ships with `r_4046` (non-growth associated maintenance) pinned at
0.7 mmol ATP/gDW/h and `prot_pool_exchange` bounded at 125 mg/gDW, which is
`p_tot * f * sigma * 1000` for the distributed defaults. Both are replaced
per condition when building the proteome-constrained models.

## yeast-GEM.yml

The conventional GEM the ecModel was built from, yeast-GEM release `v9.1.1`
(`model/yeast-GEM.yml`, sha256
`2ddfaaa73ead41ab243127d1b34fa56db2c86052f6eabcd12bb4eea93947dc17`). Used as
`conv_gem` by the geckopy adapter and as the starting model for random
sampling.

## Verifying

```bash
curl -sL https://raw.githubusercontent.com/SysBioChalmers/ecModels/880288747ad98220764ce07f8b7a396048e0537d/ecYeastGEM/models/ecYeastGEM.yml \
  | git hash-object --stdin
# d684eb01ddaf49da9ac9f5891aae372e16b34542
curl -sL https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/v9.1.1/model/yeast-GEM.yml | sha256sum
```
