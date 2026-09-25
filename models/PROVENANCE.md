# Model provenance

## ecYeastGEM.yml

The full enzyme-constrained yeast model, in GECKO 4 / RAVEN 3 YAML format.

| | |
|---|---|
| Source | [SysBioChalmers/GECKO](https://github.com/SysBioChalmers/GECKO), `tutorials/full_ecModel/models/ecYeastGEM.yml` |
| Branch | `develop4` |
| Commit | `cbc4ca3311df8ca44c28792f40008b9cbef692bd` |
| Blob | `891ad7ba769cdc7d0fd4b605b0fd625c8cdd6867` |
| Built with | GECKO 4.0.0b1, RAVEN 3.0.0b1 |

1144 enzymes over 4834 enzyme-constrained reactions. kcat sources: 3280 BRENDA,
1076 DLKcat, 251 standard, 217 custom, 6 sensitivity-tuned, 3 isozyme-derived,
1 set explicitly. Subunit stoichiometries come from Complex Portal.

The model ships with `r_4046` (non-growth associated maintenance) pinned at
0.7 mmol ATP/gDW/h and `prot_pool_exchange` bounded at 125 mg/gDW, which is
`p_tot * f * sigma * 1000` for the distributed defaults. Both are replaced
per condition when building the proteome-constrained models.

## yeast-GEM.yml

The conventional GEM the ecModel was built from, from the same commit. Used
as `conv_gem` by the geckopy adapter and as the starting model for random
sampling.

## Verifying

```bash
curl -sL https://raw.githubusercontent.com/SysBioChalmers/GECKO/cbc4ca3311df8ca44c28792f40008b9cbef692bd/tutorials/full_ecModel/models/ecYeastGEM.yml \
  | git hash-object --stdin
# 891ad7ba769cdc7d0fd4b605b0fd625c8cdd6867
```
