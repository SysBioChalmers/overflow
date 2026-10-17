# Data: glucose-limited chemostat series

Eight conditions from the chemostat series of Xia et al. (2022), *Nat Commun* 13:2819, [doi:10.1038/s41467-022-30513-2](https://doi.org/10.1038/s41467-022-30513-2), in the formats of the pipeline's `data/`.

| Condition | D [1/h] | Proteomics samples | Ptot [g/gDW] |
|---|---|---|---|
| mu003 | 0.027 | S1–S3 (μ = 0.027) | 0.299 |
| mu010 | 0.10225 | S7–S9 (μ = 0.102) | 0.371 |
| mu015 | 0.152 | S10–S12 (μ = 0.152) | 0.419 |
| mu021 | 0.21425 | S13–S15 (μ = 0.214) | 0.320 |
| mu025 | 0.254 | S16–S18 (μ = 0.254) | 0.378 |
| mu028 | 0.28425 | S19–S21 (μ = 0.284) | 0.391 |
| mu033 | 0.33375 | S22–S24 (μ = 0.334) | 0.419 |
| mu038 | 0.3785 | S25–S27 (μ = 0.379) | 0.391 |

The series' μ = 0.044 condition (S4–S6) is not used.

The source spreadsheets of the paper are not committed:

- the extracellular fluxes (sheet `Fig.1B`): D, glucose, CO2, O2 and ethanol rates. Glycerol, acetate and formate are not in it and are `NA`, which the pipeline reads as undetected; ethanol is `NA` where it is 0. With `--fit-rates`, and in sampling, an undetected byproduct may still be secreted at up to 0.01 mmol/gDW/h, as in this study.
- the proteomics, in copies per cell.
- Ptot, the total protein content, is from Supplementary Table 6 of the paper (Supplementary Information, p. 28). Its standard deviation is up to 4.8 percentage points (mu021) and 6.5 (mu038); mu021's 0.320 breaks an otherwise rising series.

`../scripts/prepare_data.py` writes the first two files from the spreadsheets; Ptot is written into its condition table.

- `fermentationData.txt`: the rates and Ptot per condition.
- `abs_proteomics.txt`: copies per cell converted per sample to mmol/gDW by scaling the sample's total protein mass to its condition's Ptot, with the molecular masses in `data/uniprot.tsv` of the repository. Protein groups (`A:B`, mostly ribosomal paralogs) are assigned to their first member with a known mass, isoform suffixes stripped.
- `ribosome_core.tsv`: the core ribosome used (see below).

Scaling to Ptot follows the paper, which derived copies per cell from each protein's share of the total protein content. Each sample's measured protein mass therefore equals Ptot, whereas the measured mass in this study is 72–83% of Ptot. Converting with the paper's 13 pg/cell dry weight instead would give 0.70–1.22× Ptot, from 1.20–1.22× at mu003 to 0.88–0.92× at mu033; the replicates of mu021 (0.80–1.13×) and mu038 (0.70–1.04×) spread the most.

Unlike the data of this study, this series quantifies ribosomal paralog groups (4–7.5% of the measured mass, rising with growth rate). They are not model enzymes, so they lower f.

## Core ribosome

The 48 subunits of this study, except at seven positions where the measured paralog takes the place of this study's member (`this_study_uniprot` gives the one it replaces). At five of them this series has the other paralog 5–295× more abundant, where this study had its own member more abundant: RPL26 (295×), RPL9 (22×), RPS1 (8.6×), RPL21 (8.1×) and RPL8 (5.2×). At the other two this study's member lacks a measurement in some samples, so the pipeline gives it no average: RPS19B is missing in 8 of 24 samples (RPS19A is 5.3× more abundant), and RPL33B is in one sample of 24, where it is 7× more abundant than RPL33A. RPL33A averages 6.9e-6 mmol/gDW, below the core threshold of 1e-5, and is the only swapped subunit the model raises (1.6–4.6×, at every condition but mu003). RPS10B is not in this series; neither RPS10 paralog was detected.

Re-selected by abundance on this series the core has 77 subunits, among them paralog groups, both paralogs of RPL6 and RPL16, and TEF4, so the two datasets would no longer model the same ribosome.
