# Data: glucose-limited chemostat series

Three conditions from the chemostat series of Xia et al. (2022), *Nat Commun* 13:2819, [doi:10.1038/s41467-022-30513-2](https://doi.org/10.1038/s41467-022-30513-2), in the formats of the pipeline's `data/`.

| Condition | D [1/h] | Proteomics samples | Ptot [g/gDW] |
|---|---|---|---|
| mu010 | 0.10225 | S7–S9 (μ = 0.102) | 0.371 |
| mu028 | 0.28425 | S19–S21 (μ = 0.284) | 0.391 |
| mu033 | 0.33375 | S22–S24 (μ = 0.334) | 0.419 |

The source spreadsheets of the paper are not committed:

- the extracellular fluxes (sheet `Fig.1B`): D, glucose, CO2, O2 and ethanol rates. Glycerol, acetate and formate are not in it and are `NA`, which the pipeline reads as undetected; ethanol is `NA` where it is 0. With `--fit-rates`, and in sampling, an undetected byproduct may still be secreted at up to 0.01 mmol/gDW/h, as in the main run.
- the proteomics, in copies per cell.
- Ptot, the total protein content, is from Supplementary Table 6 of the paper (Supplementary Information, p. 28).

`../scripts/prepare_data.py` writes the first two files from the spreadsheets.

- `fermentationData.txt`: the rates and Ptot per condition.
- `abs_proteomics.txt`: copies per cell converted per sample to mmol/gDW by scaling the sample's total protein mass to its condition's Ptot, with the molecular masses in `data/uniprot.tsv` of the repository. Protein groups (`A:B`, mostly ribosomal paralogs) are assigned to their first member with a known mass, isoform suffixes stripped.
- `ribosome_core.tsv`: the core ribosome used (see below).

Scaling to Ptot follows the paper, which derived copies per cell from each protein's share of the total protein content. Each sample's measured protein mass therefore equals Ptot, whereas the measured mass in the main data is 70–84% of Ptot. Converting with the paper's 13 pg/cell dry weight instead would give 1.07–1.11× Ptot at mu010, 0.94–0.99× at mu028 and 0.88–0.92× at mu033.

Unlike the main data, this dataset quantifies ribosomal paralog groups (4–6% of the measured mass). They are not model enzymes, so they lower f.

## Core ribosome

The 48 subunits of the main run, except at seven positions where the measured paralog takes the main run's member's place (`main_uniprot` gives the one it replaces). At five of them this dataset has the other paralog 5–290× more abundant, where the main data had the main run's member more abundant: RPL26 (288×), RPL9 (23×), RPS1 (8.4×), RPL21 (8.1×) and RPL8 (5.1×). At the other two the main run's member lacks a measurement in some samples, so the pipeline gives it no average: RPS19B is missing one sample (RPS19A is 5.1× more abundant), and RPL33B is in one sample of nine, where it is 7× more abundant than RPL33A. RPL33A averages 7.6e-6 mmol/gDW, below the core threshold of 1e-5, and is the only swapped subunit the model raises (1.6–4.1×). RPS10B is not in this dataset; neither RPS10 paralog was detected.

Re-selected by abundance on this dataset the core has 80 subunits (paralog groups, both paralogs of RPL6 and RPL16, and the ribosome-associated ZUO1, TEF4 and EGD1), and mu033 then needs 216 mg/gDW of a 199 mg/gDW pool.
