# Data

| File | Contents |
|---|---|
| `abs_proteomics.txt` | Absolute protein abundances [mmol/gDW]. 2929 rows, one header; three biological replicates per condition and four for hGR. Missing measurements are `NA`. |
| `fermentationData.txt` | Per condition: total protein content [g/gDW], dilution rate [1/h] and measured exchange rates [mmol/gDW/h]. An undetected byproduct is `NA`. |
| `ribosome.txt` | Ribosomal subunits with gene names, mass and sequence, from UniProt. |
| `selectedAnnotation.txt` | Assignment of proteins to the systems used when summarising enzyme usage. |
| `uniprot.tsv` | Yeast proteome with molecular masses, used to convert abundances from mmol/gDW to mg/gDW. |

## abs_proteomics.txt

`P61830` appears on two rows, once per gene of the duplicated histone
(`YBR010W`, `YNL031C`), with identical abundances. Readers collapse it to one
entry and fail if the duplicated rows ever disagree.

## uniprot.tsv

From [SysBioChalmers/GECKO](https://github.com/SysBioChalmers/GECKO),
`tutorials/full_ecModel/data/uniprot.tsv`, branch `develop4`, commit
`cbc4ca3311df8ca44c28792f40008b9cbef692bd`. UniProt proteome UP000002311,
reviewed entries. Covers every protein in `abs_proteomics.txt`.
