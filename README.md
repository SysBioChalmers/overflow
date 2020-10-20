[![Zenodo DOI](https://zenodo.org/badge/288161152.svg)](https://zenodo.org/badge/latestdoi/288161152)


## Instructions
The `code` folder contains a number of scripts to reconstruct the models and replicate the analysis. Note that the user should have [RAVEN 2.4.0](https://github.com/SysBioChalmers/RAVEN/releases) or later installed.

- `prepareEnvironment.m`: run this script to prepare MATLAB for model reconstruction and analysis. This includes (1) cloning the correct GECKO version; (2) if required cloning yeast-GEM and reconstructing ecYeast-GEM 8.1.3 [also already distributed with this repository]; (3) load all measured flux and proteomics data; (4) set parameters for GECKO.
- `generateProtModels.m`: generate the five condition specific proteome-constrained ec-models, and store these in the models subdirectory. In `results/modelGeneration` is written which enzyme abundances were flexibilized to allow growth at the measured dilution rate. The models are constrained for glucose uptake, growth rate is set to dilution rate and the objective function is set to minimization of the unmeasured protein pool usage. Results from this FBA is written to `results/modelSimulation`.
- `ribosome.m`: add ribosomal subunits to the ec-models. In `results/modelGeneration` is plotted the average abundances of the ribosomal subunits, to identify the "core" ribosome to be included in the model. Also written in this subdirectory is which subunit abundances were flexibilized to allow growth at the measured dilution rate.
- `analyzeUsage.m`: run same FBA as in `generateProtModels` and summarize the enzyme usages in `results/enzymeUsage`, where the capacity and absolute usages are stored separately, together with a summary for each subsystem.
