% This function prepares the environment to run the analysis, by cloning
% the correct repositories and setting appropriate parameters

code = pwd();

%% GECKO, yeast-GEM and ecModel
% Make sure that GECKO is installed
try
    GECKOInstaller.install
catch
    error('Cannot find GECKO, find installation instructions <a href = "https://github.com/SysBioChalmers/GECKO/wiki/Installation-and-upgrade">here</a>.')
end
ModelAdapter = ModelAdapterManager.setDefault('ecYeastOverflowAdapter.m');
ecModel = loadEcModel();

ecModelb = bayesianSensitivityTuning(ecModel);

%% ============== Prepare data ============== %%
%% Load proteomics data
protData = loadProtData([3,3,3,3,4],'',false);
protData = loadProtData([3,3,3,3,4],'',true);
uniprot  = loadDatabases('uniprot');
uniprot  = uniprot.uniprot;
% Change from gene to Uniprot identifiers
[a,b] = ismember(protData.uniprotIDs,uniprot.genes);
protData.uniprotIDs(a) = uniprot.ID(b(a));
% A few manual assignments
[a,~] = ismember(protData.uniprotIDs,{'YBR009C','YPR080W','YKR059W','YDR385W','YBR010W','YNL031C'});
protData.uniprotIDs(a) = {'P02309','P02994','P10081','P32324','P61830','P61830'};
% Correct from mmol/gDCW to mg/gDCW
[a,b] = ismember(protData.uniprotIDs,uniprot.ID);
protData.abundances = uniprot.MW(b).*protData.abundances;

%% Load total protein content and fermentation data
fluxData = loadFluxData();

%% Load additional parameters
cd(code)
%Set some additional parameters
oxPhos = find(startsWith(ecModel.rxns,{'r_1021','r_0439','r_0438','r_0226'}));
%Get indexes for carbon source uptake and biomass pseudoreactions
positionsEC(1) = getIndexes(ecModel,params.c_source,'rxns');
positionsEC(2) = getIndexes(ecModel,params.bioRxn,'rxns');
clear ans fID data byProds fileNames GECKO_path i fID
