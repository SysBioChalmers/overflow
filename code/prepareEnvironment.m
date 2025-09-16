% This function prepares the environment to run the analysis, by cloning
% the correct repositories and setting appropriate parameters

code = pwd();

%% GECKO and ecModel
% Make sure that GECKO is installed
try
    GECKOInstaller.install
catch
    error('Cannot find GECKO, find installation instructions <a href = https://github.com/SysBioChalmers/GECKO/wiki/Installation-and-upgrade>here</a>.')
end

adapterLocation = '../ecYeastOverflowAdapter.m';
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.params;

ecModel = loadEcModel();

%% Prepare data
%Load proteomics data
protData = loadProtData([3,3,3,3,4],'',false);
sum(protData.abundances(:,1),'omitnan')
protData = loadProtData([3,3,3,3,4],'',true);
uniprot  = loadDatabases('uniprot');
uniprot  = uniprot.uniprot;

[a,b] = ismember(protData.uniprotIDs,uniprot.genes);

protData.uniprotIDs(a) = uniprot.ID(b(a));
% A few manual assignments
[a,b] = ismember(protData.uniprotIDs,{'YBR009C','YPR080W','YKR059W','YDR385W','YBR010W','YNL031C'});
protData.uniprotIDs(a) = {'P02309','P02994','P10081','P32324','P61830','P61830'};

% Correct from mmol/gDCW to mg/gDCW
[a,b] = ismember(protData.uniprotIDs,uniprot.ID);
protData.abundances = uniprot.MW(b).*protData.abundances;


% fID       = fopen('../data/abs_proteomics.txt');
% prot.cond = textscan(fID,['%s' repmat(' %s',1,17)],1);
% prot.data = textscan(fID,['%s %s' repmat(' %f',1,16)],'TreatAsEmpty',{'NA','na','NaN'});
% prot.cond = [prot.cond{3:end}];
% prot.IDs  = prot.data{1};
% prot.data = cell2mat(prot.data(3:end));
% fclose(fID);

%Load total protein content and fermentation data
fluxData = loadFluxData();
% 
% fID       = fopen('../data/fermentationData.txt');
% byProds   = textscan(fID,['%s' repmat(' %s',1,9)],1,'Delimiter','\t');
% data      = textscan(fID,['%s' repmat(' %f',1,9)],'TreatAsEmpty',{'NA','na','NaN'});
% fclose(fID);
% 
% flux.conds     = data{1};
% flux.Ptot      = data{2};
% flux.Drate     = data{3};
% flux.GUR       = data{4};
% flux.CO2prod   = data{5};
% flux.OxyUptake = data{6};
% flux.byP_flux  = [data{7:end}];
% flux.byP_flux(isnan(flux.byP_flux))=0;
% flux.byProds   = [byProds{7:end}];

%% Load additional parameters
cd(code)

%Set some additional parameters
oxPhos = ecModel.rxns(startsWith(ecModel.rxns,params.oxPhos));
grouping=[3 3 3 3 4];
clear repl
for i=1:length(grouping)
    try
        repl.first(i)=repl.last(end)+1;
        repl.last(i)=repl.last(end)+grouping(i);
    catch
        repl.first=1;
        repl.last=grouping(1);
    end
end
%Get indexes for carbon source uptake and biomass pseudoreactions
positionsEC(1) = getIndexes(ecModel,params.c_source,'rxns');
positionsEC(2) = getIndexes(ecModel,params.bioRxn,'rxns');
clear ans fID data byProds fileNames GECKO_path i fID