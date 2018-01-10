%% lets analyze the age bin scripts
clc
clearvars

addpath('~/JOSHSTUFF/scripts/BCT/2017_01_15_BCT/')

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

% % load the data we need to analyze this ish. 
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_fit_wsbm_script_v7p3.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_templateModel_1.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_agebin_template_results_v7p3.mat');
load(loadName) ;

%% initial data stuff

% number of actual groups is thresholds + 1
thresholds = 4 ;
age_bins = thresholds + 1;

[agesSorted,ageSortIdx] = sort(datasetDemo.age) ;

low_quantile = [ 0 quantile(1:length(datasetDemo.age),4) ] ; 
high_quantile = [ quantile(1:length(datasetDemo.age),4) (length(datasetDemo.age)+1) ] ;

templateIdMat = zeros( [ size(datasetDemo.age,1) (thresholds+1) ] ) ;

for templateIdx=1:age_bins

    subjectsIdxVec = ageSortIdx > low_quantile(templateIdx) & ...
        ageSortIdx < high_quantile(templateIdx) ;
    
    templateIdMat(:,templateIdx) = ~~subjectsIdxVec ;
end

%% first align all agebin templates to the template model

[~,ca_wsbm] = community_assign(templateModel) ;
agebin_ca = zeros([ templateModel.Data.n age_bins ]);

for idx=1:age_bins

    [~,tmpCA] = community_assign(consensusCrossIdx{idx}); 
    %first align communsensus parition to template model...
    agebin_ca(:,idx) = CBIG_HungarianClusterMatch(ca_wsbm,tmpCA);
    
end

%% some analysis 
% can we comparse vi between subj and agebin template to show how
% varability in community strucutre increases over age ?

VItoAgeBin = cell([age_bins 1]) ;
nNodes = templateModel.Data.n ;

for idx = 1:age_bins
   
    ageBinData = fitWSBMAllStruct(~~templateIdMat(:,idx));
    ageBinTemplate = consensusCrossIdx{idx} ;
    
    VItoAgeBin{idx} = zeros([length(ageBinData) 1]);
    
    % now measure the VI from each subject to its agebin template
    for jdx = 1:length(ageBinData);
        
        %tmpAdj = ageBinData(jdx).centModel.Data.Raw_Data ;
        %tmpAdj(isnan(tmpAdj)) = 0;
        [~,tmpCA] = community_assign(ageBinData(jdx).centModel);
        
        VItoAgeBin{idx}(jdx) = partition_distance(agebin_ca(:,idx),tmpCA);
            
    end
end









