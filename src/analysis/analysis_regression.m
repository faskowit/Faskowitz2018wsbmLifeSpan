clc
clearvars

config_file='config_scale125.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

outIntermPrefix = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR);
outProcessPrefix = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR);

load(strcat(outIntermPrefix,'_comVecs.mat'))
load(strcat(outIntermPrefix,'_templateModel_1.mat'))
load(strcat(outProcessPrefix,'_basicData_v7p3.mat'))

%% read in data
% if we gonna use this script as a function 

% already have these defined
% dataStruct = readData.dataRaw ;
% datasetDemo = readData.demoRaw ;

%% gather data
% this is hidden in this function to make the main analysis a smoother ride
% function gatherStruct = ma_gather_data2(dataStruct , datasetDemo, ...
%   communityLabels, selectNodesFrmRaw, derivBool, thr, matDataStr)

% wsbm with ya template %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ca_wsbm = community_assign(templateModel) ;

wsbm_gatherStruct = ma_gather_data2(dataStruct, datasetDemo, ...
    ca_wsbm(:,2), selectNodesFrmRaw, 0, 1, 'countVolNormMat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% modularity matching num blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modResults = compute_mod_sweepG(templateModel) ;

idx = max(modResults(2:end,:)) == templateModel.R_Struct.k ;
modResultsSubS = modResults(:,idx) ;

% lets loop over the results to find the modular parition that is least
% distant from the WSBM parition-- for a fair comparison. 
modDist2Wsbm = zeros([ size(modResultsSubS,2) 1 ]);
for idx=1:size(modResultsSubS,2)
    
    modDist2Wsbm(idx) = partition_distance(ca_wsbm(:,2), ...
        modResultsSubS(2:end,idx));
    
end

[ ~ , minIdx ] = min(modDist2Wsbm) ; 

% align labels
ca_mod = modResultsSubS(2:end,minIdx) ;
ca_mod = CBIG_HungarianClusterMatch(ca_wsbm(:,2),ca_mod);
ca_mod = [ (1:(size(modResults,1)-1))' ca_mod ] ;

mod_gatherStruct = ma_gather_data2(dataStruct, datasetDemo, ...
    ca_mod(:,2), selectNodesFrmRaw, 0, 1, 'countVolNormMat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % yeo communities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('data/external/seven_network.mat')
% % add the temporal parietal to default mode
% mu_yeo7 = seven_network ;
% mu_yeo7(7,57) = 1 ;
% mu_yeo7(7,114) = 1 ;
% ca_yeo = community_assign(mu_yeo7);
% 
% % nope
% % % align labels
% %ca_yeo(:,2) = CBIG_HungarianClusterMatch(ca_wsbm(:,2),ca_yeo(:,2));
% 
% yeo_gatherStruct = ma_gather_data2(dataStruct, datasetDemo, ...
%     ca_yeo(:,2), selectNodesFrmRaw, 0, 1, 'countVolNormMat');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% output the comVecs file

comVecs = struct() ;
comVecs.wsbm = ca_wsbm(:,2) ;
comVecs.mod = ca_mod(:,2) ;
% comVecs.yeo = ca_yeo(:,2) ;
% 
outName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
save(outName,...
    'comVecs',...
    '-v7.3')

%% REGRESSIONS

fits = {'linear' 'quadratic' 'poisson' } ;
numBlocks = templateModel.R_Struct.k ;

% these get passed into regression func
otherArgs = {1 500 [] [] };

% predictor
X = datasetDemo.age ;
% predicted
% Y = rdens --> rdens in gather struct

% regressions with sex + totDen as covariates
wsbm_rdens_statMat = run_regression_over_yMat(X,wsbm_gatherStruct.rdens,fits,otherArgs) ;
mod_rdens_statMat = run_regression_over_yMat(X,mod_gatherStruct.rdens,fits,otherArgs) ;
% yeo_rdens_statMat = run_regression_over_yMat(X,yeo_gatherStruct.rdens,fits,otherArgs) ;

% regressions with sex + totDen + movement as covariates
wsbm_rmdens_statMat = run_regression_over_yMat(X,wsbm_gatherStruct.rmdens,fits,otherArgs) ;
mod_rmdens_statMat = run_regression_over_yMat(X,mod_gatherStruct.rmdens,fits,otherArgs) ;
% yeo_rmdens_statMat = run_regression_over_yMat(X,yeo_gatherStruct.rmdens,fits,otherArgs) ;

% %% dont run permuations for these tests to speed up
% 
otherArgs = {1 500 [] [] };

% regressions with just dens
wsbm_dens_statMat = run_regression_over_yMat(X,wsbm_gatherStruct.dens,fits,otherArgs) ;
mod_dens_statMat = run_regression_over_yMat(X,mod_gatherStruct.dens,fits,otherArgs) ;
% yeo_dens_statMat = run_regression_over_yMat(X,yeo_gatherStruct.dens,fits,otherArgs) ;

% regressions with eb density + sex + totDen as covariates
wsbm_rdenseb_statMat = run_regression_over_yMat(X,wsbm_gatherStruct.rdenseb,fits,otherArgs) ;
mod_rdenseb_statMat = run_regression_over_yMat(X,mod_gatherStruct.rdenseb,fits,otherArgs) ;
% yeo_rdenseb_statMat = run_regression_over_yMat(X,yeo_gatherStruct.rdenseb,fits,otherArgs) ;

%% save it

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_statMat_v7p3.mat');
save(outName,...
    '*_statMat',...
    '*_gatherStruct',...
    'templateModel',...
    'comVecs',...
    '-v7.3')
