function regression_analysis( readData , templateModel , selectNodesFrmRaw)
% run the regression analysis here
% make into function to be called by larger script so we can analyze
% multiple community structures

%% STEPS
% 1) fit the template to brains, measure properites across
% lifespan; assume that nodes belong to same blocks across the lifespan
% 2) fit linear vs quadratic vs something else to data for each block
% 3) compare to modularity communities

% for other script?
% 4) bin into age groups
% 5) individual fits

% add BCT path
addpath('~/JOSHSTUFF/scripts/BCT/2017_01_15_BCT/')
%addpath('~/JOSHSTUFF/projects/SBM2/aux_stuff/')

%% read in data
% if we gonna use this script as a function 

dataStruct = readData.dataRaw ;
datasetDemo = readData.demoRaw ;

%% gather data
% this is hidden in this function to make the main analysis a smoother ride
% function gatherStruct = ma_gather_data2(dataStruct , datasetDemo, ...
%   communityLabels, selectNodesFrmRaw, derivBool, thr, matDataStr)

% wsbm with ya template %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ca_wsbm = community_assign(templateModel) ;

wsbm_gatherStruct = ma_gather_data2(dataStruct, datasetDemo, ...
    ca_wsbm(:,2), selectNodesFrmRaw, 1, 1, 'countVolNormMat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we should actually only apply this if we assigned subcort to different
% communites in the fitting process
%
% % wsbm with ya template + subcort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ca_wsbm_wsubcort = vertcat(ca_wsbm(:,2), ones([ 14 1 ]) .* 8); 
% selectNodesFrmRaw_wSubcort = [ selectNodesFrmRaw 117:130 ] ;
% 
% gatherStruct_wSubcort = ma_gather_data2(dataStruct, datasetDemo, ...
%     ca_wsbm_wsubcort, selectNodesFrmRaw_wSubcort, 0, 1, 'countVolNormMat');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% sanity check
% for idx=2:size(modResultsSubS,2)
%    
%     modResultsSubS(:,idx) = CBIG_HungarianClusterMatch(modResultsSubS(:,1),...
%         modResultsSubS(:,idx));
% end

% [ ~ , idx ] = max(modResults(1,idx)) ; 

% align labels
ca_mod = modResultsSubS(2:end,minIdx) ;
ca_mod = CBIG_HungarianClusterMatch(ca_wsbm(:,2),ca_mod);
ca_mod = [ (1:(size(modResults,1)-1))' ca_mod ] ;

mod_gatherStruct = ma_gather_data2(dataStruct, datasetDemo, ...
    ca_mod(:,2), selectNodesFrmRaw, 0, 1, 'countVolNormMat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% yeo communities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('data/external/seven_network.mat')
% add the temporal parietal to default mode
mu_yeo7 = seven_network ;
mu_yeo7(7,57) = 1 ;
mu_yeo7(7,114) = 1 ;
ca_yeo = community_assign(mu_yeo7);

% nope
% % align labels
%ca_yeo(:,2) = CBIG_HungarianClusterMatch(ca_wsbm(:,2),ca_yeo(:,2));

yeo_gatherStruct = ma_gather_data2(dataStruct, datasetDemo, ...
    ca_yeo(:,2), selectNodesFrmRaw, 0, 1, 'countVolNormMat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% output the comVecs file

outName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');



%% STEP 1

fits = {'linear' 'quadratic' 'poisson' } ;
%fits = {'linear' 'quadratic' } ;

numBlocks = templateModel.R_Struct.k ;

% these get passed into regression func
otherArgs = {1 500 [] 10000 };

% predictor
X = datasetDemo.age ;
% predicted
% Y = rdens --> rdens in gather struct

wsbm_rdens_statMat = run_regression_over_yMat(X,wsbm_gatherStruct.rdens,fits,otherArgs) ;
mod_rdens_statMat = run_regression_over_yMat(X,mod_gatherStruct.rdens,fits,otherArgs) ;
yeo_rdens_statMat = run_regression_over_yMat(X,yeo_gatherStruct.rdens,fits,otherArgs) ;

comVecs = struct() ;
comVecs.wsbm = ca_wsbm(:,2) ;
comVecs.mod = ca_mod(:,2) ;
comVecs.yeo = ca_yeo(:,2) ;

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_statMat_v7p3.mat');
save(outName,...
    'wsbm_rdens_statMat',...
    'mod_rdens_statMat',...
    'yeo_rdens_statMat',...
    'wsbm_gatherStruct',...
    'mod_gatherStruct',...
    'yeo_gatherStruct',...
    'templateModel',...
    'datasetDemo',...
    'dataStruct',...
    'comVecs',...
    'templateModel',...
    '-v7.3')

%% STEP 2
%figure out which model best based on RMSE

% bonferoni threshold based on comparisons made in upper triangle
bonf_limit = 0.05 / ( sum(sum(triu(ones(numBlocks)))) * length(fits) ) ;

%% regressions on wsbm rdens 

statMapLookinAt = wsbm_rdens_statMat ;
numBlocks = templateModel.R_Struct.k ;

rdens_xValR2 = [] ;
rdens_rmse = [] ;
rdens_pval = [] ;
rdens_r2 = [] ;

for idx=1:length(fits)

    rdens_xValR2(:,:,idx) =  get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalR2'});
    rdens_rmse(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalRMSE'});
    rdens_pval(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'permStruct' 'permPvalR2'});
    rdens_r2(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'lsFitStruct' 'R2'});
    
end
   
% get the min in the third dim
[~,min_ind] = min(rdens_rmse,[],3);

% get ind where full fit is significant
rdens_pval_thresh = rdens_pval< bonf_limit ;

% get a thresholded
% at each dim, multiply by pval thr
rdens_r2_thr = zeros([numBlocks numBlocks length(fits)]);
for idx=1:length(fits)
   
    % at this dim (the fit) get ind where this fit was lowest rmse
    ind2 = min_ind == idx ;
    
    rdens_r2_thr(:,:,idx) =  rdens_pval_thresh(:,:,idx) .* ...
        rdens_r2(:,:,idx) .* ind2;
    
end

%% regressions on mod rdens 

statMapLookinAt = mod_rdens_statMat ;

rdens_rmse = [] ;
rdens_pval = [] ;
rdens_r2 = [] ;

for idx=1:length(fits)

    rdens_rmse(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalRMSE'});
    rdens_pval(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'permStruct' 'permPvalR2'});
    rdens_r2(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'lsFitStruct' 'R2'});
    
end
   
% get the min in the third dim
[~,min_ind] = min(rdens_rmse,[],3);

% get ind where full fit is significant
rdens_pval_thresh = rdens_pval< bonf_limit ;

% get a thresholded
% at each dim, multiply by pval thr
rdens_r2_thr = zeros([numBlocks numBlocks length(fits)]);
for idx=1:length(fits)
   
    % at this dim (the fit) get ind where this fit was lowest rmse
    ind2 = min_ind == idx ;
    
    rdens_r2_thr(:,:,idx) =  rdens_pval_thresh(:,:,idx) .* ...
        rdens_r2(:,:,idx) .* ind2;
    
end

%% regressions on yeo rdens 

statMapLookinAt = yeo_rdens_statMat ;
numBlocks = size(yeo_gatherStruct.rdens,1);

rdens_rmse = [] ;
rdens_pval = [] ;
rdens_r2 = [] ;

for idx=1:length(fits)

    rdens_rmse(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalRMSE'});
    rdens_pval(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'permStruct' 'permPvalR2'});
    rdens_r2(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'lsFitStruct' 'R2'});
    
end
   
% get the min in the third dim
[~,min_ind] = min(rdens_rmse,[],3);

% get ind where full fit is significant
rdens_pval_thresh = rdens_pval< bonf_limit ;

% get a thresholded
% at each dim, multiply by pval thr
rdens_r2_thr = zeros([numBlocks numBlocks length(fits)]);

for idx=1:length(fits)
   
    % at this dim (the fit) get ind where this fit was lowest rmse
    ind2 = min_ind == idx ;
    
    rdens_r2_thr(:,:,idx) =  rdens_pval_thresh(:,:,idx) .* ...
        rdens_r2(:,:,idx) .* ind2;
    
end
