%% RESULTS
% results from looking at the edge weights in static communities across the
% lifespan

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

% load up results
load('data/processed/yeo_both_normalpoisson_a0p5_statMat_v7p3.mat')

%% STEP 2
%figure out which model best based on RMSE

fits = {'linear' 'quadratic' 'poisson' } ;
numBlocks = templateModel.R_Struct.k ;

% bonferoni threshold based on comparisons made in upper triangle
bonf_limit = 0.05 / ( sum(sum(triu(ones(numBlocks)))) * length(fits) ) ;

%% regressions on wsbm rdens 

statMapLookinAt = wsbm_rdens_statMat ;

wsbm_rdens_xValR2 = zeros([ numBlocks numBlocks length(fits)]) ;
wsbm_rdens_rmse = zeros([ numBlocks numBlocks length(fits)]) ;
wsbm_rdens_pval = zeros([ numBlocks numBlocks length(fits)]) ;
wsbm_rdens_r2 = zeros([ numBlocks numBlocks length(fits)]) ;
wsbm_rdens_r2_thr = zeros([numBlocks numBlocks length(fits)]);

for idx=1:length(fits)

    wsbm_rdens_xValR2(:,:,idx) =  get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalR2'});
    wsbm_rdens_rmse(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalRMSE'});
    wsbm_rdens_pval(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'permStruct' 'permPvalR2'});
    wsbm_rdens_r2(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'lsFitStruct' 'R2'});
    
end
   
% get the min in the third dim
[~,min_ind] = min(wsbm_rdens_rmse,[],3);

% get ind where full fit is significant
rdens_pval_thresh = wsbm_rdens_pval< bonf_limit ;

% get a thresholded
% at each dim, multiply by pval thr
for idx=1:length(fits)
   
    % at this dim (the fit) get ind where this fit was lowest rmse
    ind2 = min_ind == idx ;
    
    wsbm_rdens_r2_thr(:,:,idx) =  rdens_pval_thresh(:,:,idx) .* ...
        wsbm_rdens_r2(:,:,idx) .* ind2;
    
end

%% regressions on mod rdens 

statMapLookinAt = mod_rdens_statMat ;

mod_rdens_xValR2 = zeros([ numBlocks numBlocks length(fits)]) ;
mod_rdens_rmse = zeros([ numBlocks numBlocks length(fits)]) ;
mod_rdens_pval = zeros([ numBlocks numBlocks length(fits)]) ;
mod_rdens_r2 = zeros([ numBlocks numBlocks length(fits)]) ;
mod_rdens_r2_thr = zeros([numBlocks numBlocks length(fits)]);

for idx=1:length(fits)

    mod_rdens_xValR2(:,:,idx) =  get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalR2'});
    mod_rdens_rmse(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalRMSE'});
    mod_rdens_pval(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'permStruct' 'permPvalR2'});
    mod_rdens_r2(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'lsFitStruct' 'R2'});
    
end
   
% get the min in the third dim
[~,min_ind] = min(mod_rdens_rmse,[],3);

% get ind where full fit is significant
rdens_pval_thresh = mod_rdens_pval< bonf_limit ;

% get a thresholded
% at each dim, multiply by pval thr
for idx=1:length(fits)
   
    % at this dim (the fit) get ind where this fit was lowest rmse
    ind2 = min_ind == idx ;
    
    mod_rdens_r2_thr(:,:,idx) =  rdens_pval_thresh(:,:,idx) .* ...
        mod_rdens_r2(:,:,idx) .* ind2;
    
end
