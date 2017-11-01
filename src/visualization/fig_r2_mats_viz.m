clc
clearvars

%% load up the data...

load('data/processed/yeo_both_normalpoisson_a0p5_statMat_v7p3.mat')
load('data/interim/yeo_both_normalpoisson_a0p5_templateModel_1.mat')

%%

fits = { 'linear' 'quadratic' 'possion' } ;
numBlocks = templateModel.R_Struct.k ;
bonf_limit = 0.05 / ( sum(sum(triu(ones(numBlocks)))) * length(fits) ) ;

%% regressions on wsbm rdens 

statMapLookinAt = yeo_rdens_statMat ;
%statMapLookinAt = mod_rdens_statMat ;
%numBlocks = templateModel.R_Struct.k ;
numBlocks = size(yeo_rdens_statMat,1);

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

%% viz it

% function [] = viz_statMat(mat,colorbar_bool,colorbar_axis_vals,colorbar_title,plot_title,midValue)
sp1 = subplot(1,3,1) ;
viz_statMat(rdens_r2_thr(:,:,1),0,[0 0.3],[],'Linear R-squared',0.15) ;
axis square

sp2 = subplot(1,3,2) ;
viz_statMat(rdens_r2_thr(:,:,2),0,[0 0.3],[],'Quadratic R-squared',0.15) ;
axis square

sp3 = subplot(1,3,3) ; 
viz_statMat(rdens_r2_thr(:,:,3),1,[0 0.3],[],'Poisson R-squared',0.15) ;
axis square

% take care of the resize on the third plot
sz2 = get(sp2,'position');
sz3 = get(sp3,'position');
sz3(3:4) = sz2(3:4) ;
set(sp3,'position',sz3)

% [left, bottom, width, height]
p2 = get(sp2,'pos') ;
p2(1) = p2(1) * 0.91 ;
set(sp2,'pos',p2)
p3 = get(sp3,'pos') ;
p3(1) = p3(1) * 0.895 ;
set(sp3,'pos',p3)

%% viz individual connections...

block1 = 4 ;
block2 = 8 ;
regressLayer = 2;

statMapLookinAt = wsbm_rdens_statMat ;
statEntry = statMapLookinAt{block1,block2,regressLayer};

viz_blockRegress(statEntry,1);

xlabel('age')
ylabel('adjusted streamline density')



