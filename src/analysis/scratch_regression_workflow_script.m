% regression workflow script

%load('workspaces/fit_wsbm_script_workspace.mat')

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

%% read in data
% if we gonna use this script as a function 

% dataStruct = readData.dataRaw ;
% datasetDemo = readData.demoRaw ;

%% gather data
% this is hidden in this function to make the main analysis a smoother ride
% function gatherStruct = ma_gather_data2(dataStruct , datasetDemo, ...
%   communityLabels, selectNodesFrmRaw, derivBool, thr, matDataStr)

% wsbm with ya template %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ca_wsbm = community_assign(templateModel) ;

wsbm_gatherStruct = ma_gather_data2(dataStruct, datasetDemo, ...
    ca_wsbm(:,2), selectNodesFrmRaw, 1, 1, 'countVolNormMat');

wsbm_gatherStruct_othr = ma_gather_data2(dataStruct, datasetDemo, ...
    ca_wsbm(:,2), selectNodesFrmRaw, 1, 1, 'countVolLenNormMat');

% wsbm_gatherStruct_othr2 = ma_gather_data2(dataStruct, datasetDemo, ...
%     ca_wsbm(:,2), selectNodesFrmRaw, 1, 1, 'countMat');

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
[ ~ , modResults ] = compute_mod(templateModel,1.0) ;
idx = max(modResults(2:end,:)) == templateModel.R_Struct.k ;
modResultsSubS = modResults(:,idx) ;
[ ~ , idx ] = max(modResults(1,idx)) ; 
ca_mod = [ (1:(size(modResults,1)-1))' modResultsSubS(2:end,idx) ] ;

mod_gatherStruct = ma_gather_data2(dataStruct, datasetDemo, ...
    ca_mod(:,2), selectNodesFrmRaw, 0, 1, 'countVolNormMat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% yeo communities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('aux_stuff/seven_network.mat')
% add the temporal parietal to default mode
mu_yeo7 = seven_network ;
mu_yeo7(7,57) = 1 ;
mu_yeo7(7,114) = 1 ;
ca_yeo = community_assign(mu_yeo7);
yeo_gatherStruct = ma_gather_data2(dataStruct, datasetDemo, ...
    ca_yeo(:,2), selectNodesFrmRaw, 0, 1, 'countVolNormMat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
wsbm_rdens_statMat_othr = run_regression_over_yMat(X,wsbm_gatherStruct_othr.rdens,fits,otherArgs) ;

% save this stuff
addNote = '' ;
outName = strcat(projectDir,'/analysis/', addNote, 'regression_statMats.mat');
save(outName, ...
    'wsbm_rdens_statMat', ...
    'mod_rdens_statMat', ...
    'yeo_rdens_statMat', ...
    'wsbm_rdens_statMat_othr' ...
)

%% STEP 2
%figure out which model best based on RMSE

% bonferoni threshold based on comparisons made in upper triangle
bonf_limit = 0.05 / sum(sum(triu(ones(numBlocks)))) ;

%% regressions on wsbm rdens 

statMapLookinAt = wsbm_rdens_statMat ;
numBlocks = templateModel.R_Struct.k ;

rdens_rmse = [] ;
rdens_pval = [] ;
rdens_r2 = [] ;
rdens_xvalr2 = [] ;

for idx=1:length(fits)

    rdens_rmse(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalRMSE'});
    rdens_pval(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'permStruct' 'permPvalR2'});
    rdens_r2(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'lsFitStruct' 'R2'});
    
    rdens_xvalr2(:,:,idx) = get_frm_statMat(statMapLookinAt(:,:,idx),{'xvalR2'});
    
end
   
% get the min in the third dim
[~,min_ind] = min(rdens_rmse,[],3);

% get ind where full fit is significant
rdens_pval_thresh = rdens_pval< bonf_limit ;

% get a thresholded
% at each dim, multiply by pval thr
rdens_r2_thr = zeros([numBlocks numBlocks length(fits)]);
rdens_xvalr2_thr = zeros([numBlocks numBlocks length(fits)]);

for idx=1:length(fits)
   
    % at this dim (the fit) get ind where this fit was lowest rmse
    ind2 = min_ind == idx ;
    
    rdens_r2_thr(:,:,idx) =  rdens_pval_thresh(:,:,idx) .* ...
        rdens_r2(:,:,idx) .* ind2;
    
    rdens_xvalr2_thr(:,:,idx) =  rdens_pval_thresh(:,:,idx) .* ...
        rdens_xvalr2(:,:,idx) .* ind2;
    
end

%% viz it
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

%% look at a specific block vs. block interaction
block1 = 3;
block2 = 4 ;
regressLayer = 1;

statMapLookinAt = wsbm_rdens_statMat ;
statEntry = statMapLookinAt{block1,block2,regressLayer};

viz_blockRegress(statEntry,0);

% 
% xVec = min(statEntry.coef.x):0.25:...
%     max(statEntry.coef.x);
% 
% yhat = polyval(statEntry.coef.full,xVec);
% 
% for kk = 1:size(statEntry.coef.boot,1)
%     bootCI(kk,:) = polyval(statEntry.coef.boot(kk,:),xVec);
% %     % and confidence interval for vertex
%     vCI(kk) = -(statEntry.coef.boot(kk,2)./(2*statEntry.coef.boot(kk,1)));
% end
% 
% p95 = prctile(bootCI,[2.5,97.5]);
% v95 = prctile(vCI,[2.5,97.5]);
% 
% figure;
% 
% %plot(statEntry.coef.x, statEntry.coef.y, 'o','color',[0 0 0],'markersize',6);
% plot(statEntry.coef.x, statEntry.coef.y, 'o','markersize',6);
% hold on;
% 
% fill([xVec fliplr(xVec)],[p95(1,:) fliplr(p95(2,:))],[ 0.75 0.75 0.75 ],'facealpha',.5,'edgealpha',0);
% 
% plot(xVec,yhat, 'color',[ 0.1 0.1 0.1 ],'linewidth',2);
% %plot(xVec,yhat,'linewidth',2);
% 
% % vertex pos
% % put it at the bottom essentially 
% ylims = ylim() ;
% plot(v95,repmat(ylims(1) .* 1.05,1,2),'linewidth',3,'color',[ 0.9 0.9 0.9 ]);
% % set ythe ylim back to befor
% ylim([ylims(1)*1.075 ylims(2)]);

%% regressions on wsbm with countVolLenNorMat rdens 

statMapLookinAt = wsbm_rdens_statMat_othr ;
numBlocks = templateModel.R_Struct.k ;

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

%% viz it
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

%% look at a specific block vs. block interaction
block1 = 8;
block2 = 8 ;
regressLayer = 2;

statMapLookinAt = wsbm_rdens_statMat_othr ;
statEntry = statMapLookinAt{block1,block2,regressLayer};

viz_blockRegress(statEntry,1);

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

%% viz it
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

%% regressions on yeo rdens 

statMapLookinAt = yeo_rdens_statMat ;
numBlocks = size(yeo_rdens_statMat,1);

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

%% viz it
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
