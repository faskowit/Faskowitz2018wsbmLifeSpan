clc 
clearvars

config_file='config_scale125.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_templateModel_1.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
load(loadName) ;

% the stats maps for wsbm, mod, yeo
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_statMat_v7p3.mat');
load(loadName) ;

FIGURE_NAME = 'figB' ;

outputdir = strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/');
mkdir(outputdir)

% let's make picture of regression stuff
%% setup data

comStructNames = { 'wsbm' 'mod' 'yeo' } ;
fitNames = { 'linear' 'quadratic' 'possion' } ;
% 
% statMapAcrossMethod = cell([ 6 1 ]);
% statMapAcrossMethod{1} = wsbm_rdens_statMat ;
% statMapAcrossMethod{2} = mod_rdens_statMat ;
% statMapAcrossMethod{3} = yeo_rdens_statMat ;
% statMapAcrossMethod{4} = wsbm_rmdens_statMat ;
% statMapAcrossMethod{5} = mod_rmdens_statMat ;
% statMapAcrossMethod{6} = yeo_rmdens_statMat ;

% results w/ pvals
statMap_resultsXValR2 = cell([ 6 1]);

%% extract r2
for sm = 1:length(statMapAcrossMethod)
    
    currStatMap = statMapAcrossMethod{sm} ;

    nBlocks = size(currStatMap,1);
    bonf_limit = 0.05 / ( sum(sum(triu(ones(nBlocks)))) * length(fitNames) ) ;

    
    tmp_xValR2 = zeros([ nBlocks nBlocks ]) ;
    tmp_rmse = zeros([ nBlocks nBlocks ]) ;
    tmp_rdens_pval = zeros([ nBlocks nBlocks ]) ;
    tmp_rdens_r2 = zeros([ nBlocks nBlocks ]) ;

    for idx=1:length(fitNames)

        tmp_xValR2(:,:,idx) =  get_frm_statMat(currStatMap(:,:,idx),{'xvalR2'});
        tmp_rmse(:,:,idx) = get_frm_statMat(currStatMap(:,:,idx),{'xvalRMSE'});
        tmp_rdens_pval(:,:,idx) = get_frm_statMat(currStatMap(:,:,idx),{'permStruct' 'permPvalR2'});
        tmp_rdens_r2(:,:,idx) = get_frm_statMat(currStatMap(:,:,idx),{'lsFitStruct' 'R2'});

    end

    % get the min in the third dim
    [~,min_ind] = min(tmp_rmse,[],3);

    % get ind where full fit is significant
    tmp_pval_thresh = tmp_rdens_pval < bonf_limit ;

    % get a thresholded
    % at each dim, multiply by pval thr
    tmp_r2_thr = zeros([nBlocks nBlocks length(fitNames)]);
    for idx=1:length(fitNames)

        % at this dim (the fit) get ind where this fit was lowest rmse
        ind2 = min_ind == idx ;

        tmp_r2_thr(:,:,idx) =  tmp_pval_thresh(:,:,idx) .* ...
            tmp_xValR2(:,:,idx) .* ind2;

    end

    statMap_resultsXValR2{sm} = tmp_r2_thr ;
    
end

%% also gather stats from other regressions

% otherStatMapAcrossMethods = cell([6 1]) ;
% otherStatMapAcrossMethods{1} = wsbm_dens_statMat ;
% otherStatMapAcrossMethods{2} = mod_dens_statMat ;
% otherStatMapAcrossMethods{3} = yeo_dens_statMat ;
% otherStatMapAcrossMethods{4} = wsbm_rdenseb_statMat ;
% otherStatMapAcrossMethods{5} = mod_rdenseb_statMat ;
% otherStatMapAcrossMethods{6} = yeo_rdenseb_statMat ;

otherStatMapAcrossMethods = cell([2 1]) ;
otherStatMapAcrossMethods{1} = wsbm_rdens_statMat ;
otherStatMapAcrossMethods{2} = mod_rdens_statMat ;

% results no pvals
% otherStatMap_resultsXValR2 = cell([ 6 1]);
otherStatMap_resultsXValR2 = cell([ 2 1]);

for sm = 1:length(otherStatMapAcrossMethods)
    
    currStatMap = otherStatMapAcrossMethods{sm} ;

    nBlocks = size(currStatMap,1);
    
    % no bonf for this one
    %bonf_limit = 0.05 / ( sum(sum(triu(ones(nBlocks)))) * length(fits) ) ;
    
    tmp_xValR2 = zeros([ nBlocks nBlocks ]) ;
    tmp_rmse = zeros([ nBlocks nBlocks ]) ;

    for idx=1:length(fitNames)

        tmp_xValR2(:,:,idx) =  get_frm_statMat(currStatMap(:,:,idx),{'xvalR2'});
        tmp_rmse(:,:,idx) = get_frm_statMat(currStatMap(:,:,idx),{'xvalRMSE'});

    end

    % get the min in the third dim
    [~,min_ind] = min(tmp_rmse,[],3);

    tmp_thresh = tmp_xValR2 > 1 ;
    
    % get a thresholded
    % at each dim, multiply by pval thr
    tmp_r2_thr = zeros([nBlocks nBlocks length(fitNames)]);
    for idx=1:length(fitNames)

        % at this dim (the fit) get ind where this fit was lowest rmse
        ind2 = min_ind == idx ;

        tmp_r2_thr(:,:,idx) =  tmp_thresh(:,:,idx) .* ...
            tmp_xValR2(:,:,idx) .* ind2;

    end

    otherStatMap_resultsXValR2{sm} = tmp_r2_thr ;
    
end

%% look at them now
% rdens
comStructNames = {'wsbm' 'mod'} ;

for fig = 1:length(comStructNames)
%for fig = 1:1

    plotStats = otherStatMap_resultsXValR2{fig} ;
    
    plotNames = {'Linear R-squared' 'Quadratic R-squared' 'Poisson R-squared'};
    
    % because we are reporting the x-val R2 and these are 0-100
    plotStats = plotStats ./ 100 ;

    %% viz it
    figure
    
    subp = tight_subplot(1,3,[.01 .03],[.1 .01],[.02 .02]) ;
    % make it full screen
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);

    commVec = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.(comStructNames{fig}));
    
    for tsb_idx = 1:3

        axes(subp(tsb_idx))
%         viz_statsInMat_wComm(plotStats(:,:,tsb_idx),...
%             commVec,'Paired',10,...
%             0,[0 0.3],'Greys',plotNames{tsb_idx},0.15) ;

        viz_statsInMat_wComm(plotStats(:,:,tsb_idx),...
            commVec,'Paired',11,...
            0,[0 0.3],'Greys',plotNames{tsb_idx},0.15) ;
        
        set(gca,'XColor',[ 0 0 0 0.0001 ])
        set(gca,'YColor',[ 0 0 0 0.0001 ])
        ax = gca ;
        % ax.Color = [ 1 1 0.8 ] ;
        axis square

    end

    hold off
    
%     % save it
%     fileName = strcat(comStructNames{fig},'_rdens_regressR2.png');
%     ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
%     %set(gcf,'paperpositionmode','auto');
%     print(gcf,'-dpng','-r500',ff);
%     close(gcf)

end

%% rmdens
for fig = 1:length(comStructNames)
%for fig = 1:1

    plotStats = statMap_resultsXValR2{(fig+3)} ;
    
    plotNames = {'Linear R-squared' 'Quadratic R-squared' 'Poisson R-squared'};
    
    % because we are reporting the x-val R2 and these are 0-100
    plotStats = plotStats ./ 100 ;

    %% viz it
    figure
    
    subp = tight_subplot(1,3,[.01 .03],[.1 .01],[.02 .02]) ;
    % make it full screen
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);

    commVec = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.(comStructNames{fig}));
    
    for tsb_idx = 1:3

        axes(subp(tsb_idx))
        viz_statsInMat_wComm(plotStats(:,:,tsb_idx),...
            commVec,'Paired',10,...
            0,[0 0.3],'Greys',plotNames{tsb_idx},0.15) ;


        set(gca,'XColor',[ 0 0 0 0.0001 ])
        set(gca,'YColor',[ 0 0 0 0.0001 ])
        ax = gca ;
        % ax.Color = [ 1 1 0.8 ] ;
        axis square

    end

    hold off
    
    % save it
    fileName = strcat(comStructNames{fig},'_rmdens_regressR2.png');
    ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
    %set(gcf,'paperpositionmode','auto');
    print(gcf,'-dpng','-r500',ff);
    close(gcf)

end

%% and look at 'other results

otherNames = {'dens' 'dens' 'dens' 'rdenseb' 'rdenseb' 'rdenseb'} ;
comStructNames2 = {'wsbm' 'mod' 'yeo' 'wsbm' 'mod' 'yeo'};

for fig = 1:length(comStructNames2)
    
    plotStats = otherStatMap_resultsXValR2{fig} ;

    plotNames = {'Linear R-squared' 'Quadratic R-squared' 'Poisson R-squared'};

    % because we are reporting the x-val R2 and these are 0-100
    plotStats = plotStats ./ 100 ;

    %% viz it
    figure

    subp = tight_subplot(1,3,[.01 .03],[.1 .01],[.02 .02]) ;
    % make it full screen
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);

    commVec = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.(comStructNames2{fig}));

    for tsb_idx = 1:3

        axes(subp(tsb_idx))
        viz_statsInMat_wComm(plotStats(:,:,tsb_idx),...
            commVec,'Paired',10,...
            0,[0 0.3],'Greys',plotNames{tsb_idx},0.15) ;


        set(gca,'XColor',[ 0 0 0 0.0001 ])
        set(gca,'YColor',[ 0 0 0 0.0001 ])
        ax = gca ;
        % ax.Color = [ 1 1 0.8 ] ;
        axis square

    end

    hold off

    % save it
    fileName = strcat(comStructNames2{fig},'_',otherNames{fig},'_regressR2.png');
    ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
    %set(gcf,'paperpositionmode','auto');
    print(gcf,'-dpng','-r500',ff);
    close(gcf)

    
end

%% plot some regressions

regressNames = {'linear' 'quadratic' 'poisson'} ;
cmap = brewermap(11,'paired') ;

for commStrucIdx = 1:length(comStructNames)
%for commStrucIdx = 1:1

    coms = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.((comStructNames{commStrucIdx}))) ;
    coms = unique(coms) ;
    
    tmpStatMap = statMapAcrossMethod{commStrucIdx} ;
    tmpResultsMap = otherStatMap_resultsXValR2{commStrucIdx} ;

    [~,b] = sort(tmpResultsMap(:),'descend') ;
    [i,j,k] = ind2sub(size(tmpResultsMap),b) ;

    % make a subplot to plot all in
    subp = tight_subplot(2,4,[.01 .03],[.1 .01],[.02 .02]) ;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);
    
    for idx = 1:8

        axes(subp(idx))
        
        % divide by 100 because R^2 reported as 0-100
        xvalR2 = tmpResultsMap(i(idx),j(idx),k(idx)) ./ 100 ;

        % function [ pl ] = viz_StatMapEntryRegression(statEntry,vertexBool,colorVec)
        pl = viz_StatMapEntryRegression(tmpStatMap{i(idx),j(idx),k(idx)},0,[ 0 0 0 ]) ;
        pl.MarkerEdgeAlpha = 0.1 ;
        xlim([ (min(pl.XData)-2) (max(pl.XData)+2)]);

        axis square

        % add some text
        ytext = get(subp(idx),'ylim') ;
        ytext = ytext(2) * 0.9;

        color1 = num2str(cmap(coms(i(idx)),:)) ;
        color2 = num2str(cmap(coms(j(idx)),:)) ;
        
        annotText = { strcat(regressNames{k(idx)},' regression') ... 
            strcat('block interaction:',32, ...
            strcat('\color[rgb]{',color1,'}', num2str(coms(i(idx)))   ),...
            32 , '\color{black}-' , 32 , ...
            strcat('\color[rgb]{',color2,'}', num2str(coms(j(idx)))   ) ),...
            strcat('\color{black}R^2:', 32 , num2str(round(xvalR2, 3)) ) } ;
            
        text((min(pl.XData)+5),ytext, annotText,'FontSize',10,'VerticalAlignment','cap')
    end

    % save it
    fileName = strcat(comStructNames{commStrucIdx},'_interRegress.png');
    ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
    %set(gcf,'paperpositionmode','auto');
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
    
end

