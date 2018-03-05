%% lets analyze the age bin scripts
clc
clearvars

% addpath('~/JOSHSTUFF/scripts/BCT/2017_01_15_BCT/')

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

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
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

%% vi confusion

vi_confusion = zeros(age_bins);
mi_confusion = zeros(age_bins);

for idx = 1:age_bins 
   for jdx = 1:age_bins 
    
        [vi_confusion(idx,jdx),mi_confusion(idx,jdx)] = ...
            partition_distance(agebin_ca(:,idx),agebin_ca(:,jdx));
   end
end

% DOESNT WORK
% %% lets fit the young model to old...
% % and the old model to the young
% 
% oldestModel = consensusCrossIdx{2} ;
% 
% % get the harsh_mu, with binary probabilities
% [ ~ , old_harsh_mu ] = make_WSBM_prior(oldestModel , 1) ;
% % make NULL model input
% % (need to have at least one trial...but it wont do anything) 
% old_nullModelStruct = indivModelInputStruct ;
% old_nullModelStruct.mu_0 = old_harsh_mu ;
% old_nullModelStruct.numTrials = 1 ;
% old_nullModelStruct.muMaxIter = 0 ;
% old_nullModelStruct.mainMaxIter = 0 ;
% old_nullModelStruct.verbosity = 0 ;
% % make this struct into a cell list that the wsbm likes
% a = struct2nv(old_nullModelStruct) ;
% b = struct2cell(old_nullModelStruct) ;
% c = [ a b ]';
% old_nullModelInput = c(:)' ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% youngModel = consensusCrossIdx{1} ;
% 
% % get the harsh_mu, with binary probabilities
% [ ~ , young_harsh_mu ] = make_WSBM_prior(youngModel , 1) ;
% % make NULL model input
% % (need to have at least one trial...but it wont do anything) 
% young_nullModelStruct = indivModelInputStruct ;
% young_nullModelStruct.mu_0 = young_harsh_mu ;
% young_nullModelStruct.numTrials = 1 ;
% young_nullModelStruct.muMaxIter = 0 ;
% young_nullModelStruct.mainMaxIter = 0 ;
% young_nullModelStruct.verbosity = 0 ;
% % make this struct into a cell list that the wsbm likes
% a = struct2nv(young_nullModelStruct) ;
% b = struct2cell(young_nullModelStruct) ;
% c = [ a b ]';
% young_nullModelStruct = c(:)' ;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% youngData = fitWSBMAllStruct(~~templateIdMat(:,1));
% youngDataLogEvid = zeros([ length(youngData) 1 ]);
% youngDataFitOld = cell([length(youngData) 1]);
% youngDataFitYoung = cell([length(youngData) 1]);
% 
% for idx = 1:length(youngData)
%     
%     [~,youngDataFitOld{idx}] = wsbm(youngData(idx).Raw_Data,...
%         sym_RStruct(11),...
%         old_nullModelInput{:}) ;
%      
%     [~,youngDataFitYoung{idx}] = wsbm(youngData(idx).Raw_Data,...
%         youngData(idx).Model(1).R_Struct.R,...
%         young_nullModelStruct{:}) ;
%     
% end

%% subj data

nSubj = length(dataStruct) ;

subjDataMat = zeros([nNodes nNodes nSubj]) ;

for idx = 1:nSubj
    
    tmpAdj = dataStruct(idx).countVolNormMat(selectNodesFrmRaw, selectNodesFrmRaw);
    % get rid of the diagonal
    %n=size(tmpAdj,1);
    tmpAdj(1:nNodes+1:end) = 0; 
    % mask out AdjMat entries below mask_thr
    tmpAdj_mask = dataStruct(idx).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > 1 ;    
    tmpAdj_mask(tmpAdj_mask > 0) = 1 ;   
    tmpAdj = tmpAdj .* tmpAdj_mask ;

    subjDataMat(:,:,idx) = tmpAdj ;

end

%% create the regressions at each agebin

wsbm_avgBlock_cell = cell([age_bins 1]) ;
wsbm_vec_dist_cos_cell = cell([age_bins 1]) ;
wsbm_vec_dist_cb_cell = cell([age_bins 1]) ;

for idx = 1:age_bins
    
    %ageBinData = fitWSBMAllStruct(~~templateIdMat(:,idx));
    ageBinTemplate = consensusCrossIdx{idx} ;
    [~,abt_comVec] = community_assign(ageBinTemplate) ;
   
    nBlocks = ageBinTemplate.R_Struct.k ;
    nNodes = ageBinTemplate.Data.n ;
    getIdx = ~~triu(ones(nBlocks));
    nBlockInteract = sum(sum(getIdx)) ;

    wsbm_weiVec_predict =  ageBinTemplate.Para.predict_e .*  ageBinTemplate.Para.predict_w; 
    wsbm_weiVec_predict(isnan(wsbm_weiVec_predict)) = 0;
    [~,wsbm_weiVec_predict_sqr] = make_square(wsbm_weiVec_predict);
    
    wsbm_avgBlock_cell{idx} = zeros([nBlockInteract nSubj]);

    wsbm_vec_dist_cos_cell{idx} = zeros([nSubj 1]) ;
    wsbm_vec_dist_cb_cell{idx} = zeros([nSubj 1]);
    
    %% iterate
    for sub = 1:nSubj

        tmpAdj = subjDataMat(:,:,sub) ;

        % wsbm
        [~,avgtmpbl] = get_block_mat(tmpAdj,abt_comVec);
        wsbm_avgBlock_cell{idx}(:,sub) = avgtmpbl(getIdx);

        wsbm_vec_dist_cb_cell{idx}(sub) = pdist([wsbm_weiVec_predict' ; wsbm_avgBlock_cell{idx}(:,sub)'],'cityblock') ;
        wsbm_vec_dist_cos_cell{idx}(sub) = 1 - pdist([wsbm_weiVec_predict' ; wsbm_avgBlock_cell{idx}(:,sub)'],'cosine') ;
        
    end

end

funcArgs = {1 500 [] 0 } ;
%funcArgs = {1 2 [] 0 } ;

fits = {'linear' 'quadratic' 'poisson'} ;

xVec = datasetDemo.age;

ageBin_regress_cos = cell([age_bins 1]) ;
ageBin_regress_cb = cell([age_bins 1]) ;

for idx = 1:age_bins
        
    % setup results for wsbm, mod, yeo....
    ageBin_regress_cos{idx} = cell([ length(fits) 1]);
    ageBin_regress_cb{idx} = cell([ length(fits) 1]);
        
    for jdx = 1:length(fits)
    
        disp(jdx)
        
        Y = wsbm_vec_dist_cos_cell{idx} ;
        
        regressResults = struct() ; 
        [ regressResults.xvalR2 , ...
            regressResults.xvalsqErr, ... 
            regressResults.xvalYhat, ...
            regressResults.coef, ...
            regressResults.lsFitStruct,...
            regressResults.permStruct ] ...
            = nc_FitAndEvaluateModels(Y,xVec,fits{jdx},funcArgs{:}) ;
        ageBin_regress_cos{idx}{jdx} = regressResults ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Y = wsbm_vec_dist_cb_cell{idx} ;
        
        regressResults = struct() ; 
        [ regressResults.xvalR2 , ...
            regressResults.xvalsqErr, ... 
            regressResults.xvalYhat, ...
            regressResults.coef, ...
            regressResults.lsFitStruct,...
            regressResults.permStruct ] ...
            = nc_FitAndEvaluateModels(Y,xVec,fits{jdx},funcArgs{:}) ;
        ageBin_regress_cb{idx}{jdx} = regressResults ;
        
    end
end

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_ageBin_regResults.mat');
save(outName,...
    'ageBin_regress*',...
    ...
    '-v7.3')

%% VIEW IT

figure


% subp = tight_subplot(1,5,[.10 .05],[.1 .05],[.1 .05]) ;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.95 0.55]);      

cmap = brewermap(5,'purples') ;

for idx = 1:length(ageBin_regress_cb)
    
    disp(idx)
    
    currentRegResults = ageBin_regress_cb{idx} ;

    [~,wsbm_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults));

    wsbm_trend = currentRegResults{wsbm_minIdx};
    %fig1 = figure ;
    
%     minus for because we indexing 3-6
%     axes(subp(idx / 3)) 

    % wsbm
    [ pl , ln3 ] = viz_blockRegress(wsbm_trend,0, cmap(idx,:)) ;
    set(pl,'MarkerFaceColor',cmap(idx,:),...
        'MarkerFaceColor',cmap(idx,:),...
        'MarkerFaceAlpha',.15,...
        'MarkerEdgeAlpha',.1)
   hold
    
%     xlim([min(pl.XData)-3 max(pl.XData)+3])
%     
%     ylabel(reg_ylabel_names{idx})
%     xlabel('Age')
%     
%     title(reg_result_names{idx})
    
%     yrange = ylim ;
%     yrangeAbs = yrange(2) - yrange(1) ;
%     % add the R2!! 
%     if idx < 4
%         ypos = yrangeAbs * 0.15;
%         ypos = yrange(1) + ypos ;
%         ll = legend([ln3 ln2 ln1 ], {'WSBM' 'Modular' 'Yeo'},'Location','SouthEast','FontSize',12) ;
%         set(ll,'Units','inches')
%                 legend('boxoff')
%     else
%         ypos = yrangeAbs * 0.97;
%         ypos = yrange(1) + ypos ;
%         ll = legend([ln3 ln2 ln1 ], {'WSBM' 'Modular' 'Yeo'},'Location','NorthEast','FontSize',12);
%         set(ll,'Units','inches')
%         legend('boxoff')
%     end

%     annotText = { strcat('WSBM R^2:',32,num2str(round(wsbm_trend.xvalR2 / 100,3))) ...
%         strcat('Modular R^2:',32,num2str(round(mod_trend.xvalR2 / 100,3))) ...
%         strcat('Yeo R^2:',32,num2str(round(yeo_trend.xvalR2 / 100,3))) ...
%         } ;
% 
%     text((min(pl.XData)),ypos, annotText,'FontSize',12,'VerticalAlignment','cap')    
        
end














