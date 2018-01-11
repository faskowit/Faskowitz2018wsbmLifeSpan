%% ANALYSIS LOOKING AT MODEL PREDICTION

clc
clearvars

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

outIntermPrefix = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR);
outProcessPrefix = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR);

load(strcat(outProcessPrefix,'_fit_wsbm_script_v7p3.mat'))
load(strcat(outProcessPrefix,'_fit_mod_v7p3.mat'))
load(strcat(outIntermPrefix,'_comVecs.mat'))
%load('data/interim/subj_dataStruct.mat')
load(strcat(outIntermPrefix,'_templateModel_1.mat'))
load(strcat(outProcessPrefix,'_basicData_v7p3.mat'))

%% lets look at how some metrics look across lifespan, using changing ca
% definition of paritions

nSubj = length(dataStruct) ;
nBlocks = templateModel.R_Struct.k ;
nNodes = templateModel.Data.n ;

nBlockInteract = (nBlocks^2 + nBlocks) / 2 ;

getIdx = ~~triu(ones(nBlocks));

subjDataMat = zeros([ nNodes nNodes nSubj ]);
subjWsbmCA = zeros([ nNodes nSubj ]);
subjModCA = zeros([ nNodes nSubj ]);

%% unroll block vec vars

blockVec_wsbm = zeros([nBlockInteract nSubj]);
blockVec_mod = zeros([nBlockInteract nSubj]) ;

avgBlockVec_wsbm = zeros([nBlockInteract nSubj]);
avgBlockVec_mod = zeros([nBlockInteract nSubj]) ;

binBlockVec_wsbm = zeros([nBlockInteract nSubj]);
binBlockVec_mod = zeros([nBlockInteract nSubj]) ;

%% also community convec
com_avgBlockVec_wsbm = zeros([nBlocks nBlocks nSubj]);
com_avgBlockVec_mod = zeros([nBlocks nBlocks nSubj]) ;

%% gather some of the data
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
    
    % get the fit community assignments 
    [~,tmpCA] = community_assign(fitWSBMAllStruct(idx).centModel) ;
    subjWsbmCA(:,idx) = CBIG_HungarianClusterMatch(comVecs.wsbm,tmpCA);
    subjModCA(:,idx) = CBIG_HungarianClusterMatch(comVecs.mod,...
        fitModAllStruct(idx).caMod_subj);
    
end

%% some runs of modularity did not render paritions with correct num coms
% and therefore have community assignment of all 0's
% we need to record these to skip over them

modExclude = (sum(subjModCA) == 0) ;

%% quick analysis... lets look at variability of nodal assignment

wsbmCnsns = get_nodal_versatility(subjWsbmCA) ;
modCnsns = get_nodal_versatility(subjModCA(:,~modExclude));

%% iterate
for idx = 1:nSubj
       
    % wsbm
    [tmpBl,avgtmpbl,~,bintmpbl] = get_block_mat(tmpAdj,subjWsbmCA(:,idx));
    blockVec_wsbm(:,idx) = tmpBl(getIdx);
    avgBlockVec_wsbm(:,idx) = avgtmpbl(getIdx);
    binBlockVec_wsbm(:,idx) = bintmpbl(getIdx);
   
    com_avgBlockVec_wsbm(:,:,idx) = avgtmpbl ;
    
    if modExclude(idx) == 0
    
        % mod
        [tmpBl,avgtmpbl,~,bintmpbl] = get_block_mat(tmpAdj,subjModCA(:,idx));
        blockVec_mod(:,idx) = tmpBl(getIdx);
        avgBlockVec_mod(:,idx) = avgtmpbl(getIdx);
        binBlockVec_mod(:,idx) = bintmpbl(getIdx);
        com_avgBlockVec_mod(:,:,idx) = avgtmpbl;  
        
    else
        
        % set zeros
        blockVec_mod(:,idx) = zeros([nBlockInteract 1]);
        avgBlockVec_mod(:,idx) = zeros([nBlockInteract 1]);
        binBlockVec_mod(:,idx) = zeros([nBlockInteract 1]);
        com_avgBlockVec_mod(:,:,idx) = zeros(nBlocks);  
          
    end
    
end
  
%% correlate each block pattern with model

% now different is predicted from actual?
templateData = templateModel.Data.Raw_Data;
templateData(isnan(templateData)) = 0 ;

% model_edgeVec = templateModel.Para.predict_e ;
% model_weiVec = templateModel.Para.predict_w ;
wsbmPred_weiVec =  templateModel.Para.predict_e .*  templateModel.Para.predict_w; 
wsbmPred_weiVec(isnan(wsbmPred_weiVec)) = 0;
[~,wsbm_weiVec_sqr] = make_square(wsbmPred_weiVec);

% and now creat the 'model' for the modular parition
[~,mod_weiVec_empir] = get_block_mat(templateData,comVecs.mod);
mod_weiVec_empir = mod_weiVec_empir(getIdx);
[~,mod_weiVec_empir_sqr] = make_square(mod_weiVec_empir);

%% community pattern vars

wsbm_comms_weiVec_corr = zeros([nBlocks nSubj]) ;
mod_comms_weiVec_corr = zeros([nBlocks nSubj]) ;

wsbm_comms_weiVec_cb = zeros([nBlocks nSubj]) ;
mod_comms_weiVec_cb = zeros([nBlocks nSubj]) ;

wsbm_comms_weiVec_cos = zeros([nBlocks nSubj]) ;
mod_comms_weiVec_cos = zeros([nBlocks nSubj]) ;

%% iterate
for idx = 1:nSubj
   
    % correlation 
    corrStr='pearson';
    wsbm_weiVec_corr(idx) = corr(wsbmPred_weiVec,avgBlockVec_wsbm(:,idx),...
        'type',corrStr) ;
    mod_weiVec_corr(idx) = corr(mod_weiVec_empir,avgBlockVec_mod(:,idx),...
        'type',corrStr) ;
    
    % other distances   
    wsbm_weiVec_cb(idx) = pdist([wsbmPred_weiVec' ; avgBlockVec_wsbm(:,idx)'],'cityblock') ;
    mod_weiVec_cb(idx) = pdist([mod_weiVec_empir' ; avgBlockVec_mod(:,idx)'],'cityblock') ;

    wsbm_weiVec_eud(idx) = pdist([wsbmPred_weiVec' ; avgBlockVec_wsbm(:,idx)'],'euclidean') ;
    mod_weiVec_eud(idx) = pdist([mod_weiVec_empir' ; avgBlockVec_mod(:,idx)'],'euclidean') ;
    
    wsbm_weiVec_cos(idx) = 1 - pdist([wsbmPred_weiVec' ; avgBlockVec_wsbm(:,idx)'],'cosine') ;
    mod_weiVec_cos(idx) = 1 - pdist([mod_weiVec_empir' ; avgBlockVec_mod(:,idx)'],'cosine') ;
    
    % distances of each community
    for jdx=1:nBlocks
        
        % correlation 
        wsbm_comms_weiVec_corr(jdx,idx) = corr(wsbm_weiVec_sqr(:,jdx),...
            com_avgBlockVec_wsbm(:,jdx,idx),...
            'type',corrStr) ;
        mod_comms_weiVec_corr(jdx,idx) = corr(mod_weiVec_empir_sqr(:,jdx),...
            com_avgBlockVec_mod(:,jdx,idx),...
            'type',corrStr) ;
        
        % other distances
        wsbm_comms_weiVec_cb(jdx,idx) = pdist([wsbm_weiVec_sqr(:,jdx)' ;...
            com_avgBlockVec_wsbm(:,jdx,idx)' ],...
            'cityblock') ;
        mod_comms_weiVec_cb(jdx,idx) = pdist([mod_weiVec_empir_sqr(:,jdx)' ; ...
            com_avgBlockVec_mod(:,jdx,idx)'],...
            'cityblock') ;
 
        wsbm_comms_weiVec_cos(jdx,idx) = 1 - pdist([wsbm_weiVec_sqr(:,jdx)' ;...
            com_avgBlockVec_wsbm(:,jdx,idx)' ],...
            'cosine') ;
        mod_comms_weiVec_cos(jdx,idx) = 1 - pdist([mod_weiVec_empir_sqr(:,jdx)' ; ...
            com_avgBlockVec_mod(:,jdx,idx)'],...
            'cosine') ; 
    end
    
end

%% run the regression
funcArgs = {1 500 [] 0 } ;
fits = {'linear' 'quadratic' 'poisson'} ;

xVec = datasetDemo.age;

wsbm_regResult_cos = cell([length(fits) 1]) ;
mod_regResult_cos = cell([length(fits) 1]) ;

for idx = 1:length(fits)

    Y = wsbm_weiVec_cos';
    res = struct() ; 
    [ res.xvalR2 , ...
        res.xvalsqErr, ... 
        res.xvalYhat, ...
        res.coef, ...
        res.lsFitStruct,...
        res.permStruct ] ...
        = nc_FitAndEvaluateModels(Y,xVec,fits{idx},funcArgs{:}) ;
    wsbm_regResult_cos{idx} = res ;
    
    Y =  mod_weiVec_cos' ;
    res = struct() ; 
    [ res.xvalR2 , ...
        res.xvalsqErr, ... 
        res.xvalYhat, ...
        res.coef, ...
        res.lsFitStruct,...
        res.permStruct ] ...
        = nc_FitAndEvaluateModels(Y,xVec,fits{idx},funcArgs{:}) ;
    mod_regResult_cos{idx} = res ;
    
end

%% figure out the fits that are best

[~,wsbm_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), wsbm_regResult_cos));
[~,mod_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), mod_regResult_cos));

wsbm_trend = wsbm_regResult_cos{wsbm_minIdx};
mod_trend = mod_regResult_cos{mod_minIdx};

%% plot it

fig1 = figure ;
pl = viz_blockRegress(wsbm_trend,0, [0    0.4470    0.7410] ) ;
set(pl,'MarkerFaceColor',[0    0.4470    0.7410  ],...
    'MarkerFaceColor',[0    0.4470    0.7410  ],...
    'MarkerFaceAlpha',.1,...
    'MarkerEdgeAlpha',.2)
hold
pl = viz_blockRegress(mod_trend,0,[0.8500    0.3250    0.0980]);
set(pl,'MarkerFaceColor', [0.8500    0.3250    0.0980],...
    'MarkerEdgeColor',[0.8500    0.3250    0.0980],...
    'MarkerFaceAlpha',.1,...
    'MarkerEdgeAlpha',.2) 

%% for the city block too

funcArgs = {1 500 [] 0 } ;
fits = {'linear' 'quadratic' 'poisson'} ;

xVec = datasetDemo.age;

wsbm_regResult_cb = cell([length(fits) 1]) ;
mod_regResult_cb = cell([length(fits) 1]) ;

for idx = 1:length(fits)

    Y = wsbm_weiVec_cb';
    res = struct() ; 
    [ res.xvalR2 , ...
        res.xvalsqErr, ... 
        res.xvalYhat, ...
        res.coef, ...
        res.lsFitStruct,...
        res.permStruct ] ...
        = nc_FitAndEvaluateModels(Y,xVec,fits{idx},funcArgs{:}) ;
    wsbm_regResult_cb{idx} = res ;
    
    Y =  mod_weiVec_cb' ;
    res = struct() ; 
    [ res.xvalR2 , ...
        res.xvalsqErr, ... 
        res.xvalYhat, ...
        res.coef, ...
        res.lsFitStruct,...
        res.permStruct ] ...
        = nc_FitAndEvaluateModels(Y,xVec,fits{idx},funcArgs{:}) ;
    mod_regResult_cb{idx} = res ;
    
end

[~,wsbm_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), wsbm_regResult_cb));
[~,mod_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), mod_regResult_cb));

%%

wsbm_trend_cb = wsbm_regResult_cb{wsbm_minIdx};
mod_trend_cb = mod_regResult_cb{mod_minIdx};

fig1 = figure ;
pl = viz_blockRegress(wsbm_trend_cb,0.12, [0    0.4470    0.7410] ) ;
set(pl,'MarkerFaceColor',[0    0.4470    0.7410  ],...
    'MarkerFaceColor',[0    0.4470    0.7410  ],...
    'MarkerFaceAlpha',.1,...
    'MarkerEdgeAlpha',.2)
hold
pl = viz_blockRegress(mod_trend_cb,0.15,[0.8500    0.3250    0.0980]);
set(pl,'MarkerFaceColor', [0.8500    0.3250    0.0980],...
    'MarkerEdgeColor',[0.8500    0.3250    0.0980],...
    'MarkerFaceAlpha',.1,...
    'MarkerEdgeAlpha',.2) 

%% extra viz

for idx=1:10
   subplot(2,5,idx) 
   plot(fitlm(datasetDemo.age,wsbm_comms_weiVec_cos(idx,:),'quadratic'));
   ylim([0 0.5])
end
suptitle('wsbm comm cos')

figure
for idx=1:10
   subplot(2,5,idx) 
   plot(fitlm(datasetDemo.age,mod_comms_weiVec_cos(idx,:),'quadratic'));
   ylim([0 0.5])
end
suptitle('mod comm cos')

% figure;
% for idx=1:10
%    subplot(2,5,idx) 
%    plot(fitlm(datasetDemo.age,wsbm_comms_weiVec_eud(idx,:),'quadratic'));
%    
% end



for idx=1:10
   subplot(2,5,idx) 
   plot(fitlm(datasetDemo.age,wsbm_comms_weiVec_cb(idx,:),'quadratic'));
   ylim([0 0.5])
end
suptitle('wsbm comm eud')

figure
for idx=1:10
   subplot(2,5,idx) 
   plot(fitlm(datasetDemo.age,mod_comms_weiVec_cb(idx,:),'quadratic'));
   ylim([0 0.5])
end
suptitle('mod comm eud')

