
% clc
% clearvars
% 
% load('data/processed/yeo_both_normalpoisson_a0p5_fit_wsbm_script_v7p3.mat')
% load('data/interim/yeo_both_normalpoisson_a0p5_comVecs.mat')

clc
clearvars

%% load the necessary data

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_templateModel_1.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
load(loadName) ;

%% look at cost block mat

nSubj = length(dataStruct);
nNodes = templateModel.Data.n ;
nComm = templateModel.R_Struct.k ; 

nBlockInteract = (nComm^2 + nComm) / 2 ;
nBlockInteract_yeo = (49 + 7) / 2 ;

getIdx = ~~triu(ones(nComm));
getIdxYeo = ~~triu(ones(7));

% get the template data
templateData = templateModel.Data.Raw_Data ;
templateData(isnan(templateData)) = 0;

templateSubj_data = dataStruct(datasetDemo.age > 25 & datasetDemo.age <= 35) ;
[~,templateModelLens,~,templateModelAvgCount] = make_template_mat(templateSubj_data, ...
    2:58, ...
    60:116, ...
    1) ; 

%need to resize the lens
templateModelLens = templateModelLens(selectNodesFrmRaw,selectNodesFrmRaw);
% clear diagonal
templateModelLens(1:nNodes+1:end)=0; 

% get the model cost mat!
templateModel_costMat = templateModelLens .* templateModelAvgCount ;

[~,wsbm_cost_avgBlock] = get_block_mat(templateModel_costMat,comVecs.wsbm);
[~,mod_cost_avgBlock] = get_block_mat(templateModel_costMat,comVecs.mod);
[~,yeo_cost_avgBlock] = get_block_mat(templateModel_costMat,comVecs.yeo);

%% across subj

wsbm_avgCostMat = zeros([nBlockInteract nSubj]);
mod_avgCostMat = zeros([nBlockInteract nSubj]);
yeo_avgCostMat = zeros([nBlockInteract_yeo 7 nSubj]);

subjLensMat = zeros([ nNodes nNodes nSubj ]);
subjDistMat = zeros([ nNodes nNodes nSubj ]);

for idx = 1:nSubj

    % mask out AdjMat entries below mask_thr
    tmpMask = dataStruct(idx).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > MASK_THR_INIT ;    
    tmpMask(tmpMask > 0) = 1 ;   
    
    % get streamline lengths
    tmpLens = dataStruct(idx).lensMat(selectNodesFrmRaw,selectNodesFrmRaw);
    %tmpAdj(isnan(tmpAdj)) = 0 ;
    subjLensMat(:,:,idx) = tmpLens .* tmpMask ;
    
    % get eud distances
    subjDistMat(:,:,idx) = dataStruct(idx).distCoorMM(selectNodesFrmRaw,selectNodesFrmRaw) ;
    
    tmpCount = dataStruct(idx).countMat(selectNodesFrmRaw, selectNodesFrmRaw) ;    
    
    % get cost
    tmpCost = (tmpMask .* tmpCount) .* tmpLens;  
    
%     volMat(idx) = sum(tmpCost(:)) ;
    
    [~,tmpAvg] = get_block_mat(tmpCost,comVecs.wsbm) ; 
    wsbm_avgCostMat(:,idx) = tmpAvg(getIdx);
    
    [~,tmpAvg] = get_block_mat(tmpCost,comVecs.mod) ; 
    mod_avgCostMat(:,idx) = tmpAvg(getIdx);
    
    [~,tmpAvg] = get_block_mat(tmpCost,comVecs.yeo) ; 
    yeo_avgCostMat(:,idx) = tmpAvg(getIdxYeo);
    
    % and distances
    
    wsbm_costVec_cb(idx) = pdist([ wsbm_avgCostMat(:,idx)'; wsbm_cost_avgBlock(getIdx)' ] ,...
        'cityblock');
    mod_costVec_cb(idx) = pdist([ mod_avgCostMat(:,idx)'; mod_cost_avgBlock(getIdx)' ] ,...
        'cityblock');
    yeo_costVec_cb(idx) = pdist([ yeo_avgCostMat(:,idx)' ; yeo_cost_avgBlock(getIdxYeo)' ] ,...
        'cityblock');
    
    wsbm_costVec_cos(idx) = 1 - pdist([ wsbm_avgCostMat(:,idx)'; wsbm_cost_avgBlock(getIdx)' ] ,...
        'cos');
    mod_costVec_cos(idx) = 1 - pdist([ mod_avgCostMat(:,idx)'; mod_cost_avgBlock(getIdx)' ] ,...
        'cos');
    yeo_costVec_cos(idx) = 1 - pdist([ yeo_avgCostMat(:,idx)'; yeo_cost_avgBlock(getIdxYeo)' ] ,...
        'cos'); 
    
end

%% lets just measure the distance of the communities

wsbm_com_totLen = zeros([ nComm nSubj ]);
mod_com_totLen = zeros([ nComm nSubj ]);

wsbm_com_totDist = zeros([ nComm nSubj ]);
mod_com_totDist = zeros([ nComm nSubj ]);

for idx = 1:nSubj
    
    % get the length matrix
    tmpAdj = subjLensMat(:,:,idx) ;
    
    wsbm_tmpBlMat = get_block_mat(tmpAdj,comVecs.wsbm) ;
    mod_tmpBlMat = get_block_mat(tmpAdj,comVecs.mod);
    wsbm_com_totLen(:,idx) = sum(wsbm_tmpBlMat(~~eye(nComm)),2) ;
    mod_com_totLen(:,idx) = sum(mod_tmpBlMat(~~eye(nComm)),2) ; 
    
    % get the dist matrix
    tmpAdj = subjDistMat(:,:,idx) ;
    
    wsbm_tmpBlMat = get_block_mat(tmpAdj,comVecs.wsbm) ;
    mod_tmpBlMat = get_block_mat(tmpAdj,comVecs.mod);
    wsbm_com_totDist(:,idx) = sum(wsbm_tmpBlMat(~~eye(nComm)),2) ;
    mod_com_totDist(:,idx) = sum(mod_tmpBlMat(~~eye(nComm)),2) ; 
    
end

wsbm_subj_len = mean(wsbm_com_totLen,1) ;
mod_subj_len = mean(mod_com_totDist,1);

wsbm_subj_dist = mean(wsbm_com_totDist,1) ;
mod_subj_dist = mean(mod_com_totDist,1) ;


% wsbm community structure necessitates much more streamline length

% IPN_icc(wsbm_com_totLen,3,'single')
% IPN_icc(mod_com_totLen,3,'single')










