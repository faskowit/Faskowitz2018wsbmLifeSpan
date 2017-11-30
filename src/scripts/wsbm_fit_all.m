%% load a config file that sets global parameters

% need to edit config file string to match what you want!!
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
config_file='config_yeo_take3.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

%% assign couple more vars baed on config

if strcmp(HEMI_CHOICE,'both') 
    selectNodesFrmRaw = [LEFT_HEMI_NODES RIGHT_HEMI_NODES];
    selectNodesFrmAvgBH = 1:size(selectNodesFrmRaw) ; 
elseif strcmp(HEMI_CHOICE,'left')
    selectNodesFrmRaw = LEFT_HEMI_NODES ;
    selectNodesFrmAvgBH = 1:size(selectNodesFrmRaw) ;
elseif strcmp(HEMI_CHOICE,'right')
    
    selectNodesFrmRaw = RIGHT_HEMI_NODES ;
    selectNodesFrmAvgBH = ...
        (size(selectNodesFrmRaw)+1):size(avgTempBothHemi,1) ;   
else
    disp('invalide hemi choice') 
end

%% load the NKI data

read_data_mat=strcat(PROJECT_DIR,'/',DATA_DIR,'/',PARCELLATION,'datasetStruct.mat');

if(exist(read_data_mat,'file') == 2)
    disp('already packaged data, will read in mat file')
    readData = load(read_data_mat);
    readData = readData.datasetStruct ;
else
    readData = read_nki_data( strcat(PROJECT_DIR,DATA_DIR) , PARCELLATION) ; 
end

dataStruct = readData.dataRaw ;
datasetDemo = readData.demoRaw ;

%% filter sub-optimal data here
% condition the data based on exclusion criteria

nSubj = length(dataStruct) ;
sparse_cutoff = 0.25 ;
%sparse_cutoff = 0.15 ;

sparseness = zeros([length(dataStruct) 1]) ;

removeVec = zeros( [nSubj 1] );

for idx=1:nSubj
   
    tmpSubjMat = ...
    dataStruct(idx).countVolNormMat(selectNodesFrmRaw, ...
    selectNodesFrmRaw ) ;
    
    sparseness(idx) = density_und(tmpSubjMat);
    %disp(sparseness)
    
    if sparseness(idx) < sparse_cutoff
       
        removeVec(idx) = 1 ;
    end
end

% take out bad subjects
% by filters by subjects not marked for removal
dataStruct = dataStruct(removeVec == 0);
datasetDemo = datasetDemo(removeVec == 0, :);

%% pick out young subs

youngIdx = datasetDemo.age > ageLowLim & datasetDemo.age <= ageHighLim ;
templateSubj_data = dataStruct(youngIdx == 1) ;
clear youngIdx

%% setup a template to run wsbm on

avgTempBothHemi = make_template_mat(templateSubj_data, ...
    LEFT_HEMI_NODES, ...
    RIGHT_HEMI_NODES, ...
    MASK_THR_INIT) ; 

if strcmp(HEMI_CHOICE,'both')
    
    avgTemp = avgTempBothHemi .* 1;
    
elseif strcmp(HEMI_CHOICE,'left')
    
    avgTemp = avgTempBothHemi(selectNodesFrmAvgBH, ...
        selectNodesFrmAvgBH) .* 1;
    
elseif strcmp(HEMI_CHOICE,'right')
    
    avgTemp = avgTempBothHemi(selectNodesFrmAvgBH, ...
        selectNodesFrmAvgBH) .* 1 ;    
else
    disp('invalide hemi choice') 
end

% actually replace the 0's with NaN
avgTemp(avgTemp == 0) = NaN ;

%% Model Selection for the number of groups K for young subject group
% setup the modelInputs var

% Create a List of Model Inputs to test
 kIterOver = 7:1:11; %   
%kIterOver = 12:1:13;

%set up the cell struct to pass to the looper
modelInputs = cell(numel(kIterOver),1); 
    
for idx = 1:numel(kIterOver)
    
    R_STRUCT_TO_TEST = sym_RStruct(kIterOver(idx)) ;
    
%     modelInputs{idx} = { R_STRUCT_TO_TEST, ... 
%         'W_Distr', WEIGHT_DIST, ...
%         'E_Distr', EDGE_DIST, ...
%         'alpha', INIT_ALPHA, ...
%         'numTrials', LOOPER_NUM_TRIAL , ...
%         'mainMaxIter', LOOPER_MAIN_ITER , ...
%         'muMaxIter' , LOOPER_MU_ITER };
    modelInputs{idx} = { R_STRUCT_TO_TEST, ... 
        'W_Distr', WEIGHT_DIST, ...
        'E_Distr', EDGE_DIST, ...
        'alpha', INIT_ALPHA, ...
        'mainMaxIter', LOOPER_MAIN_ITER , ...
        'muMaxIter' , LOOPER_MU_ITER,...
        'mainTol',0.01, ...
        'muTol', 0.01 ,...
        'verbosity', 0};

end

numTrialPtrn = [ 250 100 100 100 100 100 100 100 100 100 100 ] ;
% this will be one shorter because there is no prior weight for first run
% just put a 0 at the end so it doesnt mess it up
prirPtrn = [ 1 1.5 2 2.5 3.0 3.5 4.0 4.5 5.0 5.5 0 ] ;

% prirPtrn = [ 1:0.25:4.5 0 ] ;
% numTrialPtrn = [ 250 repmat(100,1,15) ] ;

%% k looper part
% iterate over num communities, k, to find out which k give best evidence

[ kLooperResults , kLooperModels ] = wsbm_looper_wrapper(avgTemp, ...
    modelInputs, ...
    LOOPER_ITER, ...
    [], numTrialPtrn,...
    prirPtrn) ;

% save this, so that we can reconstruct kLooper results if needed, ~200Mb
save(strcat(OUTPUT_DIR , '/interim/', OUTPUT_STR, '_kLooperModels.mat' ), 'kLooperModels')

%% now figure out which k is best! 

% find max at each row
[ ~ , kLoopMaxIdx ] = max(kLooperResults(:,2:end)) ;
kLoopBestIdx = (mode(kLoopMaxIdx(2:end))) ;
kBest = kIterOver(kLoopBestIdx) ;

% get all the models with best k
kLoopBestModels = cell2struct(kLooperModels(kLoopBestIdx,:), ...
    'Model', LOOPER_ITER );

% and now lets find the centroid based on pairwise variation of info
[ kCentralModel , kSimVec , kSimVUpTri ] = central_model(kLoopBestModels) ; 

%% sort the kLooperModels...
% by aligning all to the 'most central' with Hungarian algo

% the ref to align to will be the kCentralModel
[~,ref] = community_assign(kCentralModel) ;

% preallocate the mat of the bestModels community assignments
kBestModels_ca_algn = zeros([ length(selectNodesFrmRaw) LOOPER_ITER ]) ;

% first stack all plausible parcellations
for idx=1:LOOPER_ITER
    
    [~,tmp] = ...
        community_assign(kLoopBestModels(idx).Model.Para.mu) ;

    kBestModels_ca_algn(:,idx) = CBIG_HungarianClusterMatch(ref,tmp) ;
    
end

%% get kSimVec and LogEvid vec sorted

% will be size: numKTested x LOOPER_ITER
kLooperLogEvids = cellfun(@(a) a.Para.LogEvidence,kLooperModels(:,:)) ;
% get the row with the best evids, and transpose to make col
kBestLogEvids = kLooperLogEvids(kLoopBestIdx,:)' ;

% sort evid vector and kSim vector
[ sort_logEvid, sort_logEvidIdx ] = sort(kBestLogEvids,'descend');
[ sort_kSimVec, sort_kSimVecIdx ] = sort(kSimVec,'descend');

% show that there is corrleation between logEvid and VI disatnce
logEvid_kSim_corr = corr(kBestLogEvids,kSimVec) ;

%% lets iterate over the 100 k best models to create consensus model
% first we have to turn the hungarian reordered ca_align mat into a wsbm
% prior to be input into wsbm inference fit

% lets cuttoff the worst 5% of the runs to clear up the prior a little...
% TODO... maybe no need for this...
% if no cutting of... then sorting doesn't matter-but Hungarian match still
% does matter
onlyKeep = 1; %0.95 ;

% first figure out the cuttoff
cutoff = floor(LOOPER_ITER * onlyKeep) ; 
kiter_prior = zeros([ kCentralModel.R_Struct.k kCentralModel.Data.n ]) ;

% sort the ca_align by similarity
kBestModels_ca_algn_sort = kBestModels_ca_algn(:,sort_kSimVecIdx);

% iterate over each node
for idx=1:(kCentralModel.Data.n)

    kiter_prior(:,idx) = ...
        sum(bsxfun(@eq,kBestModels_ca_algn_sort(idx,1:cutoff), ...
        [1:(kCentralModel.R_Struct.k)]'),2) ./ cutoff ;
        
end    

%% loop for consensus

CONSENSUS_ITER = 10 ;
C = zeros([CONSENSUS_ITER 1]) ;

for idx=1:CONSENSUS_ITER

    runNTimes = 100 ;
    
    rr = sym_RStruct(kBest) ;
    modIn = { ... 
        'W_Distr', WEIGHT_DIST, ...
        'E_Distr', EDGE_DIST, ...
        'alpha', INIT_ALPHA, ...
        'mainMaxIter', LOOPER_MAIN_ITER , ...
        'muMaxIter' , LOOPER_MU_ITER,...
        'mainTol',0.01, ...
        'muTol', 0.01 ,...
        'verbosity', 0, ...
        'numTrials', 50 ,...
        'mu_0', kiter_prior(:,:,idx)};
    
    % function [ allModels ] = wsbmFitNTimes( adjMat, rStruct , modelInputs , numFits , numCores)
    cnsnsusModels = wsbmFitNTimes(kCentralModel.Data.Raw_Data,...
        rr,...
        modIn,...
        runNTimes, 16) ;
        
    tmpCnsnsusCa = zeros([ kCentralModel.Data.n runNTimes ]) ;
    
    for jdx=1:runNTimes
        [~,tmpCnsnsusCa(:,jdx)] = community_assign(cnsnsusModels(jdx).Model) ;
    end
   
    tmpAgreeMat = agreement(tmpCnsnsusCa) ./ runNTimes ;
    
    % get the consensus consitency
    C(idx) = consensus_consistency(tmpAgreeMat) ;
    
    % make new kiter_prior for new loop
    for kdx=1:(kCentralModel.Data.n)

        kiter_prior(:,kdx,idx+1) = ...
            sum(bsxfun(@eq,tmpCnsnsusCa(kdx,:), ...
            [1:(kCentralModel.R_Struct.k)]'),2) ./ cutoff ;
        
    end    
    
    % have we converged? or are we at the end of loop?
    if C(idx) >= 0.95 || idx == CONSENSUS_ITER
        consensus_kiter_prior = kiter_prior(:,:,idx+1) ;
        [~,consensus_ca] = community_assign(consensus_kiter_prior) ; 
        consensus_kCentralModel = central_model(cnsnsusModels) ;
        % also add the data back into consensus_kCentral...
        consensus_kCentralModel.Data = kCentralModel.Data ;
        break 
    end
    
end

%% get centroid of best k and define it as templateModel

templateModel = consensus_kCentralModel ; 

plotMUnWSBM(templateModel)

% save the inital best model
save(strcat(OUTPUT_DIR , '/interim/', OUTPUT_STR, '_templateModel_1.mat'), 'templateModel')

%% lets make the prior assignment matrix

[ muPrior , ~ ] = make_WSBM_prior(templateModel , PRIOR_WEIGHT_MORE) ;

% save the inital mu
save(strcat(OUTPUT_DIR , '/interim/', OUTPUT_STR, '_muPrior_1.mat'), 'muPrior')

%% second loop for alpha level
% only if we are testing a single hemisphere, because then we can estimate
% alpha by comparing to other hemi

if ~strcmp(HEMI_CHOICE,'both')

    % save the current templateModel as old model
    templateModel_old = templateModel ;
    
    % get the other hemisphere
    if strcmp(HEMI_CHOICE,'left')
        other_hemi_select = (length(LEFT_HEMI_NODES) + 1):(length(avgTempBothHemi));
    else %right
        other_hemi_select = 1:(length(LEFT_HEMI_NODES));
    end

    avgTempOtherHemi = avgTempBothHemi(other_hemi_select,other_hemi_select) ;

    %actually replace the 0's with NaN
    avgTempOtherHemi(avgTempOtherHemi == 0) = NaN ;

    %% select alpha
        
    scoreVI = @(Model) -varInfo(templateModel, Model);
    alphaIterOver = 0.1:0.1:0.9 ;

    %reinitialize alpha_vals
    modelInputs = cell(numel(alphaIterOver),1); 

    for idx = 1:length(alphaIterOver),
        modelInputs{idx} = {sym_RStruct(kBest), ...
            'W_Distr', WEIGHT_DIST, ...
            'E_Distr', EDGE_DIST, ...
            'alpha', alphaIterOver(idx), ...
            'numTrials', LOOPER_NUM_TRIAL ...
            'verbosity', 0, ...
            'mainMaxIter', LOOPER_MAIN_ITER };
    end

    [ alphaLooperResults , alphaLooperModels ] = wsbm_looper_wrapper( ...
        avgTempOtherHemi, ...
        modelInputs, ...
        LOOPER_ITER, ...
        scoreVI ) ;
    
    % find max at each row
    [ ~ , kLoopMaxIdx ] = max(alphaLooperResults) ;
    alphaLoopBestIdx = (mode(kLoopMaxIdx(2:end))) ;
    alphaBest = alphaIterOver(alphaLoopBestIdx) ;

    %% need to learn best model one more time and again learn the community assignments 

    % get all the models with best k
    alphaLooperModels = cell2struct(alphaLooperModels(alphaLoopBestIdx,:), ...
        'Model', LOOPER_ITER );

    % and now lets find the median
    [ alplhaCentralModel , alphaSimVec ] = central_model(alphaLooperModels) ; 
  
    templateModel = alplhaCentralModel ;

else
 
    alphaBest = INIT_ALPHA ; 
    
end % if we doing the hemi comparison 

%% commnity assign again 

[ muPrior , harsh_mu ] = make_WSBM_prior(templateModel , PRIOR_WEIGHT_MORE) ;

save(strcat(OUTPUT_DIR , '/interim/', OUTPUT_STR, '_muPrior_2.mat'), 'muPrior')

%% get comunity for each of the verticies 
% this is just for display stuffs
% actually should do this again after the second Best_Model fit

communityLabels = community_assign(templateModel.Para.mu) ; 

%% all the subjects

disp('now working on fitting all subjects together')

% get number of subjects to run
datasetSize = length(dataStruct) ; 

% initialize struct
fitWSBMAllStruct = struct() ; 

% parameters for WSMB individual fit
k_to_use = templateModel.R_Struct.k  ;
alpha_to_use = alphaBest ;
w_dist_to_use = templateModel.W_Distr ;
e_dist_to_use = templateModel.E_Distr ;
r_struct_to_use = sym_RStruct(kBest) ;


%% with the iterative fits

% % lets not do this with the grad descent exactly right now....
% maxIters = 5 ; 
% %minIters = 5 ; 
% modelFits = 5 ; 
% %convergence = -0.00001 ; 
% 
% %for i = 1:length(raw_data)
% parallel_pool = gcp ; 
% parfor subj = 1:datasetSize
%     
%     disp('working on')
%     disp(subj)
%     
%     %% read in individual data
%         
%     % loop through all the data and fit the blockmodel with the prior
%     subjAdjMat = dataStruct(subj).countVolNormMat(selectNodesFrmRaw, selectNodesFrmRaw);
%     
%     % get rid of the diagonal
%     n=size(subjAdjMat,1);
%     subjAdjMat(1:n+1:end) = 0;
%     
%     % mask out AdjMat entries below mask_thr
%     subjAdjMat_mask = dataStruct(subj).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > MASK_THR ;    
%     subjAdjMat_mask(subjAdjMat_mask > 0) = 1 ;   
%     subjAdjMat = subjAdjMat .* subjAdjMat_mask ;
%     
%     % replace the 0's with NaN
%     subjAdjMat(subjAdjMat == 0) = NaN ;
% 
%     %% record stuff for the output struct
%     
%     % reord subj id in results struct
%     fitWSBMAllStruct(subj).id = dataStruct(subj).id ;
%     
%     % save the raw data only once
%     fitWSBMAllStruct(subj).Raw_Data = subjAdjMat ; 
%     
%     % initialize the variable prior
%     variableMuPrior = muPrior ;
%     
%     %% iterative fitting the wsbm
%     for idx=1:maxIters
%     
%         % keep this in the loop, might need to for variableMu
%         indivModelInputs = {'W_Distr', w_dist_to_use, ...
%                     'E_Distr', e_dist_to_use, ...
%                     'numTrials', INDIV_WSBM_NUM_TRIAL , ...
%                     'mu_0', variableMuPrior , ...
%                     'alpha',alpha_to_use , ...
%                     'mainMaxIter', INDIV_WSBM_MAIN_ITER , ...
%                     'muMaxIter', INDIV_WSBM_MU_ITER,  ...
%                     'verbosity' , 0 }; 
% 
%         indivCentModel = wsbmCentralFit( subjAdjMat, ...
%             (sym_RStruct(k_to_use)), ...
%             indivModelInputs, ...
%             modelFits, ...
%             variableMuPrior );
%     
%         [ variableMuPrior , ~ ] = make_WSBM_prior(indivCentModel, PRIOR_WEIGHT_MORE) ;
% 
%         % record this iter's central model
%         fitWSBMAllStruct(subj).Model(idx) = indivCentModel ;
% 
%     end
%    
% end % looping over each subject

%% perhaps changing this up...

modelFits = 5 ; 
maxIters = 1 ; 

[ muPrior , ~ ] = make_WSBM_prior(templateModel , 3) ;

indiv_numTrialPtrn = [ 250 100 100 100 100 100 100 100 100 100 100 ] ;
% this will be one shorter because there is no prior weight for first run
% just put a 0 at the end so it doesnt mess it up
% prirPtrn = [ 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 0 ] ;
indiv_prirPtrn = [ 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 0 ] ;

test_times = zeros([10 1]) ;

%for subj = 1:3
parallel_pool = gcp ; 
parfor subj = 1:datasetSize
    
    t1=tic ;

    disp('working on')
    disp(subj)
    
    %% read in individual data
        
    % loop through all the data and fit the blockmodel with the prior
    subjAdjMat = dataStruct(subj).countVolNormMat(selectNodesFrmRaw, selectNodesFrmRaw);
    
    % get rid of the diagonal
    n=size(subjAdjMat,1);
    subjAdjMat(1:n+1:end) = 0;
    
    % mask out AdjMat entries below mask_thr
    subjAdjMat_mask = dataStruct(subj).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > MASK_THR ;    
    subjAdjMat_mask(subjAdjMat_mask > 0) = 1 ;   
    subjAdjMat = subjAdjMat .* subjAdjMat_mask ;
    
    % replace the 0's with NaN
    subjAdjMat(subjAdjMat == 0) = NaN ;

    %% record stuff for the output struct
    
    % reord subj id in results struct
    fitWSBMAllStruct(subj).id = dataStruct(subj).id ;
    
    % save the raw data only once
    fitWSBMAllStruct(subj).Raw_Data = subjAdjMat ; 
    
    % initialize the variable prior
    %variableMuPrior = muPrior ;
    
    tempModelStruct = struct() ; 

    for idx=1:modelFits
         
        % keep this in the loop, might need to for variableMu
        indivModelInputs = {'W_Distr', w_dist_to_use, ...
                    'E_Distr', e_dist_to_use, ...
                    'alpha',alpha_to_use , ...
                    'mainMaxIter', INDIV_WSBM_MAIN_ITER , ...
                    'muMaxIter', INDIV_WSBM_MU_ITER,  ...
                    'verbosity' , 0 }; 
                
%         indivCentModel = wsbmCentralFit( subjAdjMat, ...
%             (sym_RStruct(k_to_use)), ...
%             indivModelInputs, ...
%             modelFits, ...
%             variableMuPrior );

%       [ variableMuPrior , ~ ] = make_WSBM_prior(indivCentModel, PRIOR_WEIGHT_MORE) ;

        %function [ Model ] = wsbmFitWPttrn( adjMat, rStruct , modelInputs , initMu, numTrialPttrn, priorWeightPttrn)
        tempModelStruct(idx).Model = wsbmFitWPttrn( subjAdjMat, r_struct_to_use , ...
            indivModelInputs , muPrior, indiv_numTrialPtrn, indiv_prirPtrn) ;
        
        fitWSBMAllStruct(subj).Model(idx) = tempModelStruct(idx).Model ;

    end
   
    %recrd the central model
    fitWSBMAllStruct(subj).centModel = central_model(tempModelStruct) ;
    
    test_times(subj) = toc(t1) ;
    
end % looping over each subject

% filter sub-optimal data here
% condition the data based on exclusion criteria
% 
% nSubj = length(dataStruct) ;
% sparse_cutoff = 0.25 ;
% removeVec = zeros( [nSubj 1] );
% 
% for idx=1:nSubj
%    
%     tmpSubjMat = ...
%     dataStruct(idx).countVolNormMat(selectNodesFrmRaw, ...
%     selectNodesFrmRaw ) ;
%     
%     sparseness = density_und(tmpSubjMat);
%     
%     if sparseness < sparse_cutoff
%        
%         removeVec(idx) = 1 ;
%     end
% end
% 
% % take out bad subjects
% % by filters by subjects not marked for removal
% dataStruct = dataStruct(removeVec == 0);
% datasetDemo = datasetDemo(removeVec == 0, :);
% fitWSBMAllStruct = fitWSBMAllStruct(removeVec == 0) ;

% % lets look real quick
% 
% kBestLogEvids = zeros([ size(fitWSBMAllStruct,2) 1 ]) ;
% 
% for idx=1:(size(fitWSBMAllStruct,2))
%    
%     kBestLogEvids(idx) = fitWSBMAllStruct(idx).Model(5).Para.LogEvidence ;
% end

% end of fitting WSBM on brains
% now time to do analysis...

% save all the variables we need
indivModelInputStruct = struct();
indivModelInputStruct.('W_Distr') = WEIGHT_DIST ; 
indivModelInputStruct.('E_Distr') = EDGE_DIST ;
indivModelInputStruct.('numTrials') = INDIV_WSBM_NUM_TRIAL;
indivModelInputStruct.('mu_0') = muPrior ;
indivModelInputStruct.('alpha') = alpha_to_use ;
indivModelInputStruct.('mainMaxIter') = INDIV_WSBM_MAIN_ITER ;
indivModelInputStruct.('muMaxIter') = INDIV_WSBM_MU_ITER ;

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_fit_wsbm_script.mat');
save(outName,...
    'dataStruct',...
    'datasetDemo',...
    'muPrior',...
    'fitWSBMAllStruct',...
    'templateModel',...
    'indivModelInputStruct',...
    'PROJECT_DIR',...
    'DATA_DIR',...
    'selectNodesFrmRaw',...
    'selectNodesFrmAvgBH'...
    )

%% yoy

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_fit_wsbm_script_v7p3.mat');
save(outName,...
    'dataStruct',...
    'datasetDemo',...
    'muPrior',...
    'fitWSBMAllStruct',...
    'templateModel',...
    'indivModelInputStruct',...
    'PROJECT_DIR',...
    'DATA_DIR',...
    'selectNodesFrmRaw',...
    'selectNodesFrmAvgBH',...
    '-v7.3')

%% 

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
save(outName,...
    'dataStruct',...
    'datasetDemo',...
    'selectNodesFrmRaw',...
    'selectNodesFrmAvgBH',...
    '-v7.3')



