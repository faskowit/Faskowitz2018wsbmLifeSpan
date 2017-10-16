%% load a config file that sets global parameters

% need to edit config file string to match what you want!!
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

%% 

workspaceIWant = strcat(PROJECT_DIR, ...
    '/data/processed/',....
    OUTPUT_STR,...
    '_fit_wsbm_script.mat') ;

load(workspaceIWant)

% %% pick out young subs
% 
% youngIdx = datasetDemo.age > 25 & datasetDemo.age <= 35 ;
% templateSubj_data = dataStruct(youngIdx == 1) ;
% clear youngIdx

%% split dataset into age ranges!
% split into approx equal sized parts! 

% 620 subjects / 5 groups = 124 per group...

%we should make some age bins 
thresholds = 4 ; % meaning 5 groups! 

[agesSorted,ageSortIdx] = sort(datasetDemo.age) ;

low_quantile = [ 0 quantile(1:length(datasetDemo.age),4) ] ; 
high_quantile = [ quantile(1:length(datasetDemo.age),4) (length(datasetDemo.age)+1) ] ;

templateIdMat = zeros( [ size(datasetDemo.age,1) (thresholds+1) ] ) ;

for templateIdx=1:(thresholds+1)

    subjectsIdxVec = ageSortIdx > low_quantile(templateIdx) & ...
        ageSortIdx < high_quantile(templateIdx) ;
    
    templateIdMat(:,templateIdx) = subjectsIdxVec ;
end

s = length(LEFT_HEMI_NODES) + length(RIGHT_HEMI_NODES);
avgTempMat = zeros([ s s (thresholds+1)]);

%% setup the structurs that will save all the data

kLooperResultsCrossIdx = cell([(thresholds+1) 1]) ;
kLooperModelsCrossIdx = cell([(thresholds+1) 1]) ;
kCentralModelCrossIdx = cell([(thresholds+1) 1]) ;
consensusCrossIdx =  cell([(thresholds+1) 1]) ;

%% now the loop!!!!

for templateIdx=1:(thresholds+1)

    templateSubj_data = dataStruct(templateIdMat(:,templateIdx) == 1) ;

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

    %set up the cell struct to pass to the looper
    modelInputs = cell(numel(kIterOver),1); 

    for idx = 1:numel(kIterOver)

        R_STRUCT_TO_TEST = sym_RStruct(kIterOver(idx)) ;

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

    %% k looper part
    % iterate over num communities, k, to find out which k give best evidence

    [ kLooperResults , kLooperModels ] = wsbm_looper_wrapper(avgTemp, ...
        modelInputs, ...
        LOOPER_ITER, ...
        [], numTrialPtrn,...
        prirPtrn) ;
    
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

    %% lets save results from this looped run!!!
    
    kLooperResultsCrossIdx{templateIdx} = kLooperResults ;
    kLooperModelsCrossIdx{templateIdx} =  kLooperModels ;
    kCentralModelCrossIdx{templateIdx} = kCentralModel ;
    consensusCrossIdx{templateIdx} = consensus_kCentralModel ;
    
end

%% make a data struct to save info...

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_agebin_template_results_v7p3.mat');
save(outName,...
    'kLooperResultsCrossIdx',...
    'kLooperModelsCrossIdx',...
    'kCentralModelCrossIdx',...
    'consensusCrossIdx',...
    '-v7.3')

