%% load data

load('/home/jfaskowi/JOSHSTUFF/projects/sbm3/reports/reproduce_KCentral_script_data.mat')
addpath('/home/jfaskowi/JOSHSTUFF/projects/sbm3/src/external/CBIG')

%% run stuff

% get all the models with best k
kLoopBestModels_reproduce = cell2struct(kLooperModels3(1,:), ...
    'Model', LOOPER_ITER );

% and now lets find the centroid based on pairwise variation of info
[ kCentralModel_repo , kSimVec_repo , kSimVUpTri_repo ] = central_model(kLoopBestModels_reproduce) ; 

%% sort the kLooperModels...
% by aligning all to the 'most central' with Hungarian algo

% the ref to align to will be the kCentralModel
[~,ref2] = community_assign(kCentralModel_repo) ;

% preallocate the mat of the bestModels community assignments
kBestModels_ca_algn_repo = zeros([ length(selectNodesFrmRaw) LOOPER_ITER ]) ;

% first stack all plausible parcellations
for idx=1:LOOPER_ITER
    
    [~,tmp] = ...
        community_assign(kLoopBestModels_reproduce(idx).Model.Para.mu) ;

    kBestModels_ca_algn_repo(:,idx) = CBIG_HungarianClusterMatch(ref2,tmp) ;
    
end

%% get kSimVec and LogEvid vec sorted
% 
% % will be size: numKTested x LOOPER_ITER
% kLooperLogEvids = cellfun(@(a) a.Para.LogEvidence,kLooperModels(:,:)) ;
% % get the row with the best evids, and transpose to make col
% kBestLogEvids = kLooperLogEvids(kLoopBestIdx,:)' ;
% 
% % sort evid vector and kSim vector
% [ sort_logEvid, sort_logEvidIdx ] = sort(kBestLogEvids,'descend');
[ sort_kSimVec2, sort_kSimVecIdx2 ] = sort(kSimVec_repo,'descend');
% 
% % show that there is corrleation between logEvid and VI disatnce
% logEvid_kSim_corr = corr(kBestLogEvids,kSimVec) ;

%% lets iterate over the 100 k best models to create consensus model
% first we have to turn the hungarian reordered ca_align mat into a wsbm
% prior to be input into wsbm inference fit

% % lets cuttoff the worst 5% of the runs to clear up the prior a little...
% % TODO... maybe no need for this...
onlyKeep = 1; %0.95 ;

% first figure out the cuttoff
cutoff = floor(LOOPER_ITER * onlyKeep) ; 
kiter_prior_repo = zeros([ kCentralModel_repo.R_Struct.k kCentralModel_repo.Data.n ]) ;

% sort the ca_align by similarity
kBestModels_ca_algn_repo = kBestModels_ca_algn_repo(:,sort_kSimVecIdx2);

% iterate over each node
for idx=1:(kCentralModel_repo.Data.n)

    kiter_prior_repo(:,idx) = ...
        sum(bsxfun(@eq,kBestModels_ca_algn_repo(idx,1:cutoff), ...
        [1:(kCentralModel_repo.R_Struct.k)]'),2) ./ cutoff ;
        
end    

%% compare

addpath('~/JOSHSTUFF/scripts/BCT/2017_01_15_BCT/')
[~,ca1] = community_assign(kiter_prior(:,:,1)) ;
[~,ca2] = community_assign(kiter_prior_repo) ;

[a,b] = partition_distance(ca1,ca2) ;

%% also if we run the consensus on this new kiter_prior_repo...
% will be we close to the previous consensus

CONSENSUS_ITER = 10 ;
C_repo = zeros([CONSENSUS_ITER 1]) ;

for idx=1:CONSENSUS_ITER

    runNTimes = 100 ;
    
    rr = sym_RStruct(kCentralModel_repo.R_Struct.k) ;
    modIn = { ... 
        'W_Distr', kCentralModel_repo.W_Distr.name, ...
        'E_Distr', kCentralModel_repo.E_Distr.name, ...
        'alpha', kCentralModel_repo.Options.alpha, ...
        'mainMaxIter', 100 , ...
        'muMaxIter' , 50,...
        'mainTol',0.01, ...
        'muTol', 0.01 ,...
        'verbosity', 0, ...
        'numTrials', 50 ,...
        'mu_0', kiter_prior_repo(:,:,idx)};
    
    % function [ allModels ] = wsbmFitNTimes( adjMat, rStruct , modelInputs , numFits , numCores)
    cnsnsusModels = wsbmFitNTimes(kCentralModel_repo.Data.Raw_Data,...
        rr,...
        modIn,...
        runNTimes, 16) ;
        
    tmpCnsnsusCa = zeros([ kCentralModel_repo.Data.n runNTimes ]) ;
    
    for jdx=1:runNTimes
        [~,tmpCnsnsusCa(:,jdx)] = community_assign(cnsnsusModels(jdx).Model) ;
    end
   
    tmpAgreeMat = agreement(tmpCnsnsusCa) ./ runNTimes ;
    
    % get the consensus consitency
    C_repo(idx) = consensus_consistency(tmpAgreeMat) ;
    
    % make new kiter_prior for new loop
    for kdx=1:(kCentralModel_repo.Data.n)

        kiter_prior_repo(:,kdx,idx+1) = ...
            sum(bsxfun(@eq,tmpCnsnsusCa(kdx,:), ...
            [1:(kCentralModel_repo.R_Struct.k)]'),2) ./ cutoff ;
        
    end    
    
    % have we converged? or are we at the end of loop?
    if C_repo(idx) >= 0.95 || idx == CONSENSUS_ITER
        consensus_kiter_prior_repo = kiter_prior_repo(:,:,idx+1) ;
        [~,consensus_ca_repo] = community_assign(consensus_kiter_prior_repo) ; 
        consensus_kCentralModel_repo = central_model(cnsnsusModels) ;
        % also add the data back into consensus_kCentral...
        consensus_kCentralModel_repo.Data = kCentralModel_repo.Data ;
        break 
    end
    
end

%% now lets check how this consensus looks compared to first consensus

consensus_kCentralModel = load('/home/jfaskowi/JOSHSTUFF/projects/sbm3/data/interim/yeo_both_normalpoisson_a0p5_templateModel_1.mat');
consensus_kCentralModel = consensus_kCentralModel.templateModel ;

[~,con_ca1] = community_assign(consensus_kCentralModel) ;
[~,con_ca2] = community_assign(consensus_kCentralModel_repo) ;

[a,b] = partition_distance(con_ca1,con_ca2) ;




