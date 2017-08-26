%% load data

load('/home/jfaskowi/JOSHSTUFF/projects/sbm3/reports/reproduce_KCentral_script_data.mat')
addpath('/home/jfaskowi/JOSHSTUFF/projects/sbm3/src/external/CBIG')

%% run stuff

% get all the models with best k
kLoopBestModels2 = cell2struct(kLooperModels3(1,:), ...
    'Model', LOOPER_ITER );

% and now lets find the centroid based on pairwise variation of info
[ kCentralModel2 , kSimVec2 , kSimVUpTri2 ] = central_model(kLoopBestModels2) ; 

%% sort the kLooperModels...
% by aligning all to the 'most central' with Hungarian algo

% the ref to align to will be the kCentralModel
[~,ref2] = community_assign(kCentralModel2) ;

% preallocate the mat of the bestModels community assignments
kBestModels_ca_algn2 = zeros([ length(selectNodesFrmRaw) LOOPER_ITER ]) ;

% first stack all plausible parcellations
for idx=1:LOOPER_ITER
    
    [~,tmp] = ...
        community_assign(kLoopBestModels2(idx).Model.Para.mu) ;

    kBestModels_ca_algn2(:,idx) = CBIG_HungarianClusterMatch(ref2,tmp) ;
    
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
[ sort_kSimVec2, sort_kSimVecIdx2 ] = sort(kSimVec2,'descend');
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
kiter_prior2 = zeros([ kCentralModel2.R_Struct.k kCentralModel2.Data.n ]) ;

% sort the ca_align by similarity
kBestModels_ca_algn_sort2 = kBestModels_ca_algn2(:,sort_kSimVecIdx2);

% iterate over each node
for idx=1:(kCentralModel2.Data.n)

    kiter_prior2(:,idx) = ...
        sum(bsxfun(@eq,kBestModels_ca_algn_sort2(idx,1:cutoff), ...
        [1:(kCentralModel2.R_Struct.k)]'),2) ./ cutoff ;
        
end    

%% compare

addpath('~/JOSHSTUFF/scripts/BCT/2016_01_16_BCT/')
[~,ca1] = community_assign(kiter_prior(:,:,1)) ;
[~,ca2] = community_assign(kiter_prior2) ;

[a,b] = partition_distance(ca1,ca2) ;


