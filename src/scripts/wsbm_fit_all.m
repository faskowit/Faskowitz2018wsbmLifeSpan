%% load a config file that sets global parameters

% need to edit config file string to match what you want!!
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

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

%% pick out young subs

youngIdx = datasetDemo.age > 25 & datasetDemo.age <= 35 ;
templateSubj_data = dataStruct(youngIdx == 1) ;
clear youngIdx

%% setup a template to run wsbm on

avgTempBothHemi = make_template_mat(templateSubj_data, ...
    LEFT_HEMI_NODES, ...
    RIGHT_HEMI_NODES, ...
    MASK_THR_INIT) ; 

if strcmp(HEMI_CHOICE,'both')
    
    selectNodesFrmRaw = [LEFT_HEMI_NODES RIGHT_HEMI_NODES];
    selectNodesFrmAvgBH = 1:size(selectNodesFrmRaw) ;
    avgTemp = avgTempBothHemi .* 1;
    
elseif strcmp(HEMI_CHOICE,'left')
    
    selectNodesFrmRaw = LEFT_HEMI_NODES ;
    selectNodesFrmAvgBH = 1:size(selectNodesFrmRaw) ;
    avgTemp = avgTempBothHemi(selectNodesFrmAvgBH, ...
        selectNodesFrmAvgBH) .* 1;
    
elseif strcmp(HEMI_CHOICE,'right')
    
    selectNodesFrmRaw = RIGHT_HEMI_NODES ;
    selectNodesFrmAvgBH = ...
        (size(selectNodesFrmRaw)+1):size(avgTempBothHemi,1) ;
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
kIterOver = 5:1:12; %   

%set up the cell struct to pass to the looper
modelInputs = cell(numel(kIterOver),1); 
    
for idx = 1:numel(kIterOver)
    
    R_STRUCT_TO_TEST = sym_RStruct(kIterOver(idx)) ;
    
    modelInputs{idx} = { R_STRUCT_TO_TEST, ... 
        'W_Distr', WEIGHT_DIST, ...
        'E_Distr', EDGE_DIST, ...
        'alpha', INIT_ALPHA, ...
        'numTrials', LOOPER_NUM_TRIAL , ...
        'mainMaxIter', LOOPER_MAIN_ITER , ...
        'muMaxIter' , LOOPER_MU_ITER };
end

%% looper part

[ kLooperResults , kLooperModels ] = wsbm_looper_wrapper(avgTemp, ...
    modelInputs, ...
    LOOPER_ITER ) ;

% save this, so that we can reconstruct kLooper results if needed, ~200Mb
save(strcat(OUTPUT_DIR , '/interim/', OUTPUT_STR, '_kLooperModels.mat' ), 'kLooperModels')

%% figure our what k is best...get all these models

% find max at each row
[ ~ , maxIdx ] = max(kLooperResults) ;
kLoopBestIdx = (mode(maxIdx(2:end))) ;
kBest = kIterOver(kLoopBestIdx) ;

% get all the models with best k
kLoopBestModels = cell2struct(kLooperModels(kLoopBestIdx,:), ...
    'Model', LOOPER_ITER );

% and now lets find the median
[ kCentralModel , kSimVec ] = central_model(kLoopBestModels) ; 

%% temp

all_kSims = zeros([11 5050]);

for idx=1:11

    % get all the models with best k
    tmp = cell2struct(kLooperModels(idx,:), ...
        'Model', LOOPER_ITER );

    % and now lets find the median
    [ ~ , ~, tmp_kSimVec ] = central_model(tmp) ; 

    all_kSims(idx,:) = tmp_kSimVec';

end

%% consensus

kLooperModels_ca = zeros([ length(selectNodesFrmRaw) (size(kLooperModels,2))]) ;

% first stack all plausible parcellations
for idx=1:(size(kLooperModels,2))
    
    [~,kLooperModels_ca(:,idx)] = ...
        community_assign(kLoopBestModels(idx).Model.Para.mu) ;

end

% get agreement mat
kLooperModels_agreement = agreement(kLooperModels_ca) ;
kLooperModels_agreement = kLooperModels_agreement ./ (size(kLooperModels,2));
kLooperModels_consensus = consensus_und(kLooperModels_agreement,0.02,1000) ;

% versatility
cm = kLooperModels_agreement .* 1 ;
cs = sin(pi*cm) ;
versatility = sum(cs,1)./sum(cm,1) ;

%% 

logEvid = zeros([1 (size(kLooperModels,2))]) ;

for idx=1:(size(kLooperModels,2))
    
    logEvid(idx) = kLoopBestModels(idx).Model.Para.LogEvidence ;

end

figure
histogram(logEvid)
hold on
plot(templateModel.Para.LogEvidence * ones(1,50), 0:49);

weight_agreement = agreement_weighted(kLooperModels_ca,mat2gray(logEvid)) ;
weight_agreement(1:(size(weight_agreement,1)+1):end) = 0 ;

% show that there is corrleation between logEvid and VI disatnce
logEvid_kSim_corr = corr(logEvid',kSimVec) ;

%% what do the models close to the central model
% look like?

[ sort_kSimVec, sortIdx ] = sort(kSimVec,'descend');

% get top five
top5Model = kLoopBestModels(sortIdx(1:5)) ; 
top5ca = zeros([size(kCentralModel.Data.Raw_Data,1) 5]);
[~,ref] = community_assign(kCentralModel);

addpath('aux_stuff/')

for idx=1:5
   
    disp(top5Model(idx).Model.Para.LogEvidence) 
    [~,tmp] = community_assign(top5Model(idx).Model);
    top5ca(:,idx) = CBIG_HungarianClusterMatch(ref,tmp) ;
    
end

%% get centroid of best k and define it as templateModel

templateModel = kCentralModel ; 

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
    [ ~ , maxIdx ] = max(alphaLooperResults) ;
    alphaLoopBestIdx = (mode(maxIdx(2:end))) ;
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

%% with the iterative fits

% lets not do this with the grad descent exactly right now....
maxIters = 5 ; 
%minIters = 5 ; 
modelFits = 5 ; 
%convergence = -0.00001 ; 

%for i = 1:length(raw_data)
parallel_pool = gcp ; 
parfor subj = 1:datasetSize
    
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
    variableMuPrior = muPrior ;
    
    %% iterative fitting the wsbm
    for idx=1:maxIters
    
        % keep this in the loop, might need to for variableMu
        indivModelInputs = {'W_Distr', w_dist_to_use, ...
                    'E_Distr', e_dist_to_use, ...
                    'numTrials', INDIV_WSBM_NUM_TRIAL , ...
                    'mu_0', variableMuPrior , ...
                    'alpha',alpha_to_use , ...
                    'mainMaxIter', INDIV_WSBM_MAIN_ITER , ...
                    'muMaxIter', INDIV_WSBM_MU_ITER,  ...
                    'verbosity' , 0 }; 

        indivCentModel = wsbmCentralFit( subjAdjMat, ...
            (sym_RStruct(k_to_use)), ...
            indivModelInputs, ...
            modelFits, ...
            variableMuPrior );
    
        [ variableMuPrior , ~ ] = make_WSBM_prior(indivCentModel, PRIOR_WEIGHT_MORE) ;

        % record this iter's central model
        fitWSBMAllStruct(subj).Model(idx) = indivCentModel ;

    end
   
end % looping over each subject

%% filter sub-optimal data here
% condition the data based on exclusion criteria

nSubj = length(dataStruct) ;
sparse_cutoff = 0.25 ;
removeVec = zeros( [nSubj 1] );

for idx=1:nSubj
   
    tmpSubjMat = ...
    dataStruct(idx).countVolNormMat(selectNodesFrmRaw, ...
    selectNodesFrmRaw ) ;
    
    sparseness = density_und(tmpSubjMat);
    
    if sparseness < sparse_cutoff
       
        removeVec(idx) = 1 ;
    end
end

% take out bad subjects
% by filters by subjects not marked for removal
dataStruct = dataStruct(removeVec == 0);
datasetDemo = datasetDemo(removeVec == 0, :);
fitWSBMAllStruct = fitWSBMAllStruct(removeVec == 0) ;

%% lets look real quick

logEvid = zeros([ size(fitWSBMAllStruct,2) 1 ]) ;

for idx=1:(size(fitWSBMAllStruct,2))
   
    logEvid(idx) = fitWSBMAllStruct(idx).Model(5).Para.LogEvidence ;
end

%% end of fitting WSBM on brains
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





