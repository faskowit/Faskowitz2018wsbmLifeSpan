
%% INITALIZE PROJECT
% evevtually it would be nice to make this into a main function to be
% called by a bash script

projectDir = '/home/jfaskowi/JOSHSTUFF/projects/SBM2/' ;
mkdir(projectDir)
cd(projectDir)

projectFolders = { 
    'project_data' 
    'raw_data'
    'results_data'
    'scripts'
    'WSBM_v1.2/'
    } ;

for idx=1:length(projectFolders)
    mkdir(projectFolders{idx})
end

executionFolders = {
    'scripts'
    'WSBM_v1.2'
    'WSBM_v1.2/analysis_tools/'
    'analysis'
    } ;

for idx=1:length(executionFolders)
    addpath(strcat(projectDir,'/',executionFolders{idx}))
end

addpath('~/JOSHSTUFF/scripts/BCT/2016_01_16_BCT/') 

dataDir = '/raw_data/NKIen1/';

%% SETUP GLOBAL VARS
% global vars and directions for stuff below... could be read in by via an
% initalization file perhaps

% indication which parcellation being used
PARCELLATION = 'yeo' ;

if strcmp(PARCELLATION,'yeo')
% index of the nodes for each hemisphere, as we initially read in data that
% has '0' val and also subcortical structures
LEFT_HEMI_NODES = 2:58 ; %1:57 ;
RIGHT_HEMI_NODES = 60:116 ; %75:131 ; 
end

% WSBM options
% distributions 
WEIGHT_DIST = 'normal' ; 
EDGE_DIST = 'poisson' ; 

% options for wsbm.m
INIT_ALPHA = 0.5 ; 

% model wsbm run
MODEL_WSBM_NUM_TRIAL = 500 ; 
MODEL_WSBM_MAIN_ITER = 125 ;
MODEL_WSBM_MU_ITER = 50 ;

% individual wsbm run 
INDIV_WSBM_NUM_TRIAL = 150 ;
INDIV_WSBM_MAIN_ITER = 125 ;
INDIV_WSBM_MU_ITER = 50 ;

% looper vars when wsbm looped
LOOPER_ITER = 100 ; 
LOOPER_NUM_TRIAL = 125 ; 
LOOPER_MAIN_ITER = 125 ;
LOOPER_MU_ITER = 50 ;

% GLOBAL VARS FOR SCRIPT
PRIOR_WEIGHT_MORE = 1.0 ;
MASK_THR_INIT = 1 ;
MASK_THR = 1 ;

% 'both' 'left' right'
HEMI_CHOICE = 'both' ; 

outputStr = strcat('/', ...
    PARCELLATION , '_', ...
    HEMI_CHOICE, '_', ...
    WEIGHT_DIST, ...
    EDGE_DIST, '_', ...
    'a', strrep(num2str(INIT_ALPHA),'.','p'), ...
    '/') ;

EXPERIMENT = 'seedWyeo' ;
OUTPUT_DIR = strcat(projectDir , '/results_data/', EXPERIMENT, '/', outputStr) ; 
mkdir(OUTPUT_DIR)

%% load the NKI data

readData = read_nki_data( strcat(projectDir,dataDir) , PARCELLATION) ; 

dataStruct = readData.dataRaw ;
datasetDemo = readData.demoRaw ;

%%

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

%% use the Yeo 7 as the mu prior

load('aux_stuff/seven_network.mat')
% add the temporal parietal to default mode
mu_yeo7 = seven_network ;
mu_yeo7(7,57) = 1 ;
mu_yeo7(7,114) = 1 ;

%% make the prior

muPrior = make_WSBM_prior(mu_yeo7 , PRIOR_WEIGHT_MORE) ;

%% all the subjects

disp('now working on fitting all subjects together')

% get number of subjects to run
datasetSize = length(dataStruct) ; 

% initialize struct
fitWSBMAllStruct = struct() ; 

% parameters for WSMB individual fit
k_to_use = size(mu_yeo7,1)  ;
alpha_to_use = INIT_ALPHA ;
w_dist_to_use = WEIGHT_DIST ;
e_dist_to_use = EDGE_DIST ;

%% with the iterative fits

% lets not do this with the grad descent exactly right now....
maxIters = 5 ; 
%minIters = 5 ; 
modelFits = 5 ; 
%convergence = -0.00001 ; 

%for i = 1:length(raw_data)
parallel_pool = gcp ; 
parfor subj = 1:datasetSize
%for subj= 2:2
    
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
indivModelInputStruct.('alpha') = 0.5 ;
indivModelInputStruct.('mainMaxIter') = INDIV_WSBM_MAIN_ITER ;
indivModelInputStruct.('muMaxIter') = INDIV_WSBM_MU_ITER ;

addNote = 'seedYeo7' ;
outName = strcat(projectDir,'/workspaces/', addNote, 'fit_wsbm_script_workspace.mat');
save(outName,...
    'dataStruct',...
    'datasetDemo',...
    'muPrior',...
    'fitWSBMAllStruct',...
    'templateModel',...
    'indivModelInputStruct',...
    'projectDir',...
    'dataDir'...
    )
