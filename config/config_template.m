%% INITALIZE PROJECT
% evevtually it would be nice to make this into a main function to be
% called by a bash script

PROJECT_DIR = '/home/jfaskowi/JOSHSTUFF/projects/sbm3/' ;
mkdir(PROJECT_DIR)
cd(PROJECT_DIR)

projectFolders = { 
    'src' 
    'data'
    'results_data'
    'bin'
} ;

for idx=1:length(projectFolders)
    addpath(genpath(strcat(PROJECT_DIR,'/',projectFolders{idx})))
end

DATA_DIR = '/data/raw/NKIen1/';

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
elseif strcmp(PARCELLATION,'scale125')
    LEFT_HEMI_NODES = 116:226 ; %1:57 ;
    RIGHT_HEMI_NODES = 2:108 ; %75:131 ; 
end

% WSBM options
% distributions 
WEIGHT_DIST = 'normal' ; 
%EDGE_DIST = 'poisson' ; 
EDGE_DIST = setup_distr('poisson',[0.1,0.01]) ;

% options for wsbm.m
INIT_ALPHA = 0.5 ; 

% model wsbm run
MODEL_WSBM_NUM_TRIAL = 500 ; 
MODEL_WSBM_MAIN_ITER = 125 ;
MODEL_WSBM_MU_ITER = 50 ;

% individual wsbm run 
% INDIV_WSBM_NUM_TRIAL = 2000 ;
INDIV_WSBM_MAIN_ITER = 100 ;
INDIV_WSBM_MU_ITER = 50 ;

% looper vars when wsbm looped
LOOPER_ITER = 100 ; 
LOOPER_NUM_TRIAL = 2000 ; 
LOOPER_MAIN_ITER = 100 ;
LOOPER_MU_ITER = 50 ;

% GLOBAL VARS FOR SCRIPT
PRIOR_WEIGHT_MORE = 1.0 ;
MASK_THR_INIT = 1 ;
MASK_THR = 1 ;

% 'both' 'left' right'
HEMI_CHOICE = 'both' ; 

OUTPUT_STR = strcat(...
    PARCELLATION , '_', ...
    HEMI_CHOICE, '_', ...
    WEIGHT_DIST, ...
    EDGE_DIST, '_', ...
    'a', strrep(num2str(INIT_ALPHA),'.','p')) ;

%% make output dir

OUTPUT_DIR = strcat(PROJECT_DIR , '/data/') ; 
mkdir(OUTPUT_DIR)
