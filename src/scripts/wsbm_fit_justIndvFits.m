%% load a config file that sets global parameters

% need to edit config file string to match what you want!!
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
config_file='config_yeo_indivWLessPrior.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

%% load the data we previously computed

% loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
% load(loadName) ;

load('data/processed/yeo_both_normalpoisson_a0p5_basicData_v7p3.mat')
load('data/interim/yeo_both_normalpoisson_a0p5_templateModel_1.mat')

alphaBest = 0.5 ;

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
r_struct_to_use = sym_RStruct(templateModel.R_Struct.k) ;

%% perhaps changing this up...

modelFits = 5 ; 
maxIters = 1 ; 

% reduced this weight from 3 to 1
[ muPrior , ~ ] = make_WSBM_prior(templateModel , 1) ;

indiv_numTrialPtrn = [ 250 100 100 100 100 100 100 100 100 100 100 ] ;
% this will be one shorter because there is no prior weight for first run
% just put a 0 at the end so it doesnt mess it up
% prirPtrn = [ 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 0 ] ;
indiv_prirPtrn = [ 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 0 ] ;

test_times = zeros([10 1]) ;

%for subj = 1:3
parallel_pool = gcp ; 

% ppm = ParforProgMon(strWindowTitle, nNumIterations <, nProgressStepSize, nWidth, nHeight>);
ppm = ParforProgMon('indivFits',620,1) ;

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

        %function [ Model ] = wsbmFitWPttrn( adjMat, rStruct , modelInputs , initMu, numTrialPttrn, priorWeightPttrn)
        tempModelStruct(idx).Model = wsbmFitWPttrn( subjAdjMat, r_struct_to_use , ...
            indivModelInputs , muPrior, indiv_numTrialPtrn, indiv_prirPtrn) ;
        
        fitWSBMAllStruct(subj).Model(idx) = tempModelStruct(idx).Model ;

    end
   
    %recrd the central model
    fitWSBMAllStruct(subj).centModel = central_model(tempModelStruct) ;
    
    test_times(subj) = toc(t1) ;
    
    ppm.increment() ;
    
end % looping over each subject

%% save

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_fit_wsbm_script_v7p3.mat');
save(outName,...
    'dataStruct',...
    'datasetDemo',...
    'fitWSBMAllStruct',...
    '-v7.3')


