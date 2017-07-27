%% lets see how the average matrix compares to all the 'ya' subjects

tempData = templateModel.Data.Raw_Data ;
tempData(~~isnan(tempData)) = 0 ;

% templateSubj_data: data from subject who went into template

%% with the iterative fits

% lets not do this with the grad descent exactly right now....
maxIters = 1 ; 
%minIters = 5 ; 
modelFits = 100 ; 
%convergence = -0.00001 ; 

%for i = 1:length(raw_data)
parallel_pool = gcp ; 
parfor subj = 1:(size(templateSubj_data,2))
%for subj = 1:1 

    disp('working on')
    disp(subj)
    
    %% read in individual data
        
    % loop through all the data and fit the blockmodel with the prior
    subjAdjMat = templateSubj_data(subj).countVolNormMat(selectNodesFrmRaw, selectNodesFrmRaw);
    
    % get rid of the diagonal
    n=size(subjAdjMat,1);
    subjAdjMat(1:n+1:end) = 0;
    
    % mask out AdjMat entries below mask_thr
    subjAdjMat_mask = templateSubj_data(subj).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > MASK_THR ;    
    subjAdjMat_mask(subjAdjMat_mask > 0) = 1 ;   
    subjAdjMat = subjAdjMat .* subjAdjMat_mask ;
    
    % replace the 0's with NaN
    subjAdjMat(subjAdjMat == 0) = NaN ;

    %% record stuff for the output struct
    
    % reord subj id in results struct
    fitWSBMtemplateSubs(subj).id = templateSubj_data(subj).id ;
    
    % save the raw data only once
    fitWSBMtemplateSubs(subj).Raw_Data = subjAdjMat ; 
    
%     % initialize the variable prior
%     variableMuPrior = muPrior ;
    
    %% iterative fitting the wsbm
    for idx=1:maxIters
    
        % keep this in the loop, might need to for variableMu
        indivModelInputs = {'W_Distr', w_dist_to_use, ...
                    'E_Distr', e_dist_to_use, ...
                    'numTrials', INDIV_WSBM_NUM_TRIAL , ...
                    'alpha',alpha_to_use , ...
                    'mainMaxIter', INDIV_WSBM_MAIN_ITER , ...
                    'muMaxIter', INDIV_WSBM_MU_ITER,  ...
                    'verbosity' , 0 }; 

        [indivCentModel,allModels] = wsbmCentralFit( subjAdjMat, ...
            (sym_RStruct(k_to_use)), ...
            indivModelInputs, ...
            modelFits );
        
        % record this iter's central model
        fitWSBMtemplateSubs(subj).Model(idx) = indivCentModel ;
        fitWSBMtemplateSubs(subj).allModels = allModels ;
        
    end
   
end % looping over each subject

