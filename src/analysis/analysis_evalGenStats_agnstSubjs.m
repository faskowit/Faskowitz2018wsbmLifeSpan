%% clear stuff

clc
clearvars

%% load the necessary data

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

% loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_fit_wsbm_script_v7p3.mat');
% load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_templateModel_1.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
load(loadName) ;

FIGURE_NAME = 'figC2' ;

outputdir = strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/');
mkdir(outputdir)


%% actual data

templateAdj = templateModel.Data.Raw_Data ;
templateAdj(~~isnan(templateAdj)) = 0 ;

%% and the modular model 

muMod = dummyvar(comVecs.mod)' ;
[~,modularityModel] = wsbm(templateModel.Data.Raw_Data, ...
    templateModel.R_Struct.R, ...
    'W_Distr', templateModel.W_Distr, ...
    'E_Distr', templateModel.E_Distr, ...
    'alpha', templateModel.Options.alpha, ...
    'mu_0', muMod , ...
    'verbosity', 0);

tmpYeo = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.yeo) ;
muYeo = dummyvar(tmpYeo)' ;
muYeo = muYeo(sum(muYeo,2)>0,:) ;

[~,yeoModel] = wsbm(templateModel.Data.Raw_Data, ...
    sym_RStruct(7), ...
    'W_Distr', templateModel.W_Distr, ...
    'E_Distr', templateModel.E_Distr, ...
    'alpha', templateModel.Options.alpha, ...
    'mu_0', muYeo , ...
    'verbosity', 0);

%% try out the evalWSBM code

nNodes = templateModel.Data.n ;

templateSubj_data = dataStruct(datasetDemo.age > ageLowLim & datasetDemo.age <= ageHighLim) ;
[a,b,avgTemp_dist] = make_template_mat(templateSubj_data, ...
    LEFT_HEMI_NODES, ...
    RIGHT_HEMI_NODES, ...
    MASK_THR_INIT) ; 

% actually replace the 0's with NaN
%avgTemp(avgTemp == 0) = NaN ;
avgTemp_dist = avgTemp_dist(selectNodesFrmRaw,selectNodesFrmRaw);
% clear diagonal
avgTemp_dist(1:nNodes+1:end)=0; 

subjDataMat = zeros([ nNodes nNodes length(templateSubj_data) ]);
for idx = 1:length(templateSubj_data)  

    tmpAdj = templateSubj_data(idx).countVolNormMat(selectNodesFrmRaw, selectNodesFrmRaw);
    % get rid of the diagonal
    %n=size(tmpAdj,1);
    tmpAdj(1:nNodes+1:end) = 0; 
    % mask out AdjMat entries below mask_thr
    tmpAdj_mask = templateSubj_data(idx).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > MASK_THR ;    
    tmpAdj_mask(tmpAdj_mask > 0) = 1 ;   
    tmpAdj = tmpAdj .* tmpAdj_mask ;
    tmpAdj(isnan(tmpAdj)) = 0 ;
    subjDataMat(:,:,idx) = tmpAdj ;
end

%% eval gen call

% setup some structs to record results
eval_subj_wsbm_K = cell([ length(templateSubj_data)  1]);
eval_subj_mod_K = cell([ length(templateSubj_data)  1]);

parfor idx = 1:length(templateSubj_data)  
    
    currentSubj = subjDataMat(:,:,idx) ;
    
    % function [B,E,K] = eval_genWsbm_model(wsbmModel,D,numSims)]
    %           B,          n x n x numSims matrix of synthetic networks
    %           E,          energy for each synthetic network
    %           K,          Kolmogorov-Smirnov statistics for each synthetic
    %                       network.

    [~,~,eval_subj_wsbm_K{idx}] = eval_genWsbm_model1_agnstSubj(templateModel,avgTemp_dist,5000,0,currentSubj);
    [~,~,eval_subj_mod_K{idx}] = eval_genWsbm_model1_agnstSubj(modularityModel,avgTemp_dist,5000,0,currentSubj);

    disp(idx)

end

%% collate results

eval_subj_wsbm_K_mat = cell2mat(eval_subj_wsbm_K) ;
eval_subj_mod_K_mat = cell2mat(eval_subj_mod_K) ;

histogram(mean(eval_subj_wsbm_K_mat,2))
hold
histogram(mean(eval_subj_mod_K_mat,2))

% also compare histograms between each 53 subject
for idx = 1:length(templateSubj_data)
    
    subj_diff(idx) = median(mean(eval_subj_wsbm_K{idx},2)) - ...
        median(mean(eval_subj_mod_K{idx},2)) ;

end

%% but also how to display this...?

% load(strcat(PROJECT_DIR,'/data/processed/',OUTPUT_STR,'_evalGenReps_subj.mat'))
% save(strcat(PROJECT_DIR,'/data/processed/',OUTPUT_STR,'_evalGenReps_subj.mat'),...
%     'eval_subj_wsbm_K',...
%     'eval_subj_mod_K' ...
%      )










