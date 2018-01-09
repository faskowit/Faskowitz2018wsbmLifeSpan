%% clear stuff

clc
clearvars

%% load the necessary data
%% lets analyze the age bin scripts

addpath('~/JOSHSTUFF/scripts/BCT/2017_01_15_BCT/')

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

% load the data we need to analyze this ish. 
% loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_fit_wsbm_script_v7p3.mat');
% load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_templateModel_1.mat');
load(loadName) ;

% load the data we need to analyze this ish. 
loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
load(loadName) ;

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

% check the empircal edge existence
% [~,~,~,tmpMat] = get_block_mat(templateAdj,comVecs.mod);

%% setup vars for perm tests
% select some metrics that can be reduce to scalar, to get distributions
% across many permutations

numPerms = 100 ;

nNodes = templateModel.Data.n ;

% save some results
permsWSBM_Q = zeros([ numPerms 1 ]);
permsWSBM_CCoef = zeros([ numPerms 1 ]);
permsWSBM_Assort = zeros([ numPerms 1 ]);
permsWSBM_Parti = zeros([ numPerms 1 ]);

permsMOD_Q = zeros([ numPerms 1 ]);
permsMOD_CCoef = zeros([ numPerms 1 ]);
permsMOD_Assort = zeros([ numPerms 1 ]); 
permsMOD_Parti = zeros([ numPerms 1 ]); 

permsRAND_Q = zeros([ numPerms 1 ]);
permsRAND_CCoef = zeros([ numPerms 1 ]);
permsRAND_Assort = zeros([ numPerms 1 ]); 
permsRAND_Parti = zeros([ numPerms 1 ]); 

%% iterate through permuations

for idx=1:numPerms
   
    disp(idx)
    
    [~,tmpAdj] = genAdj_wsbm(templateModel);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj = tmpAdj ~= 0;
    
    permsWSBM_Q(idx) = eval_modularity_wu(tmpAdj,comVecs.wsbm) ;
    permsWSBM_CCoef(idx) = median(clustering_coef_bu(tmpAdj)) ;
    permsWSBM_Assort(idx) = assortativity_bin(tmpAdj,0) ;
    permsWSBM_Parti(idx) = median(participation_coef(tmpAdj,comVecs.wsbm));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,tmpAdj] = genAdj_wsbm(modularityModel);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj = tmpAdj ~= 0;

    permsMOD_Q(idx) = eval_modularity_wu(tmpAdj,comVecs.mod) ;
    permsMOD_CCoef(idx) = median(clustering_coef_bu(tmpAdj)) ;
    permsMOD_Assort(idx) = assortativity_bin(tmpAdj,0) ;
    permsMOD_Parti(idx) = median(participation_coef(tmpAdj,comVecs.mod));

    % randomize template
    randModel = randomize_wsbm_para(templateModel,3);
    [~,tmpAdj] = genAdj_wsbm(randModel);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj = tmpAdj ~= 0;
    
    permsRAND_Q(idx) = eval_modularity_wu(tmpAdj,comVecs.wsbm) ;
    permsRAND_CCoef(idx) = median(clustering_coef_bu(tmpAdj)) ;
    permsRAND_Assort(idx) = assortativity_bin(tmpAdj,0) ;
    permsRAND_Parti(idx) = median(participation_coef(tmpAdj,comVecs.wsbm));
    
end

%% empirical yo

tempBinData = templateAdj;
tempBinData = tempBinData ~= 0 ;

empQ = eval_modularity_wu(tempBinData,comVecs.yeo) ;
empCCoef = median(clustering_coef_bu(tempBinData)) ;
empAssort = assortativity_bin(tempBinData,0) ;
empParti = median(participation_coef(tempBinData,comVecs.wsbm));

%%

histogram(permsWSBM_Parti)
hold
histogram(permsRAND_Parti)


%%

assortVec = zeros([nSubj 1]);

for idx=1:nSubj
   
    assortVec(idx) = assortativity_bin(wsbm_gatherStruct.subjMatsArray(:,:,idx),0);
    
    
end

%% try out the evalWSBM code

nNodes = templateModel.Data.n ;

templateSubj_data = dataStruct(datasetDemo.age > 25 & datasetDemo.age <= 35) ;
[a,b,avgTemp_dist] = make_template_mat(templateSubj_data, ...
    LEFT_HEMI_NODES, ...
    RIGHT_HEMI_NODES, ...
    MASK_THR_INIT) ; 

% actually replace the 0's with NaN
%avgTemp(avgTemp == 0) = NaN ;
avgTemp_dist = avgTemp_dist(selectNodesFrmRaw,selectNodesFrmRaw);
% clear diagonal
avgTemp_dist(1:nNodes+1:end)=0; 

%% 

% function [B,E,K] = eval_genWsbm_model(wsbmModel,D,numSims)]
%           B,          n x n x numSims matrix of synthetic networks
%           E,          energy for each synthetic network
%           K,          Kolmogorov-Smirnov statistics for each synthetic
%                       network.
[evalB,evalE,evalK] = eval_genWsbm_model1(templateModel,avgTemp_dist,1000);
[evalMODB,evalMODE,evalMODK] = eval_genWsbm_model1(modularityModel,avgTemp_dist,1000);

%% 

histogram(evalE,'normalization','probability') 
hold 
histogram(evalMODE,'normalization','probability')
legend('WSBM model','Modularity model')





