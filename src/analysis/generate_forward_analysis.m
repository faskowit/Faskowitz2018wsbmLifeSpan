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
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_fit_wsbm_script_v7p3.mat');
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

%%

numPerms = 100 ;

nNodes = templateModel.Data.n ;

% save some results
permsWSBMQ = zeros([ numPerms 1 ]);
permsWSBMCCoef = zeros([ numPerms 1 ]);
permsWSBMAssort = zeros([ numPerms 1 ]);
permsWSBMParti = zeros([ numPerms 1 ]);

permsMODQ = zeros([ numPerms 1 ]);
permsMODCCoef = zeros([ numPerms 1 ]);
permsMODAssort = zeros([ numPerms 1 ]); 
permsMODParti = zeros([ numPerms 1 ]); 

permsRANDQ = zeros([ numPerms 1 ]);
permsRANDCCoef = zeros([ numPerms 1 ]);
permsRANDAssort = zeros([ numPerms 1 ]); 
permsRANDParti = zeros([ numPerms 1 ]); 

% two community defs for running modularity
[~,ca1] = community_assign(templateModel);
[~,ca2] = community_assign(modularityModel);

for idx=1:numPerms
   
    disp(idx)
    
    [~,tmpAdj] = genAdj_wsbm(templateModel);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj = tmpAdj ~= 0;
    
    permsWSBMQ(idx) = eval_modularity_wu(tmpAdj,ca1) ;
    permsWSBMCCoef(idx) = median(clustering_coef_bu(tmpAdj)) ;
    permsWSBMAssort(idx) = assortativity_bin(tmpAdj,0) ;
    permsWSBMParti(idx) = median(participation_coef(tmpAdj,ca1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     [~,tmpAdj] = genAdj_wsbm(modularityModel);
%     tmpAdj(1:nNodes+1:end)=0; %clear diagonal
%     
%     permsMODQ(idx) = eval_modularity_wu(tmpAdj,ca2) ;
%     permsMODCCoef(idx) = median(clustering_coef_bu(tmpAdj)) ;
%     permsMODAssort(idx) = assortativity_bin(tmpAdj,0) ;
%     permsMODParti(idx) = median(participation_coef(tmpAdj,ca2));

    %randomize template
    rMod = randomize_wsbm_para(templateModel,3);
    [~,tmpAdj] = genAdj_wsbm(rMod);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj = tmpAdj ~= 0;
    
    permsRANDQ(idx) = eval_modularity_wu(tmpAdj,ca1) ;
    permsRANDCCoef(idx) = median(clustering_coef_bu(tmpAdj)) ;
    permsRANDAssort(idx) = assortativity_bin(tmpAdj,0) ;
    permsRANDParti(idx) = median(participation_coef(tmpAdj,ca1));
    
end

%% empirical yo

tempBinData = templateAdj;
tempBinData = tempBinData ~= 0 ;

empQ = eval_modularity_wu(tempBinData,comVecs.yeo) ;
empCCoef = median(clustering_coef_bu(tempBinData)) ;
empAssort = assortativity_bin(tempBinData,0) ;
empParti = median(participation_coef(tempBinData,comVecs.wsbm));

%%

histogram(permsWSBMParti)
hold
histogram(permsRANDParti)


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





