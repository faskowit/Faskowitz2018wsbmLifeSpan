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

% % empircal edge existence and edge weights
% [~,modModel_w,~,modModel_e] = get_block_mat(templateAdj,comVecs.mod);
% 
% modModel_w = nonzeros(triu(modModel_w)) ;
% modModel_e = nonzeros(triu(modModel_e)) ;
% 
% modularityModel.Para.predict_w = modModel_w ;
% modularityModel.Para.predict_e = modModel_e ;

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

%% eval gen call

% function [B,E,K] = eval_genWsbm_model(wsbmModel,D,numSims)]
%           B,          n x n x numSims matrix of synthetic networks
%           E,          energy for each synthetic network
%           K,          Kolmogorov-Smirnov statistics for each synthetic
%                       network.
[eval_wsbm_B,eval_wsbm_E,eval_wsbm_K,eval_wsbm_EMD] = eval_genWsbm_model1(templateModel,avgTemp_dist,10000,0);
[eval_mod_B,eval_mod_E,eval_mod_K,eval_mod_EMD] = eval_genWsbm_model1(modularityModel,avgTemp_dist,10000,0);

[eval_wsbmRand_B,eval_wsbmRand_E,eval_wsbmRand_K] = eval_genWsbm_model1(templateModel,avgTemp_dist,10000,1);
[eval_modRand_B,eval_modRand_E,eval_modRand_K] = eval_genWsbm_model1(modularityModel,avgTemp_dist,10000,1);

%% plot energy 

cmap = [0    0.4470    0.7410 ;
    0.8500    0.3250    0.0980];

histogram(mean(eval_wsbm_K,2),...
    'normalization','probability',...
    'FaceColor',cmap(1,:),'EdgeAlpha',0.01) 
hold 
histogram(mean(eval_mod_K,2),...
    'normalization','probability',...
    'FaceColor',cmap(2,:),'EdgeAlpha',0.01)
legend('WSBM','Modular')

%[~,p,ci,stat] = ttest2(mean(eval_wsbm_K,2),mean(eval_mod_K,2))
% p =
%      0
% ci =
%    -0.0128
%    -0.0123
% stat = 
%     tstat: -88.6506
%        df: 19998
%        sd: 0.0100
%
% [~,p,ci,stat] = ttest2(mean(eval_wsbm_K,2),mean(eval_mod_K,2),'Vartype','unequal')
% p =
%      0
% ci =
%    -0.0128
%    -0.0123
% stat = 
%     tstat: -88.6506
%        df: 1.9812e+04
%        sd: [ 0.0105 0.0095 ]

%% plot each with their null

sp = tight_subplot(1,2,.12,[.05 .05],[.1 .1]);
axes(sp(1))
histogram(mean(eval_wsbm_K,2),'normalization','probability','FaceColor',cmap(1,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
hold
histogram(mean(eval_wsbmRand_K,2),'normalization','probability','FaceColor',cmap(1,:),'FaceAlpha',0.25,'EdgeAlpha',0.05)
pbaspect([1 0.8 1])
% [~,p,ci,stat] = ttest2(mean(eval_wsbmRand_K,2),mean(eval_wsbm_K,2),'Vartype','unequal')
% p =
% 
%      0
% ci =
% 
%     0.0101
%     0.0109
% stat = 
%     tstat: 50.8078
%        df: 1.6193e+04
%        sd: [2x1 double]

axes(sp(2))
histogram(mean(eval_mod_K,2),'normalization','probability','FaceColor',cmap(2,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
hold
histogram(mean(eval_modRand_K,2),'normalization','probability','FaceColor',cmap(2,:),'FaceAlpha',0.25,'EdgeAlpha',0.05)
pbaspect([1 0.8 1])
% [~,p,ci,stat] = ttest2(mean(eval_modRand_K,2),mean(eval_mod_K,2),'Vartype','unequal')
% p =
%    1.3789e-72
% ci =
%    -0.0043
%    -0.0035
% stat = 
%     tstat: -18.1206
%        df: 1.4522e+04
%        sd: [2x1 double]

%% plot each statistic

measure_names = {'Degree' 'Clustering' 'Betweenness' 'Edge length'} ;

for idx = 1:4
   
    subplot(1,4,idx)
    histogram(eval_wsbm_K(:,idx),'BinMethod','sturges',...
        'normalization','probability','EdgeAlpha',0.05,'FaceAlpha',0.6)
    hold
    histogram(eval_mod_K(:,idx),'BinMethod','sturges',...
        'normalization','probability','EdgeAlpha',0.05,'FaceAlpha',0.6)
    axis square
    title(measure_names{idx},'FontWeight','normal')
   
    [~,p,ci,stat] = ttest2(eval_mod_K(:,idx),eval_wsbm_K(:,idx),'Vartype','unequal')
    
end

% p = 
%      0
% ci =
%    -0.0425
%    -0.0415
% stat = 
%     tstat: -157.2182
%        df: 1.9990e+04
%        sd: [2x1 double]
%  
% p =
%      0
% ci =
%     0.0226
%     0.0236
% stat = 
%     tstat: 83.6199
%        df: 1.9015e+04
%        sd: [2x1 double]
% 
% p =
%      0
% ci =
% 
%     0.0955
%     0.0967
% stat = 
%     tstat: 304.2667
%        df: 1.9883e+04
%        sd: [2x1 double]
% 
% p =
%      0
% ci =
%    -0.0272
%    -0.0267
% stat =
%     tstat: -202.8450
%        df: 1.9998e+04
%        sd: [2x1 double]

%% earth mover difference

for idx = 1:4
   
    subplot(1,4,idx)
    histogram(eval_wsbm_EMD(:,idx),'normalization','probability','EdgeAlpha',0.01)
    hold
    histogram(eval_mod_EMD(:,idx),'normalization','probability','EdgeAlpha',0.01)

end
legend('WSBM model','Modular model')

















