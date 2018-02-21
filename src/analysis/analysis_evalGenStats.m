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

FIGURE_NAME = 'figC' ;

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

% % empircal edge existence and edge weights
% [~,modModel_w,~,modModel_e] = get_block_mat(templateAdj,comVecs.mod);
% 
% modModel_w = nonzeros(triu(modModel_w)) ;
% modModel_e = nonzeros(triu(modModel_e)) ;
% 
% modularityModel.Para.predict_w = modModel_w ;
% modularityModel.Para.predict_e = modModel_e ;

%% try out the evalWSBM code

% nNodes = templateModel.Data.n ;
% 
% templateSubj_data = dataStruct(datasetDemo.age > 25 & datasetDemo.age <= 35) ;
% [a,b,avgTemp_dist] = make_template_mat(templateSubj_data, ...
%     LEFT_HEMI_NODES, ...
%     RIGHT_HEMI_NODES, ...
%     MASK_THR_INIT) ; 
% 
% % actually replace the 0's with NaN
% %avgTemp(avgTemp == 0) = NaN ;
% avgTemp_dist = avgTemp_dist(selectNodesFrmRaw,selectNodesFrmRaw);
% % clear diagonal
% avgTemp_dist(1:nNodes+1:end)=0; 

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

% add yeo
[eval_yeo_B,eval_yeo_E,eval_yeo_K,eval_yeo_EMD] = eval_genWsbm_model1(yeoModel,avgTemp_dist,10000,0);
[eval_yeoRand_B,eval_yeoRand_E,eval_yeoRand_K] = eval_genWsbm_model1(yeoModel,avgTemp_dist,10000,1);


% load(strcat(PROJECT_DIR,'/data/processed/',OUTPUT_STR,'_evalGenReps.mat'))
save(strcat(PROJECT_DIR,'/data/processed/',OUTPUT_STR,'_evalGenReps.mat'),...
    'eval_wsbm_K','eval_wsbmRand_K','eval_wsbm_EMD',...
    'eval_mod_K','eval_modRand_K','eval_mod_EMD' ...
     )

%% plot energy 

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.4, 0.8]);

cmap = [0    0.4470    0.7410 ;
    0.8500    0.3250    0.0980];

histogram(mean(eval_wsbm_K,2),...
    'normalization','probability',...
    'FaceColor',cmap(1,:),'EdgeAlpha',0.01) 
hold 
histogram(mean(eval_mod_K,2),...
    'normalization','probability',...
    'FaceColor',cmap(2,:),'EdgeAlpha',0.01)
lg = legend('WSBM','Modular') ;
lg.FontSize = 16 ;

ylim([0 0.1])

axis square

xlabel('Mean K-S energy')
ylabel('Normalized frequency')

% [~,p,ci,stat] = ttest2(mean(eval_wsbm_K,2),mean(eval_mod_K,2))
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

ttest_struct = struct();
kstest_struct = struct();

[~,ttest_struct.p,ttest_struct.ci,ttest_struct.stat] = ttest2(mean(eval_wsbm_K,2),mean(eval_mod_K,2),'Vartype','unequal');
[~,kstest_struct.p,kstest_struct.stat] = kstest2(mean(eval_wsbm_K,2),mean(eval_mod_K,2));

% save it
fileName = 'mean_KS_energy.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

save(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/','mean_KS_energy_stats.mat'),...
    'ttest_struct','kstest_struct')

%% plot each with their null

sp = tight_subplot(1,2,.12,[.05 .05],[.1 .1]);
axes(sp(1))

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.8, 0.6]);

histogram(mean(eval_wsbm_K,2),'normalization','probability','FaceColor',cmap(1,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
hold
histogram(mean(eval_wsbmRand_K,2),'normalization','probability','FaceColor',cmap(1,:),'FaceAlpha',0.25,'EdgeAlpha',0.05)

ylim([0 0.1])

pbaspect([1 0.8 1])

xlabel('Mean K-S energy')
ylabel('Normalized frequency')
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

ylim([0 0.1])

pbaspect([1 0.8 1])

ylabel('Normalized frequency')
xlabel('Mean K-S energy')

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

% save it
fileName = 'mean_KS_energy_vs_null.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%% plot each statistic

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.8, 0.6]);

measure_names = {'Degree' 'Clustering' 'Betweenness' 'Edge length'} ;
measure_short_names = { 'd' 'c' 'b' 'e' } ;

for idx = 1:4
   
    subplot(1,4,idx)
    histogram(eval_wsbm_K(:,idx),'BinMethod','sturges',...
        'normalization','probability','EdgeAlpha',0.05,'FaceAlpha',0.6)
    hold
    histogram(eval_mod_K(:,idx),'BinMethod','sturges',...
        'normalization','probability','EdgeAlpha',0.05,'FaceAlpha',0.6)
    axis square
    title(measure_names{idx},'FontWeight','normal')
   
    if idx == 1
       ylabel('Normalized frequency') 
    end
    
    ylim([0 0.35])
    
    [~,p,ci,stat] = ttest2(eval_mod_K(:,idx),eval_wsbm_K(:,idx),'Vartype','unequal')
    
    xlabel(strcat('KS({\it ',measure_short_names{idx},'})'))
    
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

% save it
fileName = 'KS_of_stats.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%% earth mover difference

for idx = 1:4
   
    subplot(1,4,idx)
    histogram(eval_wsbm_EMD(:,idx),'normalization','probability','EdgeAlpha',0.01)
    hold
    histogram(eval_mod_EMD(:,idx),'normalization','probability','EdgeAlpha',0.01)

end
legend('WSBM model','Modular model')

















