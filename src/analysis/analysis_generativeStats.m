%% clear stuff

clc
clearvars

%% load the necessary data

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

% % empircal edge existence and edge weights
% [~,modModel_w,~,modModel_e] = get_block_mat(templateAdj,comVecs.mod);
% 
% modModel_w = nonzeros(triu(modModel_w)) ;
% modModel_e = nonzeros(triu(modModel_e)) ;
% 
% modularityModel.Para.predict_w = modModel_w ;
% modularityModel.Para.predict_e = modModel_e ;

%% setup vars for perm tests
% select some metrics that can be reduce to scalar, to get distributions
% across many permutations

numPerms = 1000 ;
nNodes = templateModel.Data.n ;

wsbmPerm = struct();
modPerm = struct();
randPerm = struct(); 

% save some results
wsbmPerm.q = zeros([ numPerms 1 ]);
wsbmPerm.assort = zeros([ numPerms 1 ]);
wsbmPerm.parti = zeros([ numPerms 1 ]);
wsbmPerm.gl_eff = zeros([ numPerms 1]);
wsbmPerm.di_eff = zeros([ numPerms 1]);
wsbmPerm.trans = zeros([ numPerms 1]);

modPerm.q = zeros([ numPerms 1 ]);
modPerm.assort = zeros([ numPerms 1 ]);
modPerm.parti = zeros([ numPerms 1 ]); 
modPerm.gl_eff = zeros([ numPerms 1 ]); 
modPerm.di_eff = zeros([ numPerms 1]);
modPerm.trans = zeros([ numPerms 1 ]);

randPerm.q = zeros([ numPerms 1 ]);
randPerm.assort = zeros([ numPerms 1 ]);
randPerm.parti = zeros([ numPerms 1 ]); 
randPerm.gl_eff = zeros([ numPerms 1 ]);
randPerm.di_eff = zeros([ numPerms 1]);
randPerm.trans = zeros([ numPerms 1]);

%% iterate through permuations

for idx=1:numPerms
   
    disp(idx)
    
    [~,tmpAdj] = genAdj_wsbm(templateModel);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj = tmpAdj ~= 0;
    tmpAdj = double(tmpAdj);
    
    wsbmPerm.q(idx) = eval_modularity_wu(tmpAdj,comVecs.wsbm) ;
    wsbmPerm.assort(idx) = assortativity_bin(tmpAdj,0) ;
    wsbmPerm.parti(idx) = median(participation_coef(tmpAdj,comVecs.wsbm));
    wsbmPerm.gl_eff(idx) = efficiency_bin(tmpAdj);
    wsbmPerm.di_eff(idx) = diffusion_efficiency(tmpAdj);
    wsbmPerm.trans(idx) = transitivity_bu(tmpAdj);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,tmpAdj] = genAdj_wsbm(modularityModel);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj = tmpAdj ~= 0;
    tmpAdj = double(tmpAdj);

    modPerm.q(idx) = eval_modularity_wu(tmpAdj,comVecs.wsbm) ;
    modPerm.assort(idx) = assortativity_bin(tmpAdj,0) ;
    modPerm.parti(idx) = median(participation_coef(tmpAdj,comVecs.wsbm));
    modPerm.gl_eff(idx) = efficiency_bin(tmpAdj);
    modPerm.di_eff(idx) = diffusion_efficiency(tmpAdj);
    modPerm.trans(idx) = transitivity_bu(tmpAdj);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % randomize template
    randModel = randomize_wsbm_para(templateModel,3);
    [~,tmpAdj] = genAdj_wsbm(randModel);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj = tmpAdj ~= 0;
    tmpAdj = double(tmpAdj);

    randPerm.q(idx) = eval_modularity_wu(tmpAdj,comVecs.wsbm) ;
    randPerm.assort(idx) = assortativity_bin(tmpAdj,0) ;
    randPerm.parti(idx) = median(participation_coef(tmpAdj,comVecs.wsbm));
    randPerm.gl_eff(idx) = efficiency_bin(tmpAdj);
    randPerm.di_eff(idx) = diffusion_efficiency(tmpAdj);
    randPerm.trans(idx) = transitivity_bu(tmpAdj);

end

%% empirical yo

tempBinData = templateAdj;
tempBinData = tempBinData ~= 0 ;
tempBinData = double(tempBinData);

empRes = struct();
empRes.q = eval_modularity_wu(tempBinData,comVecs.yeo) ;
empRes.assort = assortativity_bin(tempBinData,0) ;
empRes.parti = median(participation_coef(tempBinData,comVecs.wsbm));
empRes.gl_eff = efficiency_bin(tempBinData);
empRes.di_eff = diffusion_efficiency(tempBinData);
empRes.trans = transitivity_bu(tempBinData);

%% what about these measures in all of the subjs

nSubj = length(dataStruct);
nNodes = templateModel.Data.n ;
subjDataMat = zeros([ nNodes nNodes nSubj ]);

for idx = 1:nSubj

    tmpAdj = dataStruct(idx).countVolNormMat(selectNodesFrmRaw, selectNodesFrmRaw);
    % get rid of the diagonal
    %n=size(tmpAdj,1);
    tmpAdj(1:nNodes+1:end) = 0; 
    % mask out AdjMat entries below mask_thr
    tmpAdj_mask = dataStruct(idx).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > 1 ;    
    tmpAdj_mask(tmpAdj_mask > 0) = 1 ;   
    tmpAdj = tmpAdj .* tmpAdj_mask ;
    
    subjDataMat(:,:,idx) = tmpAdj ;
    
end

empSubjs = struct;
empSubjs.q = zeros([ nSubj 1 ]);
empSubjs.assort = zeros([ nSubj 1 ]);
empSubjs.parti = zeros([ nSubj 1 ]); 
empSubjs.gl_eff = zeros([ nSubj 1 ]);
empSubjs.di_eff = zeros([ nSubj 1]);
empSubjs.trans = zeros([ nSubj 1]);

for idx = 1:nSubj
   
    tmpAdj = subjDataMat(:,:,idx) ;
    
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    tmpAdj = tmpAdj ~= 0;
    tmpAdj = double(tmpAdj);
    
    %empSubjs.q(idx) = eval_modularity_wu(tmpAdj,comVecs.wsbm) ;
    empSubjs.assort(idx) = assortativity_bin(tmpAdj,0) ;
    %empSubjs.parti(idx) = median(participation_coef(tmpAdj,comVecs.wsbm));
    empSubjs.gl_eff(idx) = efficiency_bin(tmpAdj);
    %empSubjs.di_eff(idx) = diffusion_efficiency(tmpAdj);
    empSubjs.trans(idx) = transitivity_bu(tmpAdj);
   
end

%%

%suptitle('wsbm vs rand wsbm')

% measures = { 'q' 'assort' 'parti' 'gl_eff' 'di_eff' 'trans' };
% measure_names = {'Q' 'Assort.' 'Parti. Coef.' ...
%     'Global Eff.' 'Diffusion Eff.' 'Transitivity'};

% only non-community measures 
measures = { 'assort' 'gl_eff' 'trans' };
measure_names = { 'Assort.' 'Global Eff.' 'Transitivity'};

for idx = 1:3
   
    subplot(2,2,idx)
    histogram(wsbmPerm.(measures{idx}))
    hold
    histogram(modPerm.(measures{idx}))
    histogram(randPerm.(measures{idx}))
    histogram(empSubjs.(measures{idx}))
    % and plot a line for the empirical
    ylimits = ylim ;
    plot([empRes.(measures{idx}) empRes.(measures{idx})],...
        [ylimits(1) ylimits(2)*0.95],...
        'Color',[0 1 0 0.9],'LineWidth',1)
    title(measure_names{idx})
    
end

%legend(subplot(2,2,idx),{'wsbm' 'mod.' 'rand.'})

%%

subplot(1,4,1) 
histogram(permsWSBM_Q)
hold
histogram(permsRAND_Q)
histogram(permsMOD_Q)
% and plot a line for the empirical
ylimits = ylim ;
plot([emp_Q emp_Q],[ylimits(1) ylimits(2)*0.95],'Color',[0 1 0 0.9],'LineWidth',1)
title('Q')

subplot(1,4,2) 
histogram(permsWSBM_CCoef)
hold
histogram(permsRAND_CCoef)
histogram(permsMOD_CCoef)
% and plot a line for the empirical
ylimits = ylim ;
plot([emp_CCoef emp_CCoef],[ylimits(1) ylimits(2)*0.95],'Color',[0 1 0 0.9],'LineWidth',1)
title('CCoef')

subplot(1,4,3) 
histogram(permsWSBM_Assort)
hold
histogram(permsRAND_Assort)
histogram(permsMOD_Assort)
% and plot a line for the empirical
ylimits = ylim ;
plot([emp_Assort emp_Assort],[ylimits(1) ylimits(2)*0.95],'Color',[0 1 0 0.9],'LineWidth',1)
title('Assort')

subplot(1,4,4) 
histogram(permsWSBM_Parti)
hold
histogram(permsRAND_Parti)
histogram(permsMOD_Parti)
% and plot a line for the empirical
ylimits = ylim ;
plot([emp_Parti emp_Parti],[ylimits(1) ylimits(2)*0.95],'Color',[0 1 0 0.9],'LineWidth',1)
title('Parti')

%%

subplot(1,2,1) 
histogram(permsWSBM_Eff)
hold
histogram(permsRAND_Eff)
histogram(permsMOD_Eff)
% and plot a line for the empirical
ylimits = ylim ;
plot([emp_Eff emp_Eff],[ylimits(1) ylimits(2)*0.95],'Color',[0 1 0 0.9],'LineWidth',1)
title('Eff')

subplot(1,2,2) 
histogram(permsWSBM_DifEff)
hold
histogram(permsRAND_DifEff)
histogram(permsMOD_DifEff)
% and plot a line for the empirical
ylimits = ylim ;
plot([emp_DifEff emp_DifEff],[ylimits(1) ylimits(2)*0.95],'Color',[0 1 0 0.9],'LineWidth',1)
title('Eff')

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





