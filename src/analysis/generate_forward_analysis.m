%% lets analyze the age bin scripts

addpath('~/JOSHSTUFF/scripts/BCT/2017_01_15_BCT/')

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

% load the data we need to analyze this ish. 
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_fit_wsbm_script_v7p3.mat');
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

numPerms = 1000 ;

nNodes = templateModel.Data.n ;

% save some results
% permsQ = zeros([ numPerms 1 ]);
% permsCCoef = zeros([ numPerms 1 ]);
permsWSBMAssort = zeros([ numPerms 1 ]);
permsMODAssort = zeros([ numPerms 1 ]); 

for idx=1:numPerms
   
    disp(idx)
    
    [~,tmpAdj] = genAdj_wsbm(templateModel);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    
    %permsQ(idx) = eval_modularity_wu(tmpAdj) ;
    %permsCCoef(idx) = clustering_coef_wu(tmpAdj) ;
    permsWSBMAssort(idx) = assortativity_wei(tmpAdj,0) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,tmpAdj] = genAdj_wsbm(modularityModel);
    tmpAdj(1:nNodes+1:end)=0; %clear diagonal
    
    permsMODAssort(idx) = assortativity_wei(tmpAdj,0) ;
    
end




