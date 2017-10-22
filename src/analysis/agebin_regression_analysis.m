%% lets analyze the age bin scripts

addpath('~/JOSHSTUFF/scripts/BCT/2017_01_15_BCT/')

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

% load the data we need to analyze this ish. 
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_fit_wsbm_script_v7p3.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_agebin_template_results_v7p3.mat');
load(loadName) ;

thresholds = 4 ;

ageBin_wsbmGatherStruct = cell([(thresholds+1) 1]) ;
ageBin_wsbmRegression = cell([(thresholds+1) 1]) ;

%% 

for ageBinIdx=1:(thresholds +1)

    % wsbm with ya template %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ca_wsbm = community_assign(consensusCrossIdx{ageBinIdx}) ;

    ageBin_wsbmGatherStruct{ageBinIdx} = ma_gather_data2(dataStruct, datasetDemo, ...
        ca_wsbm(:,2), selectNodesFrmRaw, 1, 1, 'countVolNormMat');

    fits = {'linear' 'quadratic' 'poisson' } ;
    %fits = {'linear' 'quadratic' } ;

    numBlocks = consensusCrossIdx{ageBinIdx}.R_Struct.k ;

    % these get passed into regression func
    otherArgs = {1 500 [] 10000 };

    % predictor
    X = datasetDemo.age ;
    % predicted
    % Y = rdens --> rdens in gather struct

    ageBin_wsbmRegression{ageBinIdx} = run_regression_over_yMat(...
        X,ageBin_wsbmGatherStruct{ageBinIdx}.rdens,fits,otherArgs) ;

end

%% make a data struct to save info...

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_agebin_template_regress_v7p3.mat');
save(outName,...
    'ageBin_wsbmRegression',...
    'ageBin_wsbmGatherStruct',...
    '-v7.3')



