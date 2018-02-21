clc
clearvars

%% load the necessary data

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_templateModel_1.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
load(loadName) ;

%% get the stats....
% and then write to a csv

nSubj = length(dataStruct);
nNodes = templateModel.Data.n ;
nComm = templateModel.R_Struct.k ; 

% get the template data
templateData = templateModel.Data.Raw_Data ;
templateData(isnan(templateData)) = 0;

commSizes = histcounts(sort(comVecs.wsbm));

%% gather the subject-level data

subjDataMat = zeros([ nNodes nNodes nSubj ]);
for idx = 1:nSubj  

    tmpAdj = dataStruct(idx).countVolNormMat(selectNodesFrmRaw, selectNodesFrmRaw);
    % get rid of the diagonal
    %n=size(tmpAdj,1);
    tmpAdj(1:nNodes+1:end) = 0; 
    % mask out AdjMat entries below mask_thr
    tmpAdj_mask = dataStruct(idx).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > MASK_THR ;    
    tmpAdj_mask(tmpAdj_mask > 0) = 1 ;   
    tmpAdj = tmpAdj .* tmpAdj_mask ;
    tmpAdj(isnan(tmpAdj)) = 0 ;
    subjDataMat(:,:,idx) = tmpAdj ;
end

%% gather some stats on the dataset

% template density
density_und(templateData > 0)

% binary total connection density
totBin = zeros([nSubj 1]) ;
for idx = 1:nSubj
    totBin(idx) = sum(sum(subjDataMat(:,:,idx)>0)) ;
end

% weighted total strength
totDen = zeros([nSubj 1]) ;
for idx = 1:nSubj
    totDen(idx) = sum(sum(subjDataMat(:,:,idx))) ;
end

% brain volume
brainVol = datasetDemo.SupraTentorial ;

% movement
mov = datasetDemo.movement ;












