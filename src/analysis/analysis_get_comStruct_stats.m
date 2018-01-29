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

%% look at where the high degree nodes are in the template data

highD_level = 75 ; 

% get certain percentile of the data
highD_thr = prctile(sum(templateData,2),highD_level) ;
highD_node = sum(templateData,2) > highD_thr ;

wsbm_highDappear = histc(comVecs.wsbm(highD_node),1:nComm) ;
mod_highDappear = histc(comVecs.mod(highD_node),1:nComm) ;
yeo_highDappear = histc(comVecs.yeo(highD_node),1:nComm) ;

%%

% get the subject-level data
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

%% template measures

% strength

% wsbm
[~,wsbm_avgWeiBm] =  get_block_mat(templateData,comVecs.wsbm) ;
wsbm_avgWei_within = wsbm_avgWeiBm(~~eye(nComm));
wsbm_avgWei_btwn = (sum(wsbm_avgWeiBm,2) - wsbm_avgWei_within) ./ (nComm-1) ;
% modular
[~,mod_avgWeiBm] =  get_block_mat(templateData,comVecs.mod) ;
mod_avgWei_within = mod_avgWeiBm(~~eye(nComm));
mod_avgWei_btwn = (sum(mod_avgWeiBm,2) - mod_avgWei_within) ./ (nComm-1) ;
% yeo
[~,yeo_avgWeiBm] =  get_block_mat(templateData,comVecs.yeo) ;
yeo_avgWei_within = yeo_avgWeiBm(~~eye(7));
yeo_avgWei_btwn = (sum(yeo_avgWeiBm,2) - yeo_avgWei_within) ./ (7-1) ;

% modularity (community contribution)
 
% wsbm
[~,wsbm_q] = eval_modularity_wu(templateData,comVecs.wsbm) ;
% modular
[~,mod_q] = eval_modularity_wu(templateData,comVecs.mod) ;
% yeo
[~,yeo_q] = eval_modularity_wu(templateData,comVecs.yeo) ;

% avg participation coef

wsbm_parti = participation_coef(templateData,comVecs.wsbm,0) ;
mod_parti = participation_coef(templateData,comVecs.mod,0) ;
yeo_parti = participation_coef(templateData,comVecs.yeo,0) ;

wsbm_com_parti = zeros([nComm 1]);
mod_com_parti = zeros([nComm 1]);
yeo_com_parti = zeros([nComm 1]);

for idx = 1:nComm
   
    wsbm_com_parti(idx) = mean(wsbm_parti(comVecs.wsbm == idx));
    mod_com_parti(idx) = mean(mod_parti(comVecs.mod == idx));
    yeo_com_parti(idx) = mean(yeo_parti(comVecs.yeo == idx));

end

% avg assortativity

wsbm_assort = eval_com_assortatvity_wu(templateData,comVecs.wsbm) ;
mod_assort = eval_com_assortatvity_wu(templateData,comVecs.mod) ;
yeo_assort = eval_com_assortatvity_wu(templateData,comVecs.yeo) ;

%% more analysis on the model

predict_e = templateModel.Para.predict_e ;
predict_w = templateModel.Para.predict_w ;

predictE_zscore = zscore(predict_e) ;
predictW_zscore = zscore(predict_w) ;

%% subject measures

wsbm_avgWei_within_subjAll = zeros([nComm nSubj]) ;
wsbm_avgWei_btwn_subjAll = zeros([nComm nSubj]) ;
mod_avgWei_within_subjAll = zeros([nComm nSubj]) ;
mod_avgWei_btwn_subjAll = zeros([nComm nSubj]) ;
yeo_avgWei_within_subjAll = zeros([7 nSubj]) ;
yeo_avgWei_btwn_subjAll = zeros([7 nSubj]) ;

wsbm_q_subjAll = zeros([nComm nSubj]) ;
mod_q_subjAll = zeros([nComm nSubj]) ;
yeo_q_subjAll = zeros([7 nSubj]) ;

wsbm_parti_subjAll = zeros([nComm nSubj]) ;
mod_parti_subjAll = zeros([nComm nSubj]);
yeo_parti_subjAll = zeros([7 nSubj]);

wsbm_assort_subjAll = zeros([nComm nSubj]) ;
mod_assort_subjAll = zeros([nComm nSubj]);
yeo_assort_subjAll = zeros([7 nSubj]);

for idx = 1:nSubj

    currSubjData = subjDataMat(:,:,idx);
    
    % strength 
    
    % wsbm
    [~,tmp_avgWeiBm] =  get_block_mat(currSubjData,comVecs.wsbm) ;
    wsbm_avgWei_within_subjAll(:,idx) = tmp_avgWeiBm(~~eye(nComm));
    wsbm_avgWei_btwn_subjAll(:,idx) = (sum(tmp_avgWeiBm,2) - wsbm_avgWei_within_subjAll(:,idx)) ./ (nComm-1) ;
    % modular
    [~,tmp_avgWeiBm] =  get_block_mat(currSubjData,comVecs.mod) ;
    mod_avgWei_within_subjAll(:,idx)  = tmp_avgWeiBm(~~eye(nComm));
    mod_avgWei_btwn_subjAll(:,idx)  = (sum(tmp_avgWeiBm,2) - mod_avgWei_within_subjAll(:,idx)) ./ (nComm-1) ;
    % yeo
    [~,tmp_avgWeiBm] =  get_block_mat(currSubjData,comVecs.yeo) ;
    yeo_avgWei_within_subjAll(:,idx)  = tmp_avgWeiBm(~~eye(7));
    yeo_avgWei_btwn_subjAll(:,idx)  = (sum(tmp_avgWeiBm,2) - yeo_avgWei_within_subjAll(:,idx)) ./ (7-1) ;
    
    % modularity
    
    % wsbm
    [~,~,wsbm_q_subjAll(:,idx)] = eval_modularity_wu(currSubjData,comVecs.wsbm) ;
    % modular
    [~,~,mod_q_subjAll(:,idx)] = eval_modularity_wu(currSubjData,comVecs.mod) ;
    % yeo
    [~,~,yeo_q_subjAll(:,idx)] = eval_modularity_wu(currSubjData,comVecs.yeo) ;

    % avg parti

    wsbm_partiTmp = participation_coef(currSubjData,comVecs.wsbm);
    mod_partiTmp = participation_coef(currSubjData,comVecs.mod);
    yeo_partiTmp = participation_coef(currSubjData,comVecs.yeo);
    
    for jdx = 1:nComm
   
        wsbm_parti_subjAll(jdx,idx) = mean(wsbm_partiTmp(comVecs.wsbm == jdx));
        mod_parti_subjAll(jdx,idx) = mean(mod_partiTmp(comVecs.mod == jdx));
        yeo_parti_subjAll(jdx,idx) = mean(yeo_partiTmp(comVecs.yeo == jdx));

    end
    
    % fix this size
    yeo_parti_subjAll = yeo_parti_subjAll(1:7,:);
    
    % assort
    
    wsbm_assort_subjAll(:,idx) = eval_com_assortatvity_wu(currSubjData,comVecs.wsbm);
    mod_assort_subjAll(:,idx) = eval_com_assortatvity_wu(currSubjData,comVecs.mod);  
    yeo_assort_subjAll(:,idx) = eval_com_assortatvity_wu(currSubjData,comVecs.yeo); 

end

%% subject-level where high degree nodes show up

wsbm_highDappear_subjAll = zeros([ nComm nSubj ]);
mod_highDappear_subjAll = zeros([ nComm nSubj ]);
yeo_highDappear_subjAll = zeros([ 7 nSubj ]);

highD_level = 75 ; 

for idx = 1:nSubj

    currSubjData = subjDataMat(:,:,idx) ; 
    
    % get certain percentile of the data
    highD_thr = prctile(sum(currSubjData > 0,2),highD_level) ;
    highD_node = sum(currSubjData > 0,2) > highD_thr ;

    wsbm_highDappear_subjAll(:,idx) = histc(comVecs.wsbm(highD_node),1:nComm) ;
    mod_highDappear_subjAll(:,idx) = histc(comVecs.mod(highD_node),1:nComm) ;
    yeo_highDappear_subjAll(:,idx) = histc(comVecs.yeo(highD_node),1:7) ;
    
end
  
% ICC
wsbm_highD_icc = IPN_icc(wsbm_highDappear_subjAll,3,'single') ;
mod_highD_icc = IPN_icc(mod_highDappear_subjAll,3,'single') ;
yeo_highD_icc = IPN_icc(yeo_highDappear_subjAll,3,'single') ;

tmpICC = cell([3 1]);
tmpData = cell([3 1]);
tmpData{1} = wsbm_highDappear_subjAll ;
tmpData{2} = mod_highDappear_subjAll ;
tmpData{3} = yeo_highDappear_subjAll ;

for rep = 1:3

    % lets see if we can bootstrap these values
    nBoot = 500 ;
    tmpRes = zeros([ nBoot 1 ]);

    % get bootstrp indicies
    [~,bootInd] = bootstrp(nBoot,@(a)[],1:nSubj) ;

    for idx = 1:nBoot

        tmpRes(idx) = IPN_icc(tmpData{rep}(:,bootInd(:,idx)),3,'single');
    end

    tmpICC{rep} = prctile(tmpRes,[2.5 97.5]) ;
    
end

wsbm_highD_icc_c95 = tmpICC{1};
mod_highD_icc_c95 = tmpICC{2};
yeo_highD_icc_c95 = tmpICC{3};

%% make a table

com_names = { '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' }' ;
yeo_names = { 'Vis' 'SomMot' 'DorsAttn' 'SalVent' 'Limbic' 'Cont' 'Default' } ; 

multipleTables = cell([3 1]);

for idx = 1:3
    
    %%%%%%%%%%%%%%%%%%
    % template stats %
    %%%%%%%%%%%%%%%%%%

    % wsbm
    multipleTables{idx} = table() ;
    multipleTables{idx}.mean_within_weight = wsbm_avgWei_within ;
    multipleTables{idx}.mean_btwn_weight = wsbm_avgWei_btwn ;
    multipleTables{idx}.q_prcnt = wsbm_q ./ sum(wsbm_q) ;
    multipleTables{idx}.mean_comm_parti_coef = wsbm_com_parti;
    multipleTables{idx}.assort_of_comm = wsbm_assort ;
    
    if idx == 3
        multipleTables{idx}.Properties.RowNames = com_names ;
    else
        multipleTables{idx}.Properties.RowNakes = yeo_names ;
    end
    
end

%%%%%%%%%%%%%%%%%%
% template stats %
%%%%%%%%%%%%%%%%%%

% wsbm
wsbm_tmpDataTbl = table() ;
wsbm_tmpDataTbl.mean_within_weight = wsbm_avgWei_within ;
wsbm_tmpDataTbl.mean_btwn_weight = wsbm_avgWei_btwn ;
wsbm_tmpDataTbl.q_prcnt = wsbm_q ./ sum(wsbm_q) ;
wsbm_tmpDataTbl.mean_comm_parti_coef = wsbm_com_parti;
wsbm_tmpDataTbl.assort_of_comm = wsbm_assort ;
wsbm_tmpDataTbl.Properties.RowNames = com_names ;

% modular
mod_tmpDataTbl = table() ;
mod_tmpDataTbl.mean_within_weight = mod_avgWei_within ;
mod_tmpDataTbl.mean_btwn_weight = mod_avgWei_btwn ;
mod_tmpDataTbl.q_prcnt = mod_q ./ sum(mod_q);
mod_tmpDataTbl.mean_comm_parti_coef = mod_com_parti;
mod_tmpDataTbl.assort_of_comm = mod_assort ;
mod_tmpDataTbl.Properties.RowNames = com_names ;

% yeo
yeo_tmpDataTbl = table() ;
yeo_tmpDataTbl.mean_within_weight = yeo_avgWei_within ;
yeo_tmpDataTbl.mean_btwn_weight = yeo_avgWei_btwn ;
yeo_tmpDataTbl.q_prcnt = yeo_q ./ sum(yeo_q);
yeo_tmpDataTbl.mean_comm_parti_coef = yeo_com_parti(1:7);
yeo_tmpDataTbl.assort_of_comm = yeo_assort ;
yeo_tmpDataTbl.Properties.RowNames = yeo_names ;

%%%%%%%%%%%%%%%%
% subject-wise %
%%%%%%%%%%%%%%%%

% wsbm
wsbm_subjDataTbl = table() ;
wsbm_subjDataTbl.mean_within_weight_mean = mean(wsbm_avgWei_within_subjAll,2) ;
wsbm_subjDataTbl.mean_within_weight_median = median(wsbm_avgWei_within_subjAll,2) ;
wsbm_subjDataTbl.mean_within_weight_std = std(wsbm_avgWei_within_subjAll,[],2) ;

wsbm_subjDataTbl.mean_btwn_weight_mean = mean(wsbm_avgWei_btwn_subjAll,2) ;
wsbm_subjDataTbl.mean_btwn_weight_median = median(wsbm_avgWei_btwn_subjAll,2) ;
wsbm_subjDataTbl.mean_btwn_weight_std = std(wsbm_avgWei_btwn_subjAll,[],2) ;

wsbm_subjDataTbl.q_prcnt_mean = mean(wsbm_q_subjAll,2) ;
wsbm_subjDataTbl.q_prcnt_median = median(wsbm_q_subjAll,2) ;
wsbm_subjDataTbl.q_prcnt_std = std(wsbm_q_subjAll,[],2) ;

wsbm_subjDataTbl.mean_comm_parti_coef_mean = mean(wsbm_parti_subjAll,2);
wsbm_subjDataTbl.mean_comm_parti_coef_median = median(wsbm_parti_subjAll,2);
wsbm_subjDataTbl.mean_comm_parti_coef_std = std(wsbm_parti_subjAll,[],2);

wsbm_subjDataTbl.assort_of_comm_mean = mean(wsbm_assort_subjAll,2) ;
wsbm_subjDataTbl.assort_of_comm_median = median(wsbm_assort_subjAll,2) ;
wsbm_subjDataTbl.assort_of_comm_std = std(wsbm_assort_subjAll,[],2) ;

wsbm_subjDataTbl.Properties.RowNames = com_names ;

% modular
mod_subjDataTbl = table() ;
mod_subjDataTbl.mean_within_weight_mean = mean(mod_avgWei_within_subjAll,2) ;
mod_subjDataTbl.mean_within_weight_median = median(mod_avgWei_within_subjAll,2) ;
mod_subjDataTbl.mean_within_weight_std = std(mod_avgWei_within_subjAll,[],2) ;

mod_subjDataTbl.mean_btwn_weight_mean = mean(mod_avgWei_btwn_subjAll,2) ;
mod_subjDataTbl.mean_btwn_weight_median = median(mod_avgWei_btwn_subjAll,2) ;
mod_subjDataTbl.mean_btwn_weight_std = std(mod_avgWei_btwn_subjAll,[],2) ;

mod_subjDataTbl.q_prcnt_mean = mean(mod_q_subjAll,2) ;
mod_subjDataTbl.q_prcnt_median = median(mod_q_subjAll,2) ;
mod_subjDataTbl.q_prcnt_std = std(mod_q_subjAll,[],2) ;

mod_subjDataTbl.mean_comm_parti_coef_mean = mean(mod_parti_subjAll,2);
mod_subjDataTbl.mean_comm_parti_coef_median = median(mod_parti_subjAll,2);
mod_subjDataTbl.mean_comm_parti_coef_std = std(mod_parti_subjAll,[],2);

mod_subjDataTbl.assort_of_comm_mean = mean(mod_assort_subjAll,2) ;
mod_subjDataTbl.assort_of_comm_median = median(mod_assort_subjAll,2) ;
mod_subjDataTbl.assort_of_comm_std = std(mod_assort_subjAll,[],2) ;

mod_subjDataTbl.Properties.RowNames = com_names ;

% yeo
yeo_subjDataTbl = table() ;
yeo_subjDataTbl.mean_within_weight_mean = mean(yeo_avgWei_within_subjAll,2) ;
yeo_subjDataTbl.mean_within_weight_median = median(yeo_avgWei_within_subjAll,2) ;
yeo_subjDataTbl.mean_within_weight_std = std(yeo_avgWei_within_subjAll,[],2) ;

yeo_subjDataTbl.mean_btwn_weight_mean = mean(yeo_avgWei_btwn_subjAll,2) ;
yeo_subjDataTbl.mean_btwn_weight_median = median(yeo_avgWei_btwn_subjAll,2) ;
yeo_subjDataTbl.mean_btwn_weight_std = std(yeo_avgWei_btwn_subjAll,[],2) ;

yeo_subjDataTbl.q_prcnt_mean = mean(yeo_q_subjAll,2) ;
yeo_subjDataTbl.q_prcnt_median = median(yeo_q_subjAll,2) ;
yeo_subjDataTbl.q_prcnt_std = std(yeo_q_subjAll,[],2) ;

yeo_subjDataTbl.mean_comm_parti_coef_mean = mean(yeo_parti_subjAll,2);
yeo_subjDataTbl.mean_comm_parti_coef_median = median(yeo_parti_subjAll,2);
yeo_subjDataTbl.mean_comm_parti_coef_std = std(yeo_parti_subjAll,[],2);

yeo_subjDataTbl.assort_of_comm_mean = mean(yeo_assort_subjAll,2) ;
yeo_subjDataTbl.assort_of_comm_median = median(yeo_assort_subjAll,2) ;
yeo_subjDataTbl.assort_of_comm_std = std(yeo_assort_subjAll,[],2) ;

yeo_subjDataTbl.Properties.RowNames = yeo_names ;

%% write all of the tables to csv

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_wsbm_tmpl_stats.csv');
writetable(wsbm_tmpDataTbl,fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_mod_tmpl_stats.csv');
writetable(mod_tmpDataTbl,fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_yeo_tmpl_stats.csv');
writetable(yeo_tmpDataTbl,fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_wsbm_subjAgg_stats.csv');
writetable(wsbm_subjDataTbl,fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_mod_subjAgg_stats.csv');
writetable(mod_subjDataTbl,fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_yeo_subjAgg_stats.csv');
writetable(yeo_subjDataTbl,fileName,'WriteRowNames',true)











