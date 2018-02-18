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

modelNames = { 'wsbm' 'mod' 'yeo' } ;

% strength
tmpl_strength = struct() ;
tmpl_mod = struct() ; 
tmpl_parti = struct() ;
tmpl_assort = struct() ;

for mN = 1:length(modelNames)

    if mN == 3
        nComm = 7 ;
    else
        nComm = 10 ;
    end

    % wsbm
    [~,tmp_avgWeiBm] =  get_block_mat(templateData,comVecs.(modelNames{mN})) ;
    tmpl_strength.(modelNames{mN}).within = tmp_avgWeiBm(~~eye(nComm));
    tmpl_strength.(modelNames{mN}).between = (sum(tmp_avgWeiBm,2) - ...
        tmpl_strength.(modelNames{mN}).within) ./ (nComm-1) ;

    % modularity (community contribution)
    [~,tmpl_mod.(modelNames{mN})] = eval_modularity_wu(templateData,comVecs.(modelNames{mN})) ;

    % avg participation coef
    tmpl_parti.(modelNames{mN}) = participation_coef(templateData,comVecs.(modelNames{mN}),0) ;

    % avg coom parti coef
    tmpl_parti.comm.(modelNames{mN}) = zeros([nComm 1]);
    for jdx = 1:nComm
        tmpl_parti.comm.(modelNames{mN})(jdx) = mean(tmpl_parti.(modelNames{mN})(comVecs.(modelNames{mN}) == jdx));
    end

    % avg assortativity
    tmpl_assort.(modelNames{mN}) = eval_com_assortatvity_wu(templateData,comVecs.(modelNames{mN})) ;

end

%% more analysis on the model

predict_e = templateModel.Para.predict_e ;
predict_w = templateModel.Para.predict_w ;

predictE_zscore = zscore(predict_e) ;
predictW_zscore = zscore(predict_w) ;

%% subject measures

modelNames = { 'wsbm' 'mod' 'yeo' } ;

% strength
subj_strength = struct() ;
subj_mod = struct() ; 
subj_parti = struct() ;
subj_assort = struct() ;

for mN = 1:length(modelNames)

    if mN == 3
        nComm = 7 ;
    else
        nComm = 10 ;
    end

    for idx = 1:nSubj

        currSubjData = subjDataMat(:,:,idx);

        % strength 
%         subj_strength.(modelNames{mN}).within = zeros([nComm nSubj]) ;
%         subj_strength.(modelNames{mN}).between = zeros([nComm nSubj]) ;

        [~,tmp_avgWeiBm] =  get_block_mat(currSubjData,comVecs.(modelNames{mN})) ;
        subj_strength.(modelNames{mN}).within(:,idx) = tmp_avgWeiBm(~~eye(nComm));
        subj_strength.(modelNames{mN}).between(:,idx) = (sum(tmp_avgWeiBm,2) - ...
            subj_strength.(modelNames{mN}).within(:,idx)) ./ (nComm-1) ;

        % modularity (community contribution)
        [~,~,subj_mod.(modelNames{mN})(:,idx)] = eval_modularity_wu(currSubjData,comVecs.(modelNames{mN})) ;

        % avg participation coef
        subj_parti.(modelNames{mN})(:,idx) = participation_coef(currSubjData,comVecs.(modelNames{mN}),0) ; 

        %subj_parti.comm.(modNames{mN})
        for jdx = 1:nComm

            subj_parti.comm.(modelNames{mN})(jdx,idx) = mean( subj_parti.(modelNames{mN})(comVecs.(modelNames{mN}) == jdx,idx) ) ;
        end

        % avg assortativity
        subj_assort.(modelNames{mN})(:,idx) = eval_com_assortatvity_wu(currSubjData,comVecs.(modelNames{mN})) ;

    end

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

modelNames = { 'wsbm' 'mod' 'yeo' } ;
com_names = { '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' }' ;
yeo_names = { 'Vis' 'SomMot' 'DorsAttn' 'SalVent' 'Limbic' 'Cont' 'Default' }' ; 

colNames = { 'Mean_within_strength' ...
    'Mean_between_comm_strength' 'Q_contribution_fraction'...
    'Mean_comm_parti_coeff'    'Comm_assort' } ;

colNames2 = { ...
    'Mean_within_strength' 'Median_within_strength' 'StdDev_within_strength'...
    'Mean_between_comm_strength' 'Median_between_comm_strength' 'StdDev_between_comm_strength'...
    'Mean_Q_contribution_fraction' 'Median_Q_contribution_fraction' 'Stddev_Q_contribution_fraction'...
    'Mean_comm_parti_coeff' 'Median_comm_parti_coeff' 'StdDev_comm_parti_coeff'...
    'Mean_Comm_assort' 'Median_Comm_assort' 'StdDev_Comm_assort' } ;

multipleTablesTemplate = cell([3 1]);
multipleTablesSubj = cell([3 1]);

for mN = 1:length(modelNames)
    
    %%%%%%%%%%%%%%%%%%
    % template stats %
    %%%%%%%%%%%%%%%%%%

    roundTo = 2 ;
    
    % wsbm
    multipleTablesTemplate{mN} = table() ;
    multipleTablesTemplate{mN}.mean_within_weight = round(tmpl_strength.(modelNames{mN}).within,roundTo,'significant') ;
    multipleTablesTemplate{mN}.mean_btwn_weight = round(tmpl_strength.(modelNames{mN}).between,roundTo,'significant') ;
    multipleTablesTemplate{mN}.q_prcnt = round(tmpl_mod.(modelNames{mN})./ ...
        sum(tmpl_mod.(modelNames{mN})),roundTo,'significant') ;
    multipleTablesTemplate{mN}.mean_comm_parti_coef = round(tmpl_parti.comm.(modelNames{mN}),roundTo,'significant');
    multipleTablesTemplate{mN}.assort_of_comm = round(tmpl_assort.(modelNames{mN}),roundTo,'significant') ;
    
    if mN == 3
        multipleTablesTemplate{mN}.Properties.RowNames = yeo_names ;
    else
        multipleTablesTemplate{mN}.Properties.RowNames = com_names ;
    end
    
    multipleTablesTemplate{mN}.Properties.VariableNames = colNames;
    
    %%%%%%%%%%%%%%%%
    % subject-wise %
    %%%%%%%%%%%%%%%%
    
    multipleTablesSubj{mN} = table() ;
    multipleTablesSubj{mN}.mean_within_weight_mean = round(mean(subj_strength.(modelNames{mN}).within,2),roundTo,'significant') ;
    multipleTablesSubj{mN}.mean_within_weight_median = round(median(subj_strength.(modelNames{mN}).within,2),roundTo,'significant') ;
    multipleTablesSubj{mN}.mean_within_weight_std = round(std(subj_strength.(modelNames{mN}).within,[],2),roundTo,'significant') ;

    multipleTablesSubj{mN}.mean_btwn_weight_mean = round(mean(subj_strength.(modelNames{mN}).between,2),roundTo,'significant') ;
    multipleTablesSubj{mN}.mean_btwn_weight_median = round(median(subj_strength.(modelNames{mN}).between,2),roundTo,'significant') ;
    multipleTablesSubj{mN}.mean_btwn_weight_std = round(std(subj_strength.(modelNames{mN}).between,[],2),roundTo,'significant') ;

    multipleTablesSubj{mN}.q_prcnt_mean = round(mean(subj_mod.(modelNames{mN}),2),roundTo,'significant') ;
    multipleTablesSubj{mN}.q_prcnt_median = round(median(subj_mod.(modelNames{mN}),2),roundTo,'significant') ;
    multipleTablesSubj{mN}.q_prcnt_std = round(std(subj_mod.(modelNames{mN}),[],2),roundTo,'significant') ;

    multipleTablesSubj{mN}.mean_comm_parti_coef_mean = round(mean(subj_parti.comm.(modelNames{mN}),2),roundTo,'significant');
    multipleTablesSubj{mN}.mean_comm_parti_coef_median = round(median(subj_parti.comm.(modelNames{mN}),2),roundTo,'significant');
    multipleTablesSubj{mN}.mean_comm_parti_coef_std = round(std(subj_parti.comm.(modelNames{mN}),[],2),roundTo,'significant');

    multipleTablesSubj{mN}.assort_of_comm_mean = round(mean(subj_assort.(modelNames{mN}),2),roundTo,'significant') ;
    multipleTablesSubj{mN}.assort_of_comm_median = round(median(subj_assort.(modelNames{mN}),2),roundTo,'significant') ;
    multipleTablesSubj{mN}.assort_of_comm_std = round(std(subj_assort.(modelNames{mN}),[],2),roundTo,'significant') ;

    if mN == 3
        multipleTablesSubj{mN}.Properties.RowNames = yeo_names ;
    else
        multipleTablesSubj{mN}.Properties.RowNames = com_names ;
    end
    
    multipleTablesSubj{mN}.Properties.VariableNames = colNames2;  
    
end

%% write all of the tables to csv

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_wsbm_tmpl_stats.csv');
writetable(multipleTablesTemplate{1},fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_mod_tmpl_stats.csv');
writetable(multipleTablesTemplate{2},fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_yeo_tmpl_stats.csv');
writetable(multipleTablesTemplate{3},fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_wsbm_subjAgg_stats.csv');
writetable(multipleTablesSubj{1},fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_mod_subjAgg_stats.csv');
writetable(multipleTablesSubj{2},fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_yeo_subjAgg_stats.csv');
writetable(multipleTablesSubj{3},fileName,'WriteRowNames',true)











