clc
clearvars

%% load the necessary data

config_file='config_scale125.m';
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

% %% look at where the high degree nodes are in the template data
% 
% highD_level = 75 ; 
% 
% % get certain percentile of the data
% highD_thr = prctile(sum(templateData,2),highD_level) ;
% highD_node = sum(templateData,2) > highD_thr ;
% 
% wsbm_highDappear = histc(comVecs.wsbm(highD_node),1:nComm) ;
% mod_highDappear = histc(comVecs.mod(highD_node),1:nComm) ;
% yeo_highDappear = histc(comVecs.yeo(highD_node),1:nComm) ;

%% print info about the nodes in each community vect

% labNames = load('data/raw/NKIen1/yeo/nodeLabels.mat') ;
% labNames = labNames.nodeLabels ;
% labNames = regexprep(labNames,'.*17Networks_','') ;
% labNames = regexprep(labNames,'.label','') ;
% 
% % wsbm
% [wsbm_sort,wsbm_sort_idx] = sort(comVecs.wsbm);
% comm_tbl = table();
% comm_tbl.WSBM_community = wsbm_sort ;
% comm_tbl.WSBM_node_name =  labNames(wsbm_sort_idx) ;
% 
% % mod
% [mod_sort,mod_sort_idx] = sort(comVecs.mod);
% comm_tbl.Modular_community = mod_sort ;
% comm_tbl.Modular_node_name =  labNames(mod_sort_idx)  ;
% 
% % yeo
% % % realign to match figs
% % tmp_yeo = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.yeo) ;
% % [yeo_sort,yeo_sort_idx] = sort(tmp_yeo);
% % comm_tbl.Yeo_community = yeo_sort ;
% % comm_tbl.Yeo_node_name =  labNames(yeo_sort_idx) ;
% 
% fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_node_comm.csv');
% writetable(comm_tbl,fileName,'WriteRowNames',true)

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

%% template measures

% modelNames = { 'wsbm' 'mod' 'yeo' } ;
modelNames = { 'wsbm' 'mod' } ;

% strength
tmpl_strength = struct() ;
tmpl_mod = struct() ; 
tmpl_parti = struct() ;
tmpl_assort = struct() ;

for mN = 1:length(modelNames)

    if mN == 3
        nComm = 7 ;
    else
        nComm = templateModel.R_Struct.k ;
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

% find interaction low edge density, high edge weight
tmp = (predictE_zscore < -1 ) & (predictW_zscore > 1) ;
tmp = make_square(tmp) ;

% high edge and weight
tmp = (predictE_zscore > 1 ) & (predictW_zscore > 2.5) ;
tmp = make_square(tmp) ;

%% test to see which interactions are more common than null model


%% subject measures

% modelNames = { 'wsbm' 'mod' 'yeo' } ;
modelNames = { 'wsbm' 'mod' } ;

% strength
subj_strength = struct() ;
subj_mod = struct() ; 
subj_parti = struct() ;
subj_assort = struct() ;

for mN = 1:length(modelNames)

    if mN == 3
        nComm = 7 ;
    else
        nComm = templateModel.R_Struct.k ;
    end

    for idx = 1:nSubj

        currSubjData = subjDataMat(:,:,idx);

        % strength 
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

wsbm_highDappear_subjAll = zeros([ templateModel.R_Struct.k nSubj ]);
mod_highDappear_subjAll = zeros([ templateModel.R_Struct.k nSubj ]);
% yeo_highDappear_subjAll = zeros([ 7 nSubj ]);

wsbm_highSappear_subjAll = zeros([ templateModel.R_Struct.k nSubj ]);
mod_highSappear_subjAll = zeros([ templateModel.R_Struct.k nSubj ]);
% yeo_highSappear_subjAll = zeros([ 7 nSubj ]);

highD_level = 75 ; 

% need to algn the subj to the ref
tmp_mod = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.mod) ;
% tmp_yeo = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.yeo) ;

% need to do this because after alingment, yeo will have gap in com vec
% yeo_lab = unique(tmp_yeo) ;

for idx = 1:nSubj

    currSubjData = subjDataMat(:,:,idx) ; 
    
    % get certain percentile of the data
    highD_thr = prctile(sum(currSubjData > 0,2),highD_level) ;
    highD_node = sum(currSubjData > 0,2) > highD_thr ;
    
    wsbm_highDappear_subjAll(:,idx) = histc(comVecs.wsbm(highD_node),1:templateModel.R_Struct.k) ;
    mod_highDappear_subjAll(:,idx) = histc(tmp_mod(highD_node),1:templateModel.R_Struct.k) ;
%     yeo_highDappear_subjAll(:,idx) = histc(tmp_yeo(highD_node),yeo_lab) ;
   
    highD_thr = prctile(sum(currSubjData,2),highD_level) ;
    highD_node = sum(currSubjData,2) > highD_thr ;
    
    wsbm_highSappear_subjAll(:,idx) = histc(comVecs.wsbm(highD_node),1:templateModel.R_Struct.k) ;
    mod_highSappear_subjAll(:,idx) = histc(tmp_mod(highD_node),1:templateModel.R_Struct.k) ;
%     yeo_highSappear_subjAll(:,idx) = histc(tmp_yeo(highD_node),yeo_lab) ;
    
end
  
% ICC
wsbm_highD_icc = IPN_icc(wsbm_highDappear_subjAll,3,'single') ;
mod_highD_icc = IPN_icc(mod_highDappear_subjAll,3,'single') ;
% yeo_highD_icc = IPN_icc(yeo_highDappear_subjAll,3,'single') ;

wsbm_highS_icc = IPN_icc(wsbm_highSappear_subjAll,3,'single') ;
mod_highS_icc = IPN_icc(mod_highSappear_subjAll,3,'single') ;
% yeo_highS_icc = IPN_icc(yeo_highSappear_subjAll,3,'single') ;

tmpICC = cell([3 1]);
tmpICC_2 = cell([3 1]);

tmpData = cell([3 1]);
tmpData{1} = wsbm_highDappear_subjAll ;
tmpData{2} = mod_highDappear_subjAll ;
% tmpData{3} = yeo_highDappear_subjAll ;

tmpData_2 = cell([3 1]);
tmpData_2{1} = wsbm_highSappear_subjAll ;
tmpData_2{2} = mod_highSappear_subjAll ;
% tmpData_2{3} = yeo_highSappear_subjAll ;

nBoot = 500 ;

btspCommMed = cell([3 1]) ;
btspCommMed{1} = zeros([templateModel.R_Struct.k nBoot]) ;
btspCommMed{2} = zeros([templateModel.R_Struct.k nBoot]) ;
% btspCommMed{3} =zeros([7 nBoot]) ;
tmpMed = cell([3 1]) ;

for rep = 1:2

    % lets see if we can bootstrap these values
    tmpRes = zeros([ nBoot 1 ]);
    tmpRes_2 = zeros([ nBoot 1 ]);
    
    % get bootstrp indicies
    [~,bootInd] = bootstrp(nBoot,@(a)[],1:nSubj) ;

    for idx = 1:nBoot

        tmpRes(idx) = IPN_icc(tmpData{rep}(:,bootInd(:,idx)),3,'single');
        tmpRes_2(idx) = IPN_icc(tmpData_2{rep}(:,bootInd(:,idx)),3,'single');
%         % lets also bootstrap the median of each community 
%         btspCommMed{rep}(:,idx) = median(tmpData{rep}(:,bootInd(:,idx)),2) ;
        
    end

    tmpICC{rep} = prctile(tmpRes,[2.5 97.5]) ;
    tmpICC_2{rep} = prctile(tmpRes_2,[2.5 97.5]) ;
%     tmpMed{rep} = prctile(btspCommMed{rep},[2.5 97.5],2) ;
    
end

wsbm_highD_icc_c95 = tmpICC{1};
mod_highD_icc_c95 = tmpICC{2};
% yeo_highD_icc_c95 = tmpICC{3};

% wsbm_highD_med = [ (tmpMed{1}(:,1) - 0.1)  (tmpMed{1}(:,2) + 0.1) ] ;
% mod_highD_med = [ (tmpMed{2}(:,1) - 0.1)  (tmpMed{2}(:,2) + 0.1) ] ;
% yeo_highD_med = [ (tmpMed{3}(:,1) - 0.1)  (tmpMed{3}(:,2) + 0.1) ] ;

% %% plot it
% default_cmap = [0    0.4470    0.7410 ;
%                 0.8500    0.3250    0.0980;
%                 0.9290    0.6940    0.1250 ] ;
% 
% ddd = wsbm_highDappear_subjAll ;
% ddd(ddd==0) = NaN ;            
% 
% % make different percentiles of the data          
% fullrange = prctile(ddd,[0 100],2);  
% p50 = prctile(ddd,[25 75],2);  
% p25 = prctile(ddd,[37.5 62.5],2);  
% med = prctile(ddd,[50 50],2);  
% 
% rb = rangebar(fullrange,0.8);
% rb.FaceColor = default_cmap(1,:);
% rb.EdgeAlpha = 0 ;
% rb.FaceAlpha = 0.1 ;
% 
% ax = gca ;
% 
% rb2 = rangebar(ax,p50,0.7);
% rb2.FaceColor = default_cmap(1,:);
% rb2.EdgeAlpha = 0 ;
% rb2.FaceAlpha = 0.2 ;
% 
% rb = rangebar(ax,p25,0.6);
% rb.FaceColor = default_cmap(1,:);
% rb.EdgeAlpha = 0 ;
% rb.FaceAlpha = 0.3 ;
% 
% rb = rangebar(ax,med,0.9);
% rb.FaceColor = default_cmap(1,:);
% rb.EdgeAlpha = 0 ;
% rb.FaceAlpha = 1 ;
% 
% % wsbm_highD_med = tmpMed{1} ;
% % mod_highD_med = tmpMed{2} ;
% % yeo_highD_med = tmpMed{3} ;

%% make a table

% modelNames = { 'wsbm' 'mod' 'yeo' } ;
modelNames = { 'wsbm' 'mod' } ;

% com_names = { '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' }' ;
com_names = { '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11'}' ;

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

% multipleTablesTemplate = cell([3 1]);
% multipleTablesSubj = cell([3 1]);
multipleTablesTemplate = cell([2 1]);
multipleTablesSubj = cell([2 1]);

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

%% make a table that combines the stats

% % wsbm
% multipleTablesTemplate{mN} = table() ;
% multipleTablesTemplate{mN}.mean_within_weight = round(tmpl_strength.(modelNames{mN}).within,roundTo,'significant') ;
% multipleTablesTemplate{mN}.mean_btwn_weight = round(tmpl_strength.(modelNames{mN}).between,roundTo,'significant') ;
% multipleTablesTemplate{mN}.q_prcnt = round(tmpl_mod.(modelNames{mN})./ ...
%     sum(tmpl_mod.(modelNames{mN})),roundTo,'significant') ;
% multipleTablesTemplate{mN}.mean_comm_parti_coef = round(tmpl_parti.comm.(modelNames{mN}),roundTo,'significant');
% multipleTablesTemplate{mN}.assort_of_comm = round(tmpl_assort.(modelNames{mN}),roundTo,'significant') ;

tempNSubj_tbl = cell([ 3 1 ]) ;

for mN = 1:length(modelNames)

    tempNSubj_tbl{mN} = table() ;
        
    for idx = 1:5
    
        tempNSubj_tbl{mN}.(colNames{idx}) = multipleTablesTemplate{mN}.(colNames{idx}) ;
        tempNSubj_tbl{mN}.(strcat(colNames{idx},'_meanSubj')) =  multipleTablesSubj{mN}.(colNames2{((idx-1)*3) + 1}) ;
        tempNSubj_tbl{mN}.(strcat(colNames{idx},'_stdSubj')) =  multipleTablesSubj{mN}.(colNames2{((idx-1)*3) + 3}) ;
    end
    
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

%% write the combo tables to csv

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_wsbm_combo_stats.csv');
writetable(tempNSubj_tbl{1},fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_mod_combo_stats.csv');
writetable(tempNSubj_tbl{2},fileName,'WriteRowNames',true)

fileName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_yeo_combo_stats.csv');
writetable(tempNSubj_tbl{3},fileName,'WriteRowNames',true)

