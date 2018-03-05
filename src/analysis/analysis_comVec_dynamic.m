%% ANALYSIS LOOKING AT MODEL PREDICTION

clc
clearvars

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

outIntermPrefix = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR);
outProcessPrefix = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR);

load(strcat(outProcessPrefix,'_fit_wsbm_script_v7p3.mat'))
load(strcat(outProcessPrefix,'_fit_mod_v7p3.mat'))
load(strcat(outIntermPrefix,'_comVecs.mat'))
%load('data/interim/subj_dataStruct.mat')
load(strcat(outIntermPrefix,'_templateModel_1.mat'))
load(strcat(outProcessPrefix,'_basicData_v7p3.mat'))

%% lets look at how some metrics look across lifespan, using changing ca
% definition of paritions

nSubj = length(dataStruct) ;
nBlocks = templateModel.R_Struct.k ;
nNodes = templateModel.Data.n ;

totDensity = zeros([nSubj 1]);

nBlockInteract = (nBlocks^2 + nBlocks) / 2 ;

getIdx = ~~triu(ones(nBlocks));

subjDataMat = zeros([ nNodes nNodes nSubj ]);
subjWsbmCA = zeros([ nNodes nSubj ]);
subjModCA = zeros([ nNodes nSubj ]);

%% unroll block vec vars

blockVec_wsbm = zeros([nBlockInteract nSubj]);
blockVec_mod = zeros([nBlockInteract nSubj]) ;

avgBlockVec_wsbm = zeros([nBlockInteract nSubj]);
avgBlockVec_mod = zeros([nBlockInteract nSubj]) ;

binBlockVec_wsbm = zeros([nBlockInteract nSubj]);
binBlockVec_mod = zeros([nBlockInteract nSubj]) ;

%% also community comvec
com_avgBlockVec_wsbm = zeros([nBlocks nBlocks nSubj]);
com_avgBlockVec_mod = zeros([nBlocks nBlocks nSubj]) ;

%% gather some of the data
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
  
    totDensity(idx) = sum(tmpAdj(:));
    
    % get the fit community assignments 
    [~,tmpCA] = community_assign(fitWSBMAllStruct(idx).centModel) ;
    subjWsbmCA(:,idx) = CBIG_HungarianClusterMatch(comVecs.wsbm,tmpCA);
    subjModCA(:,idx) = CBIG_HungarianClusterMatch(comVecs.mod,...
        fitModAllStruct(idx).caMod_subj);
    
end

%% some runs of modularity did not render paritions with correct num coms
% and therefore have community assignment of all 0's
% we need to record these to skip over them

modExclude = (sum(subjModCA) == 0) ;

%% iterate
for idx = 1:nSubj
       
    % wsbm
    [tmpBl,avgtmpbl,~,bintmpbl] = get_block_mat(tmpAdj,subjWsbmCA(:,idx));
    blockVec_wsbm(:,idx) = tmpBl(getIdx);
    avgBlockVec_wsbm(:,idx) = avgtmpbl(getIdx);
    binBlockVec_wsbm(:,idx) = bintmpbl(getIdx);
   
    com_avgBlockVec_wsbm(:,:,idx) = avgtmpbl ;
    
    if modExclude(idx) == 0
    
        % mod
        [tmpBl,avgtmpbl,~,bintmpbl] = get_block_mat(tmpAdj,subjModCA(:,idx));
        blockVec_mod(:,idx) = tmpBl(getIdx);
        avgBlockVec_mod(:,idx) = avgtmpbl(getIdx);
        binBlockVec_mod(:,idx) = bintmpbl(getIdx);
        com_avgBlockVec_mod(:,:,idx) = avgtmpbl;  
        
    else
        
        % set zeros
        blockVec_mod(:,idx) = zeros([nBlockInteract 1]);
        avgBlockVec_mod(:,idx) = zeros([nBlockInteract 1]);
        binBlockVec_mod(:,idx) = zeros([nBlockInteract 1]);
        com_avgBlockVec_mod(:,:,idx) = zeros(nBlocks);  
          
    end
    
end
  
%% correlate each block pattern with model

% now different is predicted from actual?
templateData = templateModel.Data.Raw_Data;
templateData(isnan(templateData)) = 0 ;

% model_edgeVec = templateModel.Para.predict_e ;
% model_weiVec = templateModel.Para.predict_w ;
wsbmPred_weiVec =  templateModel.Para.predict_e .*  templateModel.Para.predict_w; 
wsbmPred_weiVec(isnan(wsbmPred_weiVec)) = 0;
[~,wsbm_weiVec_sqr] = make_square(wsbmPred_weiVec);

% and now creat the 'model' for the modular parition
[~,mod_weiVec_empir] = get_block_mat(templateData,comVecs.mod);
mod_weiVec_empir = mod_weiVec_empir(getIdx);
[~,mod_weiVec_empir_sqr] = make_square(mod_weiVec_empir);

%% community pattern vars

wsbm_comms_weiVec_corr = zeros([nBlocks nSubj]) ;
mod_comms_weiVec_corr = zeros([nBlocks nSubj]) ;

wsbm_comms_weiVec_cb = zeros([nBlocks nSubj]) ;
mod_comms_weiVec_cb = zeros([nBlocks nSubj]) ;

wsbm_comms_weiVec_cos = zeros([nBlocks nSubj]) ;
mod_comms_weiVec_cos = zeros([nBlocks nSubj]) ;

%% iterate
for idx = 1:nSubj
   
    % correlation 
    corrStr='pearson';
    wsbm_weiVec_corr(idx) = corr(wsbmPred_weiVec,avgBlockVec_wsbm(:,idx),...
        'type',corrStr) ;
    mod_weiVec_corr(idx) = corr(mod_weiVec_empir,avgBlockVec_mod(:,idx),...
        'type',corrStr) ;
    
    % other distances   
    wsbm_weiVec_cb(idx) = pdist([wsbmPred_weiVec' ; avgBlockVec_wsbm(:,idx)'],'cityblock') ;
    mod_weiVec_cb(idx) = pdist([mod_weiVec_empir' ; avgBlockVec_mod(:,idx)'],'cityblock') ;

    wsbm_weiVec_eud(idx) = pdist([wsbmPred_weiVec' ; avgBlockVec_wsbm(:,idx)'],'euclidean') ;
    mod_weiVec_eud(idx) = pdist([mod_weiVec_empir' ; avgBlockVec_mod(:,idx)'],'euclidean') ;
    
    wsbm_weiVec_cos(idx) = 1 - pdist([wsbmPred_weiVec' ; avgBlockVec_wsbm(:,idx)'],'cosine') ;
    mod_weiVec_cos(idx) = 1 - pdist([mod_weiVec_empir' ; avgBlockVec_mod(:,idx)'],'cosine') ;
    
    % distances of each community
    for jdx=1:nBlocks
        
        % correlation 
        wsbm_comms_weiVec_corr(jdx,idx) = corr(wsbm_weiVec_sqr(:,jdx),...
            com_avgBlockVec_wsbm(:,jdx,idx),...
            'type',corrStr) ;
        mod_comms_weiVec_corr(jdx,idx) = corr(mod_weiVec_empir_sqr(:,jdx),...
            com_avgBlockVec_mod(:,jdx,idx),...
            'type',corrStr) ;
        
        % other distances
        wsbm_comms_weiVec_cb(jdx,idx) = pdist([wsbm_weiVec_sqr(:,jdx)' ;...
            com_avgBlockVec_wsbm(:,jdx,idx)' ],...
            'cityblock') ;
        mod_comms_weiVec_cb(jdx,idx) = pdist([mod_weiVec_empir_sqr(:,jdx)' ; ...
            com_avgBlockVec_mod(:,jdx,idx)'],...
            'cityblock') ;
 
        wsbm_comms_weiVec_cos(jdx,idx) = 1 - pdist([wsbm_weiVec_sqr(:,jdx)' ;...
            com_avgBlockVec_wsbm(:,jdx,idx)' ],...
            'cosine') ;
        mod_comms_weiVec_cos(jdx,idx) = 1 - pdist([mod_weiVec_empir_sqr(:,jdx)' ; ...
            com_avgBlockVec_mod(:,jdx,idx)'],...
            'cosine') ; 
    end
    
end

%% look and how community distances are distributed

wsbm_com_totLen = zeros([ nBlocks nSubj ]);
mod_com_totLen = zeros([ nBlocks nSubj ]);

wsbm_com_totDist = zeros([ nBlocks nSubj ]);
mod_com_totDist = zeros([ nBlocks nSubj ]);

subjLensMat = zeros([ nNodes nNodes nSubj ]);
subjDistMat = zeros([ nNodes nNodes nSubj ]);

% gather some data
for idx = 1:nSubj

    % mask out AdjMat entries below mask_thr
    tmpMask = dataStruct(idx).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > MASK_THR_INIT ;    
    tmpMask(tmpMask > 0) = 1 ;   
    
    % get streamline lengths
    tmpLens = dataStruct(idx).lensMat(selectNodesFrmRaw,selectNodesFrmRaw);
    %tmpAdj(isnan(tmpAdj)) = 0 ;
    subjLensMat(:,:,idx) = tmpLens .* tmpMask ;
    
    % get eud distances
    subjDistMat(:,:,idx) = dataStruct(idx).distCoorMM(selectNodesFrmRaw,selectNodesFrmRaw) ;

end
    
for idx = 1:nSubj
    
    % get the length matrix
    tmpAdj = subjLensMat(:,:,idx) ;
    
    wsbm_tmpBlMat = get_block_mat(tmpAdj,subjWsbmCA(:,idx)) ;
    wsbm_com_totLen(:,idx) = sum(wsbm_tmpBlMat(~~eye(nBlocks)),2) ;
    
    if modExclude(idx) == 0
        mod_tmpBlMat = get_block_mat(tmpAdj,subjModCA(:,idx));
        mod_com_totLen(:,idx) = sum(mod_tmpBlMat(~~eye(nBlocks)),2) ; 
    else
        mod_com_totLen(:,idx) = zeros([nBlocks 1]) ; 
    end
    
    % get the dist matrix
    tmpAdj = subjDistMat(:,:,idx) ;
    
    wsbm_tmpBlMat = get_block_mat(tmpAdj,subjWsbmCA(:,idx)) ;
    wsbm_com_totDist(:,idx) = sum(wsbm_tmpBlMat(~~eye(nBlocks)),2) ;
    
    if modExclude(idx) == 0
        mod_tmpBlMat = get_block_mat(tmpAdj,subjModCA(:,idx));
        mod_com_totDist(:,idx) = sum(mod_tmpBlMat(~~eye(nBlocks)),2) ; 
    else
        mod_com_totDist(:,idx) = zeros([nBlocks 1]) ; 
    end
    
end

%% quick analysis... lets look at variability of nodal assignment

wsbm_vers = get_nodal_versatility(subjWsbmCA(:,~modExclude)) ;
mod_vers = get_nodal_versatility(subjModCA(:,~modExclude));

% and gets look at consensus in age bins...
% number of actual groups is thresholds + 1
thresholds = 4 ;
age_bins = thresholds + 1;

[agesSorted,ageSortIdx] = sort(datasetDemo.age) ;

low_quantile = [ 0 quantile(1:length(datasetDemo.age),4) ] ; 
high_quantile = [ quantile(1:length(datasetDemo.age),4) (length(datasetDemo.age)+1) ] ;

templateIdMat = zeros( [ size(datasetDemo.age,1) (thresholds+1) ] ) ;

for templateIdx=1:age_bins

    subjectsIdxVec = ageSortIdx > low_quantile(templateIdx) & ...
        ageSortIdx < high_quantile(templateIdx) ;
    
    templateIdMat(:,templateIdx) = ~~subjectsIdxVec ;
end

templateIdMat2 = ~~templateIdMat(~modExclude,:) ;
subjWsbmCA2 = subjWsbmCA(:,~modExclude) ;
subjModCA2 = subjModCA(:,~modExclude);

wsbm_agebin_vers = zeros([ nNodes 5]) ;
mod_agebin_vers = zeros([ nNodes 5]) ;

for idx=1:age_bins
    
    wsbm_agebin_vers(:,idx) = get_nodal_versatility(subjWsbmCA2(:,templateIdMat2(:,idx))) ;
    mod_agebin_vers(:,idx) = get_nodal_versatility(subjModCA2(:,templateIdMat2(:,idx))) ; 
end



%% some sort of comparison between the cos and cb overall

[~,tb1,anova2_stat_cb] = anova2([wsbm_weiVec_cb(~modExclude)' mod_weiVec_cb(~modExclude)' ],1,'off') ;
[~,tb2,anova2_stat_cos] = anova2([wsbm_weiVec_cos(~modExclude)' mod_weiVec_cos(~modExclude)' ],1,'off') ;

% % same as above
% [~,tb12] = anova_rm([wsbm_weiVec_cb' mod_weiVec_cb' yeo_weiVec_cb'],'off') ;
% [~,tb22] = anova_rm([wsbm_weiVec_cos' mod_weiVec_cos' yeo_weiVec_cos'],'off') ;

anova1_stat_cb_mult = multcompare(anova2_stat_cb) ;
anova1_stat_cos_mult = multcompare(anova2_stat_cos) ;

% or just a paired ttest too...
[tt_h,tt_p,tt_ci,tt_stats] = ttest2(wsbm_weiVec_cb(~modExclude)',...
    mod_weiVec_cb(~modExclude)','vartype','unequal') ;

[tt_h,tt_p,tt_ci,tt_stats] = ttest2(wsbm_weiVec_cos(~modExclude)',...
    mod_weiVec_cos(~modExclude)','vartype','unequal') ;

%% and ttest on node versatility

[tt_h,tt_p,tt_ci,tt_stats] = ttest2(wsbm_vers',...
    mod_vers','vartype','unequal') ;

%% bootstrapping & perms

emp_diff = wsbm_vers' - mod_vers' ;

% get book indicies 
nBoot = 500 ;
vers_btsp_res = zeros([ nNodes nBoot ]);
[~,bootInd] = bootstrp(nBoot,@(a)[],1:sum(~modExclude)) ;

wsbmCA = subjWsbmCA(:,~modExclude) ;
modCA = subjModCA(:,~modExclude) ;

for idx = 1:nBoot

    disp(idx)
    
    % get bootstrapped sample of wsbm and mod community strucutre 
    wsbm_btsp_smpl = wsbmCA(:,bootInd(:,idx)) ;
    mod_btsp_smpl = modCA(:,bootInd(:,idx)) ;
    
    wsbm_btsp_vers = get_nodal_versatility(wsbm_btsp_smpl) ;
    mod_btsp_vers = get_nodal_versatility(mod_btsp_smpl) ;
    
    vers_btsp_res(:,idx) = wsbm_btsp_vers - mod_btsp_vers ;
        
end

%% perms

nPerms = 10000 ;
permInd = randi([ 0 1 ] , [nPerms , sum(~modExclude)]) ; 

comboCA1 = zeros([nNodes sum(~modExclude)]) ;
comboCA2 = zeros([nNodes sum(~modExclude)]) ;

permsRes = zeros([ nNodes nPerms ]);

for idx = 1:nPerms 

    disp(idx) 
    
    comboCA1(:,~~permInd(idx,:)) = wsbmCA(:,~~permInd(idx,:)) ;
    comboCA1(:,~permInd(idx,:)) = modCA(:,~permInd(idx,:)) ;

    comboCA2(:,~permInd(idx,:)) = wsbmCA(:,~permInd(idx,:)) ;
    comboCA2(:,~~permInd(idx,:)) = modCA(:,~~permInd(idx,:)) ;
    
    tmp1 = get_nodal_versatility(comboCA1) ;
    tmp2 = get_nodal_versatility(comboCA2) ;
    
    permsRes(:,idx) = tmp1 - tmp2 ;
end

pval_emp_diff = zeros([nNodes 1]) ;
% permStruct.permPvalR2 = (sum(permDist > lsFitStruct.R2) + 1) / (permIter + 1) ;

for idx = 1:nNodes
    pval_emp_diff(idx) = (sum( abs(permsRes(idx,:)) > abs(emp_diff(idx)) ) + 1) / ( nPerms + 1) ; 
end

pval_bonf = 0.05 / nNodes ; 
passBonf = pval_emp_diff < pval_bonf ;

passBonfthr_emp_diff = emp_diff .* passBonf ;
passBonfthr_emp_diff(passBonfthr_emp_diff == 0) = NaN ;

%% gather info

vers_btsp_CI = prctile(vers_btsp_res',[2.5 97.5])' ;

% check which do not cross 0
btsp_good_ci = ((vers_btsp_CI(:,1) < 0) & (vers_btsp_CI(:,2) < 0)) | ...
    ((vers_btsp_CI(:,1) > 0) & (vers_btsp_CI(:,2) > 0));

btsp_mean_diff = btsp_good_ci .* mean(vers_btsp_res,2) ;
btsp_mean_diff(btsp_mean_diff == 0) = NaN ;

emp_diff_pass = btsp_good_ci .* emp_diff ;
emp_diff_pass(emp_diff_pass == 0) = NaN ;

%% find vals in temporal

lab = load('data/raw/NKIen1/yeo/nodeLabels.mat') ;  
temp_idx = regexp(lab.nodeLabels,'.*Temp.*','once') ;
temp_idx = ~cellfun(@isempty,temp_idx) ;

mean(emp_diff_pass(temp_idx))

[~,ind] = max(emp_diff_pass) ;

%% save it

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_comVec_dynamic_results.mat');
save(outName,...
    '*_weiVec_corr',...
    '*_weiVec_cb',...
    '*_weiVec_eud',...
    '*_weiVec_cos',...
    '*_vers',...
    'totDensity',...
    'modExclude',...
    'wsbm_vers','mod_vers',...
    'vers_btsp_CI','btsp_mean_diff',...
    'emp_diff','pval_emp_diff','btsp_good_ci','btsp_mean_diff',...
    'passBonfthr_emp_diff',...
    ...
    '-v7.3')


%% get agreement matrices

wsbm_agree = agreement(subjWsbmCA2) ./ size(subjWsbmCA2,2) ;
mod_agree = agreement(subjModCA2) ./ size(subjModCA2,2) ;

[~,wsbm_avgProb ] = get_block_mat(wsbm_agree, comVecs.wsbm) ;
[~,mod_avgProb ] = get_block_mat(mod_agree, comVecs.mod) ;


%% test viz

% inputAgreeData = cell([2 1]) ;
% inputAgreeData{1} = wsbm_agree ;
% inputAgreeData{2} = mod_agree ;
% 
% nNodes = templateModel.Data.n ;
% nComm = templateModel.R_Struct.k ; 
% 
% parcels = {'wsbm' 'mod' };
% parcelName = {'WSBM' 'Modular'} ;
% 
% % % get the template data
% % templateData = templateModel.Data.Raw_Data ;
% % templateData(isnan(templateData)) = 0;
% % 
% % % scramble
% % scrmb = randperm(nNodes);
% % templateData = templateData(scrmb,scrmb) ;
% % %comLabels = comVecs.wsbm(scrmb) ;
% 
% nComm4ColorMap = templateModel.R_Struct.k ;
% 
% for fig = 1:length(parcels)
% %for fig = 3
%  
%     figure
%     
%     comLabels = comVecs.((lower(parcels{fig}))) ;    
%     comLabels = CBIG_HungarianClusterMatch(comVecs.wsbm,comLabels);
%     comLabels = comLabels(scrmb);
%     
%     nComm = length(unique(comLabels)) ;
%     uniqueLab = unique(comLabels) ;
%     
%     [xOnDiag,yOnDiag,sortIdx] = grid_communities(comLabels);
% 
%     % code from WSBM stuffs
%     A_sort = zeros(nNodes);
%     list = zeros(1,nNodes);
%     breaks = zeros(1,nComm);
%     cur = 1;
%     for idx = 1:nComm
%         indicies = find(comLabels == uniqueLab(idx));
%         
%         if isempty(indicies)
%          
%             list(cur) = [] ; %list(cur:cur+length(indicies)-1) = indicies;
%             cur = 0 ;%cur + length(indicies);
%             breaks(idx) = cur ;%cur-1;
%             
%         end
%         
%         list(cur:cur+length(indicies)-1) = indicies;
%         cur = cur + length(indicies);
%         breaks(idx) = cur-1;
%     end
%     for idx = 1:nNodes
%         A_sort(idx,:) = inputAgreeData{fig}(list(idx),list);
%     end
% 
%     %Plot the Matrix
%     h = imagesc(A_sort,[min(A_sort(:))-.00001,max(A_sort(:))]);
%     colormap(brewermap(nComm,'BuPu'))
%     %set(h,'alphadata',(A_sort > 0) .* 0.15);
% 
%     % fix image properties 
%     ax = gca ;
%     axis square
%     % axis([0.5 (nNodes+0.5) 0.5 (nNodes+0.5)]); 
%     % ax.Box = 'on' ;
%     % set(ax,'Ydir','reverse');
%     % set(ax,'ytick',[])
%     % set(ax,'xtick',[])
% 
%     hold on
% 
%     % plot off diagonal
%     for idx = 1:(nComm-1)
% 
%         lineWidth = 1.5;
%         offDiagColor = [1 0 0 0.25] ; 
% 
%         % vertical   
%         plot([breaks(idx)+0.5,breaks(idx)+0.5],[breaks(idx+1)+0.5,breaks(nComm)+.5],...
%             'Color',offDiagColor,'LineWidth',lineWidth);
%         if idx > 1
%             plot([breaks(idx)+0.5,breaks(idx)+0.5],[-0.5,breaks(idx-1)+0.5],...
%                 'Color',offDiagColor,'LineWidth',lineWidth);
%         end
% 
%         % horizontal
%         plot([breaks(idx+1)+0.5,breaks(nComm)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
%             'Color',offDiagColor,'LineWidth',lineWidth);
%         if idx > 1
%             plot([-.5,breaks(idx-1)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
%                 'Color',offDiagColor,'LineWidth',lineWidth);
%         end
% 
%     end      
% 
%     % MAKE THIS TEN...
%     cmap_mod = brewermap(nComm4ColorMap,'paired') ;
% 
%     for idx=0:(nComm-1)
% 
%         plot(xOnDiag( (idx*6)+1:((idx+1)*6) ) ,...
%             yOnDiag(  (idx*6)+1:((idx+1)*6) ) ,...
%             'Color',cmap_mod(uniqueLab(idx+1),:),'linewidth',3.5)
%     end
% 
%     % compute some yticks
%     breaks2 = [ 0 breaks ] ;
%     midlabelpoint = zeros([nComm 1]);
%     for idx = 1:length(breaks)
%         midlabelpoint(idx) = floor( (breaks2(idx+1) - breaks2(idx)) / 2) + breaks2(idx);  
%     end
% 
%     set(ax,'xtick',midlabelpoint)
%     set(ax,'xticklabel',uniqueLab)
%     set(ax,'ticklength',[ 0 0]) 
% 
%     set(ax,'ytick',midlabelpoint)
%     set(ax,'yticklabel',uniqueLab)
%     set(ax,'ticklength',[ 0 0]) 
% 
%     set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.4, 0.8]);
% 
%     ax.Title.String = { parcelName{fig}, ' community colors'};
%     ax.TitleFontSizeMultiplier = 1.5 ;
% 
% %     % save it
% %     fileName = strcat(parcels{fig},'_commView.png');
% %     ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
% %     %set(gcf,'paperpositionmode','auto');
% %     print(gcf,'-dpng','-r500',ff);
% %     close(gcf)
%     
% end

















