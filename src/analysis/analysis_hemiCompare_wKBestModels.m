%% HEMISPHERE ANALYSIS 
% 
% clc
% clearvars

%% load the necessary data

% LOAD WORKSPACE: findK_1p1.mat

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
load(loadName) ;
 
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
load(loadName) ;

%% lets look at how some metrics look across lifespan, using static
% definition of paritions

nSubj = length(dataStruct) ;
nBlocks = templateModel.R_Struct.k ;
nNodes = templateModel.Data.n ;
nModels = length(kLoopBestModels);

nodeAssort_mat_wsbm = zeros([nModels nNodes nSubj]);
nodeAssort_mat_mod = zeros([nNodes nSubj]);

parti_mat_wsbm = zeros([nModels nNodes nSubj]) ;
parti_mat_mod = zeros([nNodes nSubj]);

zscore_mat_wsbm = zeros([nModels nNodes nSubj]) ;
zscore_mat_mod = zeros([nNodes nSubj]) ;

subjDataMat = zeros([ nNodes nNodes nSubj ]);

%% gather subject data

for idx=1:nSubj
    
    % extract a simple graph, replace NaN with 0
    tmpAdj = dataStruct(idx).countVolNormMat(selectNodesFrmRaw, selectNodesFrmRaw);
    % get rid of the diagonal
    %n=size(tmpAdj,1);
    tmpAdj(1:nNodes+1:end) = 0; 
    % mask out AdjMat entries below mask_thr
    tmpAdj_mask = dataStruct(idx).countMat(selectNodesFrmRaw, selectNodesFrmRaw) > 1 ;    
    tmpAdj_mask(tmpAdj_mask > 0) = 1 ;   
    tmpAdj = tmpAdj .* tmpAdj_mask ;
    tmpAdj(isnan(tmpAdj)) = 0 ;
    
    subjDataMat(:,:,idx) = tmpAdj ;
    
end

%% evaluate symmetry

% make this variable just for sanity check that we extracting different ca
kBestCa = zeros([nNodes nModels]) ;

% across 100 wsbm models
for mM = 1:nModels
    
    [~,tmpCA] = community_assign(kLoopBestModels(mM).Model) ;
    
    kBestCa(:,mM) = tmpCA ;
    
% across subject at each of 100 models
% get the community-based node-wise statistic
for idx=1:nSubj

    % 1 of 620 subject data
    tmpAdj = subjDataMat(:,:,idx) ; 
    
    [~,nodeAssort_mat_wsbm(mM,:,idx)] = ...
        eval_com_assortatvity_wu(tmpAdj,tmpCA);
    
    parti_mat_wsbm(mM,:,idx) = participation_coef(tmpAdj,tmpCA);

    zscore_mat_wsbm(mM,:,idx) = module_degree_zscore(tmpAdj,tmpCA,0);
  
end

    disp(mM)

end

% and collect the modular stats 
for idx=1:nSubj

    tmpAdj = subjDataMat(:,:,idx) ; 
    
    [~,nodeAssort_mat_mod(:,idx)] = ...
        eval_com_assortatvity_wu(tmpAdj,comVecs.mod);
    
    parti_mat_mod(:,idx) = participation_coef(tmpAdj,comVecs.mod);

    zscore_mat_mod(:,idx) = module_degree_zscore(tmpAdj,comVecs.mod,0);
  
end

%% hemicompare across all 100 kLoopBestModels

% l_hemi = 1:57 ;
% r_hemi = 58:114 ;
l_hemi = ~~([ ones(57,1) ; zeros(57,1) ]) ;
r_hemi = ~~([ zeros(57,1) ; ones(57,1) ]) ;

wsbm_ks_parti = zeros([nModels nSubj]);
wsbm_ks_zscore = zeros([nModels nSubj]);
wsbm_ks_assort = zeros([nModels nSubj]);

mod_ks_parti = zeros([ nSubj 1]);
mod_ks_zscore = zeros([ nSubj 1]);
mod_ks_assort = zeros([ nSubj 1]);

% across the 100 wsbm models 
for mM=1:nModels

for idx=1:nSubj
    
    [~,wsbm_ks_parti(mM,idx)] = simple_hist_dist(parti_mat_wsbm(mM,l_hemi,idx)',parti_mat_wsbm(mM,r_hemi,idx)');
    [~,wsbm_ks_zscore(mM,idx)] = simple_hist_dist(zscore_mat_wsbm(mM,l_hemi,idx)',zscore_mat_wsbm(mM,r_hemi,idx)');
    [~,wsbm_ks_assort(mM,idx)] = simple_hist_dist(nodeAssort_mat_wsbm(mM,l_hemi,idx)',nodeAssort_mat_wsbm(mM,r_hemi,idx)');
    
end

end

% and the mod
for idx=1:nSubj
    
    [~,mod_ks_parti(idx)] = simple_hist_dist(parti_mat_mod(l_hemi,idx),parti_mat_mod(r_hemi,idx));
    [~,mod_ks_zscore(idx)] = simple_hist_dist(zscore_mat_mod(l_hemi,idx),zscore_mat_mod(r_hemi,idx));
    [~,mod_ks_assort(idx)] = simple_hist_dist(nodeAssort_mat_mod(l_hemi,idx),nodeAssort_mat_mod(r_hemi,idx));
    
end

%% compare it to the histogram we originally got....
% to show that sym. is better across kLooper best models

for idx = 1:nModels

    diff(idx) = median(wsbm_ks_parti(idx,:)) - median(mod_ks_parti) ;
    [rs_p(idx),~,rsStat(idx)] = ranksum(wsbm_ks_parti(idx,:)',mod_ks_parti) ;
    zvals(idx) = rsStat(idx).zval ;
end

% zvals are significantly lower
sum(zvals(rs_p < 10e-9) & (zvals(rs_p < 10e-9) < 0))

% range
min(zvals((rs_p < 10e-9) & (zvals < 0)))
max(zvals((rs_p < 10e-9) & (zvals < 0)))

