%% HEMISPHERE ANALYSIS 

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

%% lets look at how some metrics look across lifespan, using static
% definition of paritions

nSubj = length(dataStruct) ;
nBlocks = templateModel.R_Struct.k ;
nNodes = templateModel.Data.n ;

[~,ageSortIdx] = sort(datasetDemo.age);

assort_vec = zeros([nSubj 1]) ;
totDen_vec = zeros([nSubj 1]) ;

% q_mat = zeros([nBlocks nSubj]);
% q_mat_mod = zeros([nBlocks nSubj]);
% q_mat_yeo = zeros([7 nSubj]);

comAssort_mat = zeros([nBlocks nSubj]);
nodeAssort_mat = zeros([nNodes nSubj]);
comAssort_mat_mod = zeros([nBlocks nSubj]);
nodeAssort_mat_mod = zeros([nNodes nSubj]);
% comAssort_mat_yeo = zeros([7 nSubj]);
% nodeAssort_mat_yeo = zeros([nNodes nSubj]);

parti_mat_wsbm = zeros([nNodes nSubj]) ;
parti_mat_mod = zeros([nNodes nSubj]);
% parti_mat_yeo = zeros([nNodes nSubj]);

zscore_mat_wsbm = zeros([nNodes nSubj]) ;
zscore_mat_mod = zeros([nNodes nSubj]) ;
% zscore_mat_yeo = zeros([nNodes nSubj]) ;

nNodes = templateModel.Data.n ;
subjDataMat = zeros([ nNodes nNodes nSubj ]);

% loop through subjects 
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
    
    assort_vec(idx) = assortativity_wei(tmpAdj,0);
    totDen_vec(idx) = sum(sum(tmpAdj)) ; 
    
%     [~,q_mat(:,idx)] = eval_modularity_wu(tmpAdj,comVecs.wsbm) ;
%     [~,q_mat_mod(:,idx)] = eval_modularity_wu(tmpAdj,comVecs.mod) ;
%     [~,q_mat_yeo(:,idx)] = eval_modularity_wu(tmpAdj,comVecs.yeo) ;
    
    [comAssort_mat(:,idx),nodeAssort_mat(:,idx)] = ...
        eval_com_assortatvity_wu(tmpAdj,comVecs.wsbm);
    [comAssort_mat_mod(:,idx),nodeAssort_mat_mod(:,idx)] = ...
        eval_com_assortatvity_wu(tmpAdj,comVecs.mod);  
%     [comAssort_mat_yeo(:,idx),nodeAssort_mat_yeo(:,idx)] = ...
%         eval_com_assortatvity_wu(tmpAdj,comVecs.yeo); 
    
    parti_mat_wsbm(:,idx) = participation_coef(tmpAdj,comVecs.wsbm);
    parti_mat_mod(:,idx) = participation_coef(tmpAdj,comVecs.mod);
%     parti_mat_yeo(:,idx) = participation_coef(tmpAdj,comVecs.yeo);

    zscore_mat_wsbm(:,idx) = module_degree_zscore(tmpAdj,comVecs.wsbm,0);
    zscore_mat_mod(:,idx) = module_degree_zscore(tmpAdj,comVecs.mod,0);
%     zscore_mat_yeo(:,idx) = module_degree_zscore(tmpAdj,comVecs.yeo,0);

    
end

%% look at hemispheric differences

if strcmp(PARCELLATION,'yeo')
    l_hemi = 1:57 ;
    r_hemi = 58:114 ;
elseif strcmp(PARCELLATION,'scale125')
    l_hemi = 1:111 ;
    r_hemi = 112:218 ;
end

wsbm_ks_parti = zeros([nSubj 1]);
mod_ks_parti = zeros([nSubj 1]);

wsbm_ks_zscore = zeros([nSubj 1]);
mod_ks_zscore = zeros([nSubj 1]);

wsbm_ks_assort = zeros([nSubj 1]);
mod_ks_assort = zeros([nSubj 1]);

wsbm_ks_fibLen = zeros([nSubj 1]);
mod_ks_fibLen = zeros([nSubj 1]);

for idx=1:nSubj
    
    [~,wsbm_ks_parti(idx)] = simple_hist_dist(parti_mat_wsbm(l_hemi,idx),parti_mat_wsbm(r_hemi,idx));
    [~,mod_ks_parti(idx)] = simple_hist_dist(parti_mat_mod(l_hemi,idx),parti_mat_mod(r_hemi,idx));
    
    [~,wsbm_ks_zscore(idx)] = simple_hist_dist(zscore_mat_wsbm(l_hemi,idx),zscore_mat_wsbm(r_hemi,idx));
    [~,mod_ks_zscore(idx)] = simple_hist_dist(zscore_mat_mod(l_hemi,idx),zscore_mat_mod(r_hemi,idx));

    [~,wsbm_ks_assort(idx)] = simple_hist_dist(nodeAssort_mat(l_hemi,idx),nodeAssort_mat(r_hemi,idx));
    [~,mod_ks_assort(idx)] = simple_hist_dist(nodeAssort_mat_mod(l_hemi,idx),nodeAssort_mat_mod(r_hemi,idx));    
    
end

%% plot it

cmap = [0    0.4470    0.7410 ;
    0.8500    0.3250    0.0980];

subplot(1,3,1)
histogram(wsbm_ks_parti,'FaceColor',cmap(1,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
hold
histogram(mod_ks_parti,'FaceColor',cmap(2,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
axis square
title('Participation coefficient','FontWeight','normal')
xlabel(strcat('KS({\itbtwn. hemi.})'))
ylabel('Frequency')
xlim([ 0 0.5])

subplot(1,3,2)
histogram(wsbm_ks_zscore,'FaceColor',cmap(1,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
hold
histogram(mod_ks_zscore,'FaceColor',cmap(2,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
axis square
title('Within module z-score','FontWeight','normal')
xlabel(strcat('KS({\itbtwn. hemi.})'))
ylabel('Frequency')
xlim([ 0 0.5])

subplot(1,3,3)
histogram(wsbm_ks_assort,'FaceColor',cmap(1,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
hold
histogram(mod_ks_assort,'FaceColor',cmap(2,:),'FaceAlpha',0.6,'EdgeAlpha',0.01)
axis square
title('Assortativity','FontWeight','normal')
xlabel(strcat('KS({\itbtwn. hemi.})'))
ylabel('Frequency')
xlim([ 0 0.5])

legend('WSBM','Modular')
legend('boxoff')

hold off

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.8, 0.6]);

FIGURE_NAME = 'figL' ;

fileName = 'KS_hemi_stats.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%% WHOLE MAT VIZ

figure
%[~,sortIdx] = sort(comVecs.wsbm);
[xOnDiag,yOnDiag,sortIdx] = grid_communities(comVecs.wsbm);
h = imagesc(tmpAdj(sortIdx,sortIdx)) ;
set(h,'alphadata',tmpAdj(sortIdx,sortIdx)~=0);
hold on
p = plot(xOnDiag,yOnDiag,'Color',[0.1 0.1 0.1 0.1],'linewidth',0.5);
set(gca,'ytick',[])
set(gca,'xtick',[])
axis square
axis off
%hold off
h.AlphaData = h.AlphaData * 0.15 ;

cmap = brewermap(10,'Paired');
cmap = cmap(3:4,:);      

nodeAnnotData = [ ones([57 1]) ; ones([57 1]).*2 ] ;
nodeAnnotData = nodeAnnotData(sortIdx);

startMarks = [1 ; diff(nodeAnnotData)~=0 ] ;
numChunks = sum(startMarks) ;

startInd = find(startMarks);
endInd = [startInd(2:end)-1 ; length(startMarks)];

hold on

for i = 1:length(startInd)

    mn = startInd(i) - 0.5;
    mx = endInd(i) + 0.5;
    x = [0.5  0.5 ];
    y = [mn mx ];
    plot(x,y,'Color',cmap(nodeAnnotData(startInd(i)),:),'linewidth',10);             % plot community boundaries
    x = [ length(nodeAnnotData)+0.5 length(nodeAnnotData)+0.5 ] ;
    plot(y,x,'Color',cmap(nodeAnnotData(startInd(i)),:),'linewidth',10);
    
end

ll = findobj(gca,'Type','line');
legend([ ll(end-1) ll(end-3) ],'Left hemi.','Right hemi.')

hold off

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);

fileName = 'mat_hemi_example.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%% example hist...

% both...
figure
histogram(nodeAssort_mat(:,end),'EdgeAlpha',0.001,'FaceAlpha',0.25,...
    'FaceColor',[0.1 0.1 0.1])
axis square
title('Full brain node metric','FontWeight','normal')
ylabel('Frequency')
xlabel('Node metric')

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);

fileName = 'mat_hemi_example_both.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/','figL','/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

% left
figure
histogram(nodeAssort_mat(:,end),'EdgeAlpha',0.001,'FaceAlpha',0.05,...
    'FaceColor',[0.1 0.1 0.1])
hold on
histogram(nodeAssort_mat(l_hemi,end),'EdgeAlpha',0.01,'FaceAlpha',0.6,...
    'FaceColor',cmap(1,:))
axis square
title('Left hemi. node metric','FontWeight','normal')
ylabel('Frequency')
xlabel('Node metric')
legend('Full brain dist.','Left hemi. dist.')
hold off

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);

fileName = 'mat_hemi_example_left.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/','figL','/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

% right
figure
histogram(nodeAssort_mat(:,end),'EdgeAlpha',0.001,'FaceAlpha',0.05,...
    'FaceColor',[0.1 0.1 0.1])
hold on
histogram(nodeAssort_mat(r_hemi,end),'EdgeAlpha',0.01,'FaceAlpha',0.6,...
    'FaceColor',cmap(2,:))
axis square
title('Right hemi. node metric','FontWeight','normal')
legend('Full brain dist.','Right hemi. dist.')
ylabel('Frequency')
xlabel('Node metric')
hold off

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);

fileName = 'mat_hemi_example_right.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/','figL','/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%% ttests between the two 

[~,p,ci,stat] = ttest2(wsbm_ks_parti,mod_ks_parti,'Vartype','unequal') ;
mean(wsbm_ks_parti)
std(wsbm_ks_parti)

mean(mod_ks_parti)
std(mod_ks_parti)

[~,p,ci,stat] = ttest2(wsbm_ks_zscore,mod_ks_zscore,'Vartype','unequal') ;
mean(wsbm_ks_zscore)
std(wsbm_ks_zscore)

mean(mod_ks_zscore)
std(mod_ks_zscore)

[~,p,ci,stat] = ttest2(wsbm_ks_assort,mod_ks_assort,'Vartype','unequal') ;
mean(wsbm_ks_assort)
std(wsbm_ks_assort)

mean(mod_ks_assort)
std(mod_ks_assort)



