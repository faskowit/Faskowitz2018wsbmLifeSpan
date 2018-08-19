clc 
clearvars

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_templateModel_1.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
load(loadName) ;

% % the stats maps for wsbm, mod, yeo
% loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_comVec_static_results.mat');
% load(loadName) ;

% the stats maps for wsbm, mod, yeo
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_comVec_static_regResults.mat');
load(loadName) ;

FIGURE_NAME = 'figE' ;

outputdir = strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/');
mkdir(outputdir)

%% figure out the fits that are best

all_reg_results = cell([6 1]) ;
all_reg_results{1} = regResultsCov_cos ;
all_reg_results{2} = regResultsCovWMov_cos ;
all_reg_results{3} = regResultsNoCov_cos ;
all_reg_results{4} = regResultsCov_cb ;
all_reg_results{5} = regResultsCovWMov_cb ;
all_reg_results{6} = regResultsNoCov_cb ;

reg_result_names = {'Cos sim. w/ covar.' 'Cos sim. w/ covar+mov' 'Cos similarity' ...
    'CB dist. w/ covar' 'CB dist. w/ covar+mov' 'CB distance'};

reg_ylabel_names = { 'Cos similarity adj. for covar' 'Cos similarity adj. for covar, mov' 'Cos similarity' ...
    'CB distance adj. for covar' 'CB distance adj. for covar, mov' 'CB distance'} ;

default_cmap = [0    0.4470    0.7410 ;
                0.8500    0.3250    0.0980;
                0.9290    0.6940    0.1250 ] ;

subp = tight_subplot(2,3,[.15 .1],[.1 .05],[.15 .1]) ;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);      

for idx = 1:length(all_reg_results)
    
    currentRegResults = all_reg_results{idx} ;

    [~,wsbm_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults{1}));
    [~,mod_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults{2}));
%     [~,yeo_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults{3}));

    wsbm_trend = currentRegResults{1}{wsbm_minIdx};
    mod_trend = currentRegResults{2}{mod_minIdx};
%     yeo_trend = currentRegResults{3}{yeo_minIdx};
    
%% plot it

    %fig1 = figure ;
    axes(subp(idx)) 
    
%     % yeo
%     pl = viz_blockRegress(yeo_trend,0,default_cmap(3,:));
%     set(pl,'MarkerFaceColor', default_cmap(3,:),...
%         'MarkerEdgeColor',default_cmap(3,:),...
%         'MarkerFaceAlpha',.15,...
%         'MarkerEdgeAlpha',.1) 
%     hold

    % mod
    pl = viz_blockRegress(mod_trend,0,default_cmap(2,:));
    set(pl,'MarkerFaceColor', default_cmap(2,:),...
        'MarkerEdgeColor',default_cmap(2,:),...
        'MarkerFaceAlpha',.15,...
        'MarkerEdgeAlpha',.1) 
    hold

    % wsbm
    pl = viz_blockRegress(wsbm_trend,0, default_cmap(1,:)) ;
    set(pl,'MarkerFaceColor',default_cmap(1,:),...
        'MarkerFaceColor',default_cmap(1,:),...
        'MarkerFaceAlpha',.15,...
        'MarkerEdgeAlpha',.1)
    hold
   
    xlim([min(pl.XData)-3 max(pl.XData)+3])
    
    ylabel(reg_ylabel_names{idx})
    xlabel('Age')
    
    title(reg_result_names{idx})
    
    yrange = ylim ;
    yrangeAbs = yrange(2) - yrange(1) ;
    % add the R2!! 
    if idx < 4
        ypos = yrangeAbs * 0.25;
        ypos = yrange(1) + ypos ;
    else
        ypos = yrangeAbs * 0.95;
        ypos = yrange(1) + ypos ;
    end

    set(gca, 'FontSize', 14)
    
    annotText = { strcat('WSBM R^2:',32,num2str(round(wsbm_trend.xvalR2 / 100,3))) ...
        strcat('Mod R^2:',32,num2str(round(mod_trend.xvalR2 / 100,3))) ...
%         strcat('yeo R^2:',32,num2str(round(yeo_trend.xvalR2 / 100,3))) ...
        } ;

    text((min(pl.XData)+2),ypos, annotText,'FontSize',16,'VerticalAlignment','cap')    
    
end

% % save it
fileName = 'comVec_static_distances.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%% but also just visualize the similarity and distance wihtout covar

figure
subp = tight_subplot(1,2,[.10 .10],[.1 .05],[.15 .1]) ;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);      

for idx = [ 3 6 ]
% for idx = 1:length(all_reg_results)
    
    currentRegResults = all_reg_results{idx} ;

    [~,wsbm_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults{1}));
    [~,mod_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults{2}));
%     [~,yeo_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults{3}));

    wsbm_trend = currentRegResults{1}{wsbm_minIdx};
    mod_trend = currentRegResults{2}{mod_minIdx};
%     yeo_trend = currentRegResults{3}{yeo_minIdx};
    %fig1 = figure ;
    
    % minus for because we indexing 3-6
    axes(subp(idx / 3)) 
    
%     % yeo
%     [ pl, ln1] = viz_blockRegress(yeo_trend,0,default_cmap(3,:));
%     set(pl,'MarkerFaceColor', default_cmap(3,:),...
%         'MarkerEdgeColor',default_cmap(3,:),...
%         'MarkerFaceAlpha',.15,...
%         'MarkerEdgeAlpha',.1) 
%     hold

    % mod
    [ pl , ln2 ] = viz_blockRegress(mod_trend,0,default_cmap(2,:));
    set(pl,'MarkerFaceColor', default_cmap(2,:),...
        'MarkerEdgeColor',default_cmap(2,:),...
        'MarkerFaceAlpha',.15,...
        'MarkerEdgeAlpha',.1) 
    hold

    % wsbm
    [ pl , ln3 ] = viz_blockRegress(wsbm_trend,0, default_cmap(1,:)) ;
    set(pl,'MarkerFaceColor',default_cmap(1,:),...
        'MarkerFaceColor',default_cmap(1,:),...
        'MarkerFaceAlpha',.15,...
        'MarkerEdgeAlpha',.1)
    hold
   
    xlim([min(pl.XData)-3 max(pl.XData)+3])
    
    ylabel(reg_ylabel_names{idx})
    xlabel('Age')
    
    title(reg_result_names{idx})
    
    yrange = ylim ;
    yrangeAbs = yrange(2) - yrange(1) ;
    % add the R2!! 
    
    if idx < 4
        ypos = yrangeAbs * 0.15;
        ypos = yrange(1) + ypos ;
%         ll = legend([ln3 ln2 ], {'WSBM' 'Modular'},'Location','SouthEast','FontSize',16) ;
% 
%         set(ll,'Units','inches')
%                 legend('boxoff')
    else
        ypos = yrangeAbs * 0.97;
        ypos = yrange(1) + ypos ;
%         ll = legend([ln3 ln2 ln1 ], {'WSBM' 'Modular' 'Yeo'},'Location','NorthEast','FontSize',12);
%         ll = legend([ln3 ln2 ], {'WSBM' 'Modular'},'Location','NorthEast','FontSize',16) ;
%         set(ll,'Units','inches')
%         legend('boxoff')
    end

    set(gca, 'FontSize', 16) 
    
    annotText = { strcat('WSBM R^2:',32,num2str(round(wsbm_trend.xvalR2 / 100,3))) ...
        strcat('Modular R^2:',32,num2str(round(mod_trend.xvalR2 / 100,3))) ...
%         strcat('Yeo R^2:',32,num2str(round(yeo_trend.xvalR2 / 100,3))) ...
        } ;

    text((min(pl.XData)),ypos, annotText,'FontSize',18,'VerticalAlignment','cap')    
        
end

% save it
fileName = 'comVec_static_distances_noCov.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%% extra viz for communities...

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_comVec_static_results.mat');
load(loadName) ;

y_lim = [0.5 1] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WSBM

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.55]);      
for idx=1:10
   subplot(2,5,idx) 
%    plot(fitlm(datasetDemo.age,wsbm_comms_weiVec_cos(idx,:),'quadratic'));
    scatter(datasetDemo.age,wsbm_comms_weiVec_cos(idx,:),24,default_cmap(1,:)) ;
   %ylim([0 0.5])
   ylim(y_lim)
      
   title(strcat('community',32,num2str(idx)));
   
end
suptitle('wsbm commumity cos similarity')

% save
fileName = 'comVec_static_comm_wsbm_cos.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOD

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.55]);      
for idx=1:10
   subplot(2,5,idx) 
%    plot(fitlm(datasetDemo.age,mod_comms_weiVec_cos(idx,:),'quadratic'));
    scatter(datasetDemo.age,mod_comms_weiVec_cos(idx,:),24,default_cmap(2,:))
   %ylim([0 0.5])
   ylim(y_lim)
         
   title(strcat('community',32,num2str(idx)));
   
end
suptitle('mod community cos similarity')

% save
fileName = 'comVec_static_comm_mod_cos.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YEO

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.55]);      
for idx=1:7 % change for YEO
   subplot(2,5,idx) 
%    plot(fitlm(datasetDemo.age,mod_comms_weiVec_cos(idx,:),'quadratic'));
    scatter(datasetDemo.age,yeo_comms_weiVec_cos(idx,:),24,default_cmap(3,:))
   %ylim([0 0.5])
   ylim(y_lim)
         
   title(strcat('community',32,num2str(idx)));
   
end
suptitle('yeo community cos similarity')

% save
fileName = 'comVec_static_comm_yeo_cos.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAIRWISE

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.55]);      
subp = tight_subplot(2,5,[.15 .05],[.1 .001],[.1 .05]) ;
for idx=1:10
 
    axes(subp(idx))
    
%    plot(fitlm(datasetDemo.age,mod_comms_weiVec_cos(idx,:),'quadratic'));
%    sc = scatter(wsbm_comms_weiVec_cos(idx,:),mod_comms_weiVec_cos(idx,:),24,'k') ;
   sc = scatter(wsbm_comms_weiVec_cos(idx,:),mod_comms_weiVec_cos(idx,:),24,datasetDemo.age) ;
   colormap(brewermap(50,'RdYlGn'))
   %ylim([0 0.5])
   %ylim(y_lim)
   axis square
   
   xlim([ min([min(sc.XData) min(sc.YData)]) 1 ])
   ylim([ min([min(sc.XData) min(sc.YData)]) 1 ])
   
   xlabel('wsbm')
   ylabel('mod')
   
   refline(1,0)
  
   title(strcat('community',32,num2str(idx)));
   
end
suptitle('wsbm vs mod community cos similarity')

% save
fileName = 'comVec_static_comm_paired_cos.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

figure ;

sc = scatter(datasetDemo.age,datasetDemo.age,24,datasetDemo.age) ;
cb = colorbar()
colormap(brewermap(50,'RdYlGn'))
set(gca,'Visible','off')

fileName = 'comVec_static_comm_paired_cos_colorbar.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%% CB distance by community too

y_lim = [0 0.75] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WSBM

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.55]);      
for idx=1:10
   subplot(2,5,idx) 
   scatter(datasetDemo.age,wsbm_comms_weiVec_cb(idx,:),24,default_cmap(1,:)) ;
   ylim(y_lim)
      
   title(strcat('community',32,num2str(idx)));
   
end
suptitle('wsbm commumity CB distance')

% save
fileName = 'comVec_static_comm_wsbm_cb.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOD

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.55]);      
for idx=1:10
   subplot(2,5,idx) 
   scatter(datasetDemo.age,mod_comms_weiVec_cb(idx,:),24,default_cmap(2,:))
   ylim(y_lim)
         
   title(strcat('community',32,num2str(idx)));
   
end
suptitle('mod community CB distance')

% save
fileName = 'comVec_static_comm_mod_cb.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YEO

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.55]);      
for idx=1:7 % change for YEO
   subplot(2,5,idx) 
%    plot(fitlm(datasetDemo.age,mod_comms_weiVec_cos(idx,:),'quadratic'));
    scatter(datasetDemo.age,yeo_comms_weiVec_cb(idx,:),24,default_cmap(3,:))
   %ylim([0 0.5])
   ylim(y_lim)
         
   title(strcat('community',32,num2str(idx)));
   
end
suptitle('yeo community CB distance')

% save
fileName = 'comVec_static_comm_yeo_cb.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAIRWISE

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.55]);      
subp = tight_subplot(2,5,[.15 .05],[.1 .001],[.1 .05]) ;
for idx=1:10
 
    axes(subp(idx))

   sc = scatter(wsbm_comms_weiVec_cb(idx,:),mod_comms_weiVec_cb(idx,:),24,datasetDemo.age) ;
   colormap(brewermap(50,'RdYlGn'))
   
   axis square
   
   xlim([ min([min(sc.XData) min(sc.YData)]) 1 ])
   ylim([ min([min(sc.XData) min(sc.YData)]) 1 ])
   
   xlabel('wsbm')
   ylabel('mod')
   
   refline(1,0)
  
   title(strcat('community',32,num2str(idx)));
   
end
suptitle('wsbm vs mod community CB distancre')

% save
fileName = 'comVec_static_comm_paired_cb.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

figure ;

sc = scatter(datasetDemo.age,datasetDemo.age,24,datasetDemo.age) ;
cb = colorbar()
colormap(brewermap(50,'RdYlGn'))
set(gca,'Visible','off')

fileName = 'comVec_static_comm_paired_cb_colorbar.png';
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)