
clc 
clearvars

config_file='config_scale125.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_templateModel_1.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
load(loadName) ;

FIGURE_NAME = 'figA' ;

outputdir = strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/');
mkdir(outputdir)

% lets make a nice picture of the template model

%% data stuff

nNodes = templateModel.Data.n ;
nComm = templateModel.R_Struct.k ; 

% get the template data
templateData = templateModel.Data.Raw_Data ;
templateData(isnan(templateData)) = 0;

% scramble
scrmb = randperm(nNodes);
templateData = templateData(scrmb,scrmb) ;
%comLabels = comVecs.wsbm(scrmb) ;


%% plot main adjacecny matrix
% adjacency matrix plot with on-diagonal communities outlined in bold red
% lines and off diagonal communities marked by red lines
% get community assignment stuff

% parcels = {'wsbm' 'mod' 'yeo' };
parcels = {'wsbm' 'mod' };
parcelName = {'WSBM' 'Modular'} ;

for fig = 1:length(parcels)
%for fig = 2
 
    figure
    
    comLabels = comVecs.((lower(parcels{fig}))) ;    
    comLabels = CBIG_HungarianClusterMatch(comVecs.wsbm,comLabels);
    comLabels = comLabels(scrmb);
    
    nComm = length(unique(comLabels)) ;
    uniqueLab = unique(comLabels) ;
    
    [xOnDiag,yOnDiag,sortIdx] = grid_communities(comLabels);

    % code from WSBM stuffs
    A_sort = zeros(nNodes);
    list = zeros(1,nNodes);
    breaks = zeros(1,nComm);
    cur = 1;
    for idx = 1:nComm
        indicies = find(comLabels == uniqueLab(idx));
        
        if isempty(indicies)
         
            list(cur) = [] ; %list(cur:cur+length(indicies)-1) = indicies;
            cur = 0 ;%cur + length(indicies);
            breaks(idx) = cur ;%cur-1;
            
        end
        
        list(cur:cur+length(indicies)-1) = indicies;
        cur = cur + length(indicies);
        breaks(idx) = cur-1;
    end
    for idx = 1:nNodes
        A_sort(idx,:) = templateData(list(idx),list);
    end

    %Plot the Matrix
    h = imagesc(A_sort,[min(A_sort(:))-.00001,max(A_sort(:))]);
    set(h,'alphadata',(A_sort > 0) .* 0.9);

    % fix image properties 
    ax = gca ;
    axis square
    % axis([0.5 (nNodes+0.5) 0.5 (nNodes+0.5)]); 
    % ax.Box = 'on' ;
    % set(ax,'Ydir','reverse');
    % set(ax,'ytick',[])
    % set(ax,'xtick',[])

    hold on

    % plot off diagonal
    for idx = 1:(nComm-1)

        lineWidth = 1.5;
        offDiagColor = [1 0 0 0.5] ; 

        % vertical   
        plot([breaks(idx)+0.5,breaks(idx)+0.5],[breaks(idx+1)+0.5,breaks(nComm)+.5],...
            'Color',offDiagColor,'LineWidth',lineWidth);
        if idx > 1
            plot([breaks(idx)+0.5,breaks(idx)+0.5],[-0.5,breaks(idx-1)+0.5],...
                'Color',offDiagColor,'LineWidth',lineWidth);
        end

        % horizontal
        plot([breaks(idx+1)+0.5,breaks(nComm)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
            'Color',offDiagColor,'LineWidth',lineWidth);
        if idx > 1
            plot([-.5,breaks(idx-1)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
                'Color',offDiagColor,'LineWidth',lineWidth);
        end

    end      

    %plot on diagonal 
    pp = plot(xOnDiag,yOnDiag,'Color',[1 0 0 1],'linewidth',2.5);

    % compute some yticks
    breaks2 = [ 0 breaks ] ;
    midlabelpoint = zeros([nComm 1]);
    for idx = 1:length(breaks)
        midlabelpoint(idx) = floor( (breaks2(idx+1) - breaks2(idx)) / 2) + breaks2(idx);  
    end

    set(ax,'xtick',midlabelpoint)
    set(ax,'xticklabel',uniqueLab)
    set(ax,'ticklength',[ 0 0]) 

    set(ax,'ytick',midlabelpoint)
    set(ax,'yticklabel',uniqueLab)
    set(ax,'ticklength',[ 0 0]) 

    % add the colorbar
    cb = colorbar('peer',ax);
    cb.Label.String = 'Streamline Density' ;
    cb.Label.FontSize = 12 ;
    cb.Label.FontName = 'Arial';

    set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.4, 0.8]);

    ax.Title.String = { parcelName{fig}, ' community structure'};
    ax.TitleFontSizeMultiplier = 1.5 ;

    % save it
    fileName = strcat(parcels{fig},'_adj.png');
    ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
    %set(gcf,'paperpositionmode','auto');
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
    
end

%% plot the predicted edge existence and weight

figure; 

nNodes = templateModel.Data.n ;
nComm = templateModel.R_Struct.k ; 

tsubp = tight_subplot(2,1,0.07) ;
axes(tsubp(1))

%Plot the Matrix
h = imagesc(make_square(templateModel.Para.predict_e));
set(h,'alphadata',make_square(templateModel.Para.predict_e > 0))

% fix image properties 
ax = gca ;
axis square
axis([0.5 (nComm+0.5) 0.5 (nComm+0.5)]); 
ax.Box = 'on' ;
%set(ax,'Ydir','reverse');
set(ax,'ytick',1:nComm)
set(ax,'xtick',1:nComm)

% add the colorbar
cb = colorbar('peer',ax);
cb.Label.String = 'Predicted existence' ;
cb.Label.FontSize = 12 ;
cb.Label.FontName = 'Arial';

% weights now

axes(tsubp(2))

%Plot the Matrix
h = imagesc(make_square(templateModel.Para.predict_w));
set(h,'alphadata',make_square(templateModel.Para.predict_w > 0))

% fix image properties 
ax = gca ;
axis square
axis([0.5 (nComm+0.5) 0.5 (nComm+0.5)]); 
ax.Box = 'on' ;
%set(ax,'Ydir','reverse');
set(ax,'ytick',1:nComm)
set(ax,'xtick',1:nComm)

% add the colorbar
cb = colorbar('peer',ax);
cb.Label.String = 'Predicted weight' ;
cb.Label.FontSize = 12 ;
cb.Label.FontName = 'Arial';

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.4, 0.8]);

% save it
fileName = strcat('wsbm_predicted_e_w.png');
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

%% module colors

nComm4ColorMap = templateModel.R_Struct.k ;

for fig = 1:length(parcels)
%for fig = 3
 
    figure
    
    comLabels = comVecs.((lower(parcels{fig}))) ;    
    comLabels = CBIG_HungarianClusterMatch(comVecs.wsbm,comLabels);
    comLabels = comLabels(scrmb);
    
    nComm = length(unique(comLabels)) ;
    uniqueLab = unique(comLabels) ;
    
    [xOnDiag,yOnDiag,sortIdx] = grid_communities(comLabels);

    % code from WSBM stuffs
    A_sort = zeros(nNodes);
    list = zeros(1,nNodes);
    breaks = zeros(1,nComm);
    cur = 1;
    for idx = 1:nComm
        indicies = find(comLabels == uniqueLab(idx));
        
        if isempty(indicies)
         
            list(cur) = [] ; %list(cur:cur+length(indicies)-1) = indicies;
            cur = 0 ;%cur + length(indicies);
            breaks(idx) = cur ;%cur-1;
            
        end
        
        list(cur:cur+length(indicies)-1) = indicies;
        cur = cur + length(indicies);
        breaks(idx) = cur-1;
    end
    for idx = 1:nNodes
        A_sort(idx,:) = templateData(list(idx),list);
    end

    %Plot the Matrix
    h = imagesc(A_sort,[min(A_sort(:))-.00001,max(A_sort(:))]);
    set(h,'alphadata',(A_sort > 0) .* 0.15);

    % fix image properties 
    ax = gca ;
    axis square
    % axis([0.5 (nNodes+0.5) 0.5 (nNodes+0.5)]); 
    % ax.Box = 'on' ;
    % set(ax,'Ydir','reverse');
    % set(ax,'ytick',[])
    % set(ax,'xtick',[])

    hold on

    % plot off diagonal
    for idx = 1:(nComm-1)

        lineWidth = 1.5;
        offDiagColor = [1 0 0 0.25] ; 

        % vertical   
        plot([breaks(idx)+0.5,breaks(idx)+0.5],[breaks(idx+1)+0.5,breaks(nComm)+.5],...
            'Color',offDiagColor,'LineWidth',lineWidth);
        if idx > 1
            plot([breaks(idx)+0.5,breaks(idx)+0.5],[-0.5,breaks(idx-1)+0.5],...
                'Color',offDiagColor,'LineWidth',lineWidth);
        end

        % horizontal
        plot([breaks(idx+1)+0.5,breaks(nComm)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
            'Color',offDiagColor,'LineWidth',lineWidth);
        if idx > 1
            plot([-.5,breaks(idx-1)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
                'Color',offDiagColor,'LineWidth',lineWidth);
        end

    end      

    % MAKE THIS TEN...
    cmap_mod = brewermap(nComm4ColorMap,'paired') ;

    for idx=0:(nComm-1)

        plot(xOnDiag( (idx*6)+1:((idx+1)*6) ) ,...
            yOnDiag(  (idx*6)+1:((idx+1)*6) ) ,...
            'Color',cmap_mod(uniqueLab(idx+1),:),'linewidth',3.5)
    end

    % compute some yticks
    breaks2 = [ 0 breaks ] ;
    midlabelpoint = zeros([nComm 1]);
    for idx = 1:length(breaks)
        midlabelpoint(idx) = floor( (breaks2(idx+1) - breaks2(idx)) / 2) + breaks2(idx);  
    end

    set(ax,'xtick',midlabelpoint)
    set(ax,'xticklabel',uniqueLab)
    set(ax,'ticklength',[ 0 0]) 

    set(ax,'ytick',midlabelpoint)
    set(ax,'yticklabel',uniqueLab)
    set(ax,'ticklength',[ 0 0]) 

    set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.4, 0.8]);

    ax.Title.String = { parcelName{fig}, ' community colors'};
    ax.TitleFontSizeMultiplier = 1.5 ;

    % save it
    fileName = strcat(parcels{fig},'_commView.png');
    ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
    %set(gcf,'paperpositionmode','auto');
    print(gcf,'-dpng','-r500',ff);
    close(gcf)
    
end

%% plot pred w & e
% as a zscore plot

pred_w = templateModel.Para.predict_w ;
pred_w(isnan(pred_w)) = 0 ;

z_predict_e = zscore(templateModel.Para.predict_e);
z_predict_w = zscore(pred_w);

scat = scatter(z_predict_e,z_predict_w) ;
%axis square

hold

plot([ 0 0 ], get(gca,'ylim'),'k','Color',[0 0 0 0.25])
plot(get(gca,'xlim'), [0 0],'k','Color',[0 0 0 0.25])

xlabel('Z-score predicted edge')
ylabel('Z-score predicted weight')

set(gcf, 'Units', 'Normalized', 'Position', [0.2, 0.2, 0.4, 0.4]);
pbaspect([1 1 1])
tightfig()

% ttt = lsline()
beta = polyfit(z_predict_e,z_predict_w,1);
linLine = refline(gca,beta);
linLine.Color = [ linLine.Color 0.25 ] ;

% save it
fileName = strcat('wsbm_predicted_e_vs_w.png');
ff = fullfile(strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/',OUTPUT_STR,'_',fileName)); 
%set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r500',ff);
close(gcf)

