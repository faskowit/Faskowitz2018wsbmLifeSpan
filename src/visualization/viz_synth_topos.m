%% no community

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);

R = [1,2,3,4;
     1,2,3,4;
     1,2,3,4;
     1,2,3,4] ;    

theta_w = [10,1; 10,1; 10,1; 10,1];
theta_e = [0.33; 0.33; 0.33; 0.33];
group_sizes = [25;25;25;25];
[E,~] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);
E = Edg2Adj(E);
E = triu(E,0) + triu(E,1)';
E(isnan(E)) = 0;
E(E<0) = 0;

node_annot = [ ones(group_sizes(1),1) .* 1 ;
               ones(group_sizes(2),1) .* 2 ;
               ones(group_sizes(3),1) .* 3 ;
               ones(group_sizes(4),1) .* 4 ];

viz_plot3Way_tight(E,R,theta_e,node_annot)

%% modular

% SET UP DATA
R = [1,2,3,4;
     2,1,2,3;
     3,2,1,2;
     4,3,2,1] ;    

theta_w = [10,5; 1,1; 1,1; 1,1];
theta_e = [0.5; 0.1; 0.1; 0.1];
group_sizes = [25;25;25;25];
[E,~] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);

% make symmetric
E = Edg2Adj(E);
E = triu(E) + triu(E,1)';
E(isnan(E)) = 0;
E(E<0) = 0;

node_annot = [ ones(group_sizes(1),1) .* 1 ;
               ones(group_sizes(2),1) .* 2 ;
               ones(group_sizes(3),1) .* 3 ;
               ones(group_sizes(4),1) .* 4 ];

           
viz_plot3Way_tight(E,R,theta_e,node_annot)

%% core periphery

R = [1,2,2,2;
     2,3,3,4;
     2,3,4,4;
     2,4,4,4] ;    

theta_w = [10,5; 10,1; 10,1; 10,1];
theta_e = [0.75; 0.33; 0.1; 0.1];       
group_sizes = [25;25;25;25];
[E,~] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);

% make symmetric
E = Edg2Adj(E);
E = triu(E) + triu(E,1)';
E(isnan(E)) = 0;
E(E<0) = 0;

viz_plot3Way_tight(E,R,theta_e,node_annot)

%% dissasortative 

R = [1,2,3,4;
     2,1,2,3;
     3,2,1,2;
     4,3,2,1] ;    

theta_w = [10,1; 10,5; 10,1; 10,1];
theta_e = [0.05; 0.33; 0.25; 0.15];    
group_sizes = [25;25;25;25];
[E,~] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);

% make symmetric
E = Edg2Adj(E);
E = triu(E) + triu(E,1)';
E(isnan(E)) = 0;
E(E<0) = 0;

viz_plot3Way_tight(E,R,theta_e,node_annot)

%% data we actually got...

tmpData = templateModel.Data.Raw_Data;
%tmpData(isnan(tmpData)) = 0;

% maybe lets get one hemi and fit a blockmodel...
tmpData = tmpData(1:57,1:57);
[~,hypotheticalModel] = wsbm(tmpData, ...
    4, ...
    'W_Distr', 'normal', ...
    'E_Distr', 'bernoulli', ...
    'NumTrials', 25,...
    'verbosity', 1);

% coords = load('data/external/yeo114_coords.mat') ;
% coords = coords.coords ;

% node_annot = [ ones(group_sizes(1),1) .* 1 ;
%                ones(group_sizes(2),1) .* 2 ;
%                ones(group_sizes(3),1) .* 3 ;
%                ones(group_sizes(4),1) .* 4 ];

tmpData(isnan(tmpData)) = 0;
[~,node_annot] = community_assign(hypotheticalModel) ;
[~,sortIdx] = sort(node_annot);
vvv = viz_plot3Way_tight(tmpData(sortIdx,sortIdx),hypotheticalModel.R_Struct.R,...
    hypotheticalModel.Para.predict_e,node_annot) ;

axes(vvv(2)) 

[X,Y,INDSORT] = grid_communities(node_annot); % call function
hold on;                                 % hold on to overlay community visualization
plot(X,Y,'Color',[0.1 0.1 0.1 0.1],'linewidth',2);             % plot community boundaries

%%


figure
g_obs = graph(tmpData);
cmap = brewermap(length(unique(node_annot)),'Paired');
colormap(cmap)
H = plot(g_obs,'-',...
    'MarkerSize',8, ...               % node size in log scale
    'EdgeColor',[.6 .6 .6],...
    'EdgeAlpha',0.3,...
        'NodeCData',node_annot,...
    'NodeLabel',{}) ;
   % 'XData',coords(:,1),'YData',coords(:,2));
layout(H,'force','Iterations',1000)

axis square
set(gca,'ytick',[])
set(gca,'xtick',[])

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);








