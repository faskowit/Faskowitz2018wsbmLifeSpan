%% no community

sp1 = subplot(1,2,1) ;
%pos = get(sp1,'position');
%set(sp1,'position',[pos(1:2)/4 pos(3:4)*2])

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

% PLOT THE ADJ
h = imagesc(E) ;
set(h,'alphadata',~isnan(E));
axis square
set(gca,'ytick',[])
set(gca,'xtick',[])

% PLOT THE NETWORK

node_annot = [ ones(group_sizes(1),1) .* 1 ;
               ones(group_sizes(2),1) .* 2 ;
               ones(group_sizes(3),1) .* 3 ;
               ones(group_sizes(4),1) .* 4 ];

sp2 = subplot(1,2,2) ;
 
E(isnan(E)) = 0;
g_obs = graph(Edg2Adj(E));
cmap = brewermap(length(unique(node_annot)),'PuOr');
H = plot(g_obs,'-',...
    'MarkerSize',8, ...               % node size in log scale
    'EdgeColor',[.8 .8 .8],...
    'EdgeAlpha',0.3,...
    'NodeCData',node_annot,...
    'NodeLabel',{});
layout(H,'force','Iterations',1000)

colormap(sp2,cmap)

axis square
set(gca,'ytick',[])
set(gca,'xtick',[])

%% modular

sp1 = subplot(1,2,1) ;

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);

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

% PLOT THE ADJ
h = imagesc(E) ;
set(h,'alphadata',~isnan(E));
axis square
set(gca,'ytick',[])
set(gca,'xtick',[])

% PLOT THE NETWORK

% get rid of the NaNs for the network view
E(isnan(E)) = 0;

node_annot = [ ones(group_sizes(1),1) .* 1 ;
               ones(group_sizes(2),1) .* 2 ;
               ones(group_sizes(3),1) .* 3 ;
               ones(group_sizes(4),1) .* 4 ];

sp2 = subplot(1,2,2) ;

g_obs = graph(Edg2Adj(E));
cmap = brewermap(length(unique(node_annot)),'PuOr');
H = plot(g_obs,'-',...
    'MarkerSize',8, ...               % node size in log scale
    'EdgeColor',[.6 .6 .6],...
    'EdgeAlpha',0.3,...
    'NodeCData',node_annot,...
    'NodeLabel',{});
layout(H,'force','Iterations',1000)

colormap(sp2,cmap)

axis square
set(gca,'ytick',[])
set(gca,'xtick',[])

%% core periphery

sp1 = subplot(1,2,1) ;

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);

R = [1,2,2,2;
     2,3,3,4;
     2,3,4,4;
     2,4,4,4] ;    

theta_w = [10,5; 10,1; 10,1; 10,1];
theta_e = [0.5; 0.35; 0.1; 0.1];
        
group_sizes = [25;25;25;25];

[E,~] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);

% make symmetric
E = Edg2Adj(E);
E = triu(E) + triu(E,1)';

% PLOT THE ADJ
h = imagesc(E) ;
set(h,'alphadata',~isnan(E));
axis square
set(gca,'ytick',[])
set(gca,'xtick',[])

% PLOT THE NETWORK

% get rid of the NaNs for the network view
E(isnan(E)) = 0;

node_annot = [ ones(group_sizes(1),1) .* 1 ;
               ones(group_sizes(2),1) .* 2 ;
               ones(group_sizes(3),1) .* 3 ;
               ones(group_sizes(4),1) .* 4 ];

sp2 = subplot(1,2,2) ;

g_obs = graph(Edg2Adj(E));
cmap = brewermap(length(unique(node_annot)),'PuOr');
H = plot(g_obs,'-',...
    'MarkerSize',8, ...               % node size in log scale
    'EdgeColor',[.6 .6 .6],...
    'EdgeAlpha',0.3,...
    'NodeCData',node_annot,...
    'NodeLabel',{});
layout(H,'force','Iterations',1000)

colormap(sp2,cmap)

axis square
set(gca,'ytick',[])
set(gca,'xtick',[])


%% dissasortative 

sp1 = subplot(1,2,1) ;

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);

R = [1,2,3,4;
     2,1,2,3;
     3,2,1,2;
     4,3,2,1] ;    

theta_w = [10,1; 10,5; 10,1; 10,1];
theta_e = [0.05; 0.35; 0.25; 0.15];
        
group_sizes = [25;25;25;25];

[E,~] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);

% make symmetric
E = Edg2Adj(E);
E = triu(E) + triu(E,1)';

% PLOT THE ADJ
h = imagesc(E) ;
set(h,'alphadata',~isnan(E));
axis square
set(gca,'ytick',[])
set(gca,'xtick',[]);

% PLOT THE NETWORK

% get rid of the NaNs for the network view
E(isnan(E)) = 0;

node_annot = [ ones(group_sizes(1),1) .* 1 ;
               ones(group_sizes(2),1) .* 2 ;
               ones(group_sizes(3),1) .* 3 ;
               ones(group_sizes(4),1) .* 4 ];

sp2 = subplot(1,2,2) ;

g_obs = graph(Edg2Adj(E));
cmap = brewermap(length(unique(node_annot)),'PuOr');
H = plot(g_obs,'-',...
    'MarkerSize',8, ...               % node size in log scale
    'EdgeColor',[.6 .6 .6],...
    'EdgeAlpha',0.3,...
    'NodeCData',node_annot,...
    'NodeLabel',{});
layout(H,'force','Iterations',1000)

colormap(sp2,cmap)

axis square
set(gca,'ytick',[])
set(gca,'xtick',[])


%% data we actually got...

tmpData = templateModel.Data.Raw_Data;

tmpData(isnan(tmpData)) = 0;

coords = load('data/external/yeo114_coords.mat') ;
coords = coords.coords ;

% node_annot = [ ones(group_sizes(1),1) .* 1 ;
%                ones(group_sizes(2),1) .* 2 ;
%                ones(group_sizes(3),1) .* 3 ;
%                ones(group_sizes(4),1) .* 4 ];

[~,node_annot] = community_assign(templateModel) ;

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




