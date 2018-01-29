function [ subp ] =  viz_plot3Way_tight(Edata,Rdata,thetaEdata,nodeAnnotData)
           
% setup plot       
%sp1 = subplot(1,3,1) ;
subp = tight_subplot(1,3,[.01 .03],[.1 .01],[.01 .01]) ;

% make it full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);

%% PLOT THE NETWORK

axes(subp(1))

% get rid of the NaNs for the network view
g_obs = graph(Edg2Adj(Edata));

% set a colormap for this pannel
cmap = brewermap(length(unique(nodeAnnotData)),'Spectral');

%line widths for better viz
LWidths = abs(2*g_obs.Edges.Weight/max(g_obs.Edges.Weight));

Dsizes = (8*g_obs.degree ./ max(g_obs.degree)) + 4;

H = plot(g_obs,'-',...
    'MarkerSize',Dsizes, ...     
    'EdgeColor',[.6 .6 .6],...
    'EdgeAlpha',0.3,...
    'LineWidth',LWidths,...
    'NodeCData',nodeAnnotData,...
    'NodeLabel',{});
layout(H,'force','Iterations',1000)

colormap(subp(1),cmap) ;
cb1 = colorbar();
set( cb1, 'YDir', 'reverse' );
cb1.Ticks = [] ;
ylabel(cb1,'Communities','FontSize',18);

axis square
set(gca,'ytick',[])
set(gca,'xtick',[])

%% PLOT THE ADJ

%sp2 = subplot(1,3,2) ;

axes(subp(2))

% PLOT THE ADJ
h = imagesc(Edata) ;
%colormap('parula')

set(h,'alphadata',(Edata > 0));
axis square
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'visible','off')

colormap(flipud(gray))
cb2 = colorbar() ;
cb2.Ticks = [] ;
ylabel(cb2,'Edge Strength','FontSize',18);
%colorbar(cb2,'off')

nc = max(nodeAnnotData);
c = sort(nodeAnnotData);

hold on

for i = 1:nc
    ind = find(c == i);
    if ~isempty(ind)
        mn = min(ind) - 0.5;
        mx = max(ind) + 0.5;
        x = [0.5  0.5 ];
        y = [mn mx ];
        %xGrid = [xGrid, x]; 
        %yGrid = [yGrid, y];
        plot(x,y,'Color',cmap(i,:),'linewidth',6);             % plot community boundaries
        %plot(y,x,'Color',cmap(i,:),'linewidth',8);             % plot community boundaries
        x = [ length(nodeAnnotData)+0.5 length(nodeAnnotData)+0.5 ] ;
        plot(y,x,'Color',cmap(i,:),'linewidth',6);
    end
end

%hold on;                                 % hold on to overlay community visualization
%plot(xGrid,yGrid,'r','linewidth',5);             % plot community boundaries
%hold off;

%% PLOT the affinity mat

%sp3 = subplot(1,3,3) ;
axes(subp(3))

% make edge probability mat
mat = zeros(length(unique(nodeAnnotData)));
for idx=1:(length(unique(nodeAnnotData)))
    for jdx=1:(length(unique(nodeAnnotData)))
        mat(idx,jdx) = thetaEdata(Rdata(idx,jdx));
    end
end

im = imagesc(mat);
caxis([0 1])

textStrings = num2str(mat(:),'%.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
%wantInd = mat>0; % only write text for non-zero vals
%textStrings(wantInd == 0) = {''} ; % set the 0's to null
[x,y] = meshgrid(1:(size(mat,1)));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');   

%midValue = 0.4 ; % set to make it easier to read
textColors = repmat(mat(:)> 0.25,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2),'FontSize',15);  %# Change the text colors

cmap3 = brewermap(100,'Purples');

colormap(subp(3),cmap3)
cb3 = colorbar() ;
cb3.Ticks = [] ;
ylabel(cb3,'Edge Exist Probability','FontSize',18);

axis square
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'visible','off')

gridLen = 1:(length(unique(nodeAnnotData)));

nc = max(gridLen);
c = sort(gridLen);

hold on

for i = 1:nc
    ind = find(c == i);
    if ~isempty(ind)
        mn = min(ind) - 0.5;
        mx = max(ind) + 0.5;
        x = [0.5  0.5 ];
        y = [mn mx ];
        %xGrid = [xGrid, x]; 
        %yGrid = [yGrid, y];
        plot(x,y,'Color',cmap(i,:),'linewidth',6);             % plot community boundaries
        %plot(y,x,'Color',cmap(i,:),'linewidth',8);             % plot community boundaries
        x = [ (length(unique(nodeAnnotData)))+0.5 (length(unique(nodeAnnotData)))+0.5 ] ;
        plot(y,x,'Color',cmap(i,:),'linewidth',6);
    end
end

%hold on;                                 % hold on to overlay community visualization
%plot(xGrid,yGrid,'r','linewidth',5);             % plot community boundaries
%hold off;

