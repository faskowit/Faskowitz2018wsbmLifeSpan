% lets make a nice picture of the template model

figure;

% get community assignment stuff
[~,templateModelca] = community_assign(templateModel);
[xOnDiag,yOnDiag,sortIdx] = grid_communities(templateModelca);

% code from WSBM stuffs
n = length(templateModelca);
k = length(unique(templateModelca));
A_sort = zeros(n);
list = zeros(1,n);
breaks = zeros(1,k);
cur = 1;
for idx = 1:k
    indicies = find(templateModelca == idx);
    list(cur:cur+length(indicies)-1) = indicies;
    cur = cur + length(indicies);
    breaks(idx) = cur-1;
end
for idx = 1:n
    A_sort(idx,:) = data(list(idx),list);
end

%Plot the Matrix
h = imagesc(A_sort,[min(A_sort(:))-.00001,max(A_sort(:))]);
set(h,'alphadata',~isnan(A_sort));
%h.AlphaData = h.AlphaData * 0.25 ;

cb = colorbar;
%ylabel(cb,'Streamline Density')
cb.Label.String = 'Streamline Density' ;
cb.Label.FontSize = 14 ;
cb.Label.FontName = 'Arial';

hold all; 

% plot off diagonal
for idx = 1:(k-1)
    
    lineWidth = 1;
    %onDiagColor = [1 0 0.75 0.25] ;
    offDiagColor = [1 0 0 0.35] ; 
    
    
    % vertical   
    plot([breaks(idx)+0.5,breaks(idx)+0.5],[breaks(idx+1)+0.5,breaks(k)+.5],...
        'Color',offDiagColor,'LineWidth',lineWidth);
    if idx > 1
        plot([breaks(idx)+0.5,breaks(idx)+0.5],[-0.5,breaks(idx-1)+0.5],...
            'Color',offDiagColor,'LineWidth',lineWidth);
    end

    
    % horizontal
    plot([breaks(idx+1)+0.5,breaks(k)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
        'Color',offDiagColor,'LineWidth',lineWidth);
    if idx > 1
        plot([-.5,breaks(idx-1)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
            'Color',offDiagColor,'LineWidth',lineWidth);
    end
    
end      

%plot on diagonal 
pp = plot(xOnDiag,yOnDiag,'Color',[1 0 0 1],'linewidth',2.5);
%uistack(pp,'top');

% % lets try to plot on diagonal with lines according to community
% % assignments
% 
% cmap_mod = brewermap(k,'paired') ;
% 
% for idx=0:(k-1)
%    
%     plot(xOnDiag( (idx*6)+1:((idx+1)*6) ) ,...
%         yOnDiag(  (idx*6)+1:((idx+1)*6) ) ,...
%         'Color',cmap_mod(idx+1,:),'linewidth',4.5)
%     
%     
% end

set(gca,'ytick',[])
set(gca,'xtick',[])

axis square

hold off;

%% 

fileName = 'templateModel_adj_mat.png';
ff = fullfile(strcat('reports/figures/',fileName)); 
set(gcf,'paperpositionmode','auto');
print(gcf,'-dpng','-r300',ff);

%% 


