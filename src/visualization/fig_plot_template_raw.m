%% lets make a nice picture of raw data
figure

load('data/external/seven_network.mat');
load('data/external/seventeen_network.mat');

data = templateModel.Data.Raw_Data;

ca_seven = sum(bsxfun(@times,seven_network,[1:(size(seven_network,1))]'))';
ca_seven(ca_seven == 0) = 8 ;

ca_seventeen = sum(bsxfun(@times,seventeen_network,[1:(size(seventeen_network,1))]'))';

% get community assignment stuff
[xOnDiag,yOnDiag,sortIdx] = grid_communities(ca_seven);

% % code from WSBM stuffs
% n = length(ca_seven);
% %k = length(unique(ca_seven));
% k = max(ca_seven);
% A_sort = zeros(n);
% list = zeros(1,n);
% breaks = zeros(1,k);
% cur = 1;
% for idx = (unique(ca_seven))'
%     indicies = find(ca_seven == idx); 
%     list(cur:cur+length(indicies)-1) = indicies;
%     cur = cur + length(indicies);
%     breaks(idx) = cur-1;
% end
% 
% %breaks(breaks == 0) = 114 ;
% breaks(6) = 114;
% %breaks(9) = -100;
% 
% for idx = 1:n
%     A_sort(idx,:) = data(list(idx),list);
% end

%Plot the Matrix
%h = imagesc(A_sort,[min(A_sort(:))-.00001,max(A_sort(:))]);
h = imagesc(data);
set(h,'alphadata',~isnan(data));
%h.AlphaData = h.AlphaData * 0.25 ;

cb = colorbar;
%ylabel(cb,'Streamline Density')
cb.Label.String = 'Streamline Density' ;
cb.Label.FontSize = 14 ;
cb.Label.FontName = 'Arial';

hold all; 

% % plot off diagonal
% for idx = 1:(max(ca_seven)-1)
%     
%     lineWidth = 1;
%     %onDiagColor = [1 0 0.75 0.25] ;
%     offDiagColor = [1 0 0 0.35] ; 
%     
%     
%     % vertical   
%     plot([breaks(idx)+0.5,breaks(idx)+0.5],[breaks(idx+1)+0.5,breaks(k)+.5],...
%         'Color',offDiagColor,'LineWidth',lineWidth);
%     if idx > 1
%         plot([breaks(idx)+0.5,breaks(idx)+0.5],[-0.5,breaks(idx-1)+0.5],...
%             'Color',offDiagColor,'LineWidth',lineWidth);
%     end
%         
%     plot([breaks(idx+1)+0.5,breaks(k)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
%         'Color',offDiagColor,'LineWidth',lineWidth);
%     if idx > 1
%         plot([-.5,breaks(idx-1)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
%             'Color',offDiagColor,'LineWidth',lineWidth);
%     end
%     
% end      
% 
% %plot on diagonal 
pp = plot(xOnDiag,yOnDiag,'Color',[1 0 0 1],'linewidth',1.5);
uistack(pp,'top');

% lets try to plot on diagonal with lines according to community
% assignments

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

% %% %Plot the Matrix
% h = imagesc(A_sort,[min(A_sort(:))-.00001,max(A_sort(:))]);
% set(h,'alphadata',~isnan(A_sort));
% h.AlphaData = h.AlphaData * 0.11 ;
% 
% %cb = colorbar;
% %ylabel(cb,'Streamline Density')
% %cb.Label.String = 'Streamline Density' ;
% %cb.Label.FontSize = 16 ;
% %cb.Label.FontName = 'Arial';
% 
% hold all; 
% 
% % plot off diagonal
% for idx = 1:(k-1)
%     
%     lineWidth = 1;
%     %onDiagColor = [1 0 0.75 0.25] ;
%     offDiagColor = [1 0 0 0.15] ; 
%     
%     
%     % vertical   
%     plot([breaks(idx)+0.5,breaks(idx)+0.5],[breaks(idx+1)+0.5,breaks(k)+.5],...
%         'Color',offDiagColor,'LineWidth',lineWidth);
%     if idx > 1
%         plot([breaks(idx)+0.5,breaks(idx)+0.5],[-0.5,breaks(idx-1)+0.5],...
%             'Color',offDiagColor,'LineWidth',lineWidth);
%     end
% 
%     
%     % horizontal
%     plot([breaks(idx+1)+0.5,breaks(k)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
%         'Color',offDiagColor,'LineWidth',lineWidth);
%     if idx > 1
%         plot([-.5,breaks(idx-1)+.5],[breaks(idx)+.5,breaks(idx)+.5],...
%             'Color',offDiagColor,'LineWidth',lineWidth);
%     end
%     
% end      
% 
% %plot on diagonal 
% %plot(xOnDiag,yOnDiag,'r','linewidth',2.5)
% 
% % lets try to plot on diagonal with lines according to community
% % assignments
% 
% cmap_mod = brewermap(k,'paired') ;
% 
% for idx=0:(k-1)
%    
%     plot(xOnDiag( (idx*6)+1:((idx+1)*6) ) ,...
%         yOnDiag(  (idx*6)+1:((idx+1)*6) ) ,...
%         'Color',cmap_mod(idx+1,:),'linewidth',3)
%     
%     
% end
% 
% set(gca,'ytick',[])
% set(gca,'xtick',[])
% 
% axis square
% 
% hold off;



%%

% fileName = 'templateModel_community_mat.png';
% ff = fullfile(strcat('reports/figures/',fileName)); 
% set(gcf,'paperpositionmode','auto');
% print(gcf,'-dpng','-r300',ff);


