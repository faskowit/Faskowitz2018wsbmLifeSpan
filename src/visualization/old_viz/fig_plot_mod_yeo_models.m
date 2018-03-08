% lets make a nice picture of the template model

figure;

% algin just incase
modular_comVec = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.mod);

% get community assignment stuff
%[~,comVecs] = community_assign(comVecs.mod);
[xOnDiag,yOnDiag,sortIdx] = grid_communities(modular_comVec);

% code from WSBM stuffs
n = length(comVecs.mod);
k = length(unique(modular_comVec));
A_sort = zeros(n);
list = zeros(1,n);
breaks = zeros(1,k);
cur = 1;
for idx = 1:k
    indicies = find(modular_comVec == idx);
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

%plot on diagonal 
pp = plot(xOnDiag,yOnDiag,'Color',[1 0 0 1],'linewidth',2.5);

set(gca,'ytick',[])
set(gca,'xtick',[])

axis square
hold off;

%%
figure;

% algin just incase
yeo_comVec = CBIG_HungarianClusterMatch(comVecs.wsbm,comVecs.yeo);

% get community assignment stuff
%[~,comVecs] = community_assign(comVecs.mod);
[xOnDiag,yOnDiag,sortIdx] = grid_communities(yeo_comVec);

% code from WSBM stuffs
n = length(comVecs.mod);
k = length(unique(yeo_comVec));
A_sort = zeros(n);
list = zeros(1,n);
breaks = zeros(1,k);
cur = 1;
for idx = unique(yeo_comVec)' %1:k
    indicies = find(yeo_comVec == idx);
    list(cur:cur+length(indicies)-1) = indicies;
    cur = cur + length(indicies);
    breaks(idx) = cur-1;
end
for idx = 1:n
    A_sort(idx,:) = data(list(idx),list);
end

breaks(2) = 30 ;
breaks(6) = 80 ;
breaks(9) = 114 ;

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

%plot on diagonal 
pp = plot(xOnDiag,yOnDiag,'Color',[1 0 0 1],'linewidth',2.5);

set(gca,'ytick',[])
set(gca,'xtick',[])

axis square
hold off;

%%

pictureWSBM(yeo_comVec,'both')


%% 

% fileName = 'templateModel_adj_mat.png';
% ff = fullfile(strcat('reports/figures/',fileName)); 
% set(gcf,'paperpositionmode','auto');
% print(gcf,'-dpng','-r300',ff);

%% 


