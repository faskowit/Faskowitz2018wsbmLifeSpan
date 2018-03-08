addpath('figures/DrosteEffect-BrewerMap-04533de/')
addpath('figures/')

sp1 = subplot(1,2,1) ;
imagesc(caMat_newAlign)
set(gca,'xtick',1:10)

%color
map = brewermap(10,'Spectral'); 
colormap(sp1,map)

title('Age bin template comm. assignments')
xlabel('Comunities at each bin')
ylabel('Nodes')

sp2 = subplot(1,2,2) ;
imagesc(agreementMat)
set(gca,'ytick',[])
% [left, bottom, width, height]
p2 = get(sp2,'pos') ;
p2(1) = p2(1) * 0.83 ;
set(sp2,'pos',p2)

%colorbar
map = brewermap(11,'GnBu'); 
map(1,:) = [1 1 1]; % optionally force first color to white 
colormap(sp2,map)
cb = colorbar(); 
% caxis([0 20]) % sets colorbar limits 
set(cb,'xtick',(0:10))
set(cb,'ticklength',0.01)
caxis([0 11])

xlabel('Nodes')
title('Agreement matrix across templates')




