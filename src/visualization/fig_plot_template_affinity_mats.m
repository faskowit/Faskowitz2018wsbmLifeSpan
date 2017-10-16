%% FIGURE
% visualize the parameter matricies of the template model

%%

predicted_w = make_square(templateModel.Para.predict_w) ;

subplot(1,2,1)
h = imagesc(predicted_w) ;
set(h,'alphadata',predicted_w > 0.0001);
%h.AlphaData = h.AlphaData * 0.25 ;

axis square

cb = colorbar;
%ylabel(cb,'Streamline Density')
cb.Label.String = 'Predicted Weight' ;
cb.Label.FontSize = 14 ;
cb.Label.FontName = 'Arial';

set(gca,'ytick',1:10)
set(gca,'xtick',1:10)

%% 

predicted_e = make_square(templateModel.Para.predict_e) ;

subplot(1,2,2)
h = imagesc(predicted_e) ;
set(h,'alphadata',predicted_e > 0.0001);
%h.AlphaData = h.AlphaData * 0.25 ;

axis square

cb = colorbar;
%ylabel(cb,'Streamline Density')
cb.Label.String = 'Predicted Edge' ;
cb.Label.FontSize = 14 ;
cb.Label.FontName = 'Arial';

set(gca,'ytick',1:10)
set(gca,'xtick',1:10)













