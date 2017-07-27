function [] = plotMUnWSBM(Model)
% super simple

% plot
subplot(1,2,1) ; 
colormap('default')
plotWSBM(Model.Data.Raw_Data,Model.Para.mu)
subplot(1,2,2) ; 
plotMu(Model)

subplot(1,2,1) ; 


subplot(1,2,1) ; 
