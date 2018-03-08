function [ pl , ln ] = viz_blockRegress(statEntry,vertexBool,colorVec)

% take care of possible empty vars
if ~exist('vertexBool','var') || isempty(vertexBool)
    vertexBool = 0;
end

if ~exist('colorVec','var') || isempty(colorVec)
    colorVec = [];
end
    
xVec = min(statEntry.coef.x):0.1:...
    max(statEntry.coef.x);

fitString = statEntry.coef.name ;

% get the predicted yVector to draw
switch(fitString)
    case {'linear' 'quadratic'}
        
        yhat = polyval(statEntry.coef.full,xVec); 
   
        for idx = 1:size(statEntry.coef.boot,1)
            bootCI(idx,:) = polyval(statEntry.coef.boot(idx,:),xVec);  
            % and confidence interval for vertex
            vCI(idx) = -(statEntry.coef.boot(idx,2)./(2*statEntry.coef.boot(idx,1)));
        end
        
    case {'poisson'}

        yhat = evalPoissonCurve(statEntry.coef.full,xVec);
        
        for idx = 1:size(statEntry.coef.boot,1)
            bootCI(idx,:) = evalPoissonCurve(statEntry.coef.boot(idx,:),xVec);  
            % and confidence interval for vertex
            vCI(idx) = 1./statEntry.coef.boot(idx,2);
        end
end

% get the confidence intervals
p95 = prctile(bootCI,[2.5,97.5]);
v95 = prctile(vCI,[2.5,97.5]);

%figure stuffs
%figure;

%plot(statEntry.coef.x, statEntry.coef.y, 'o','color',[0 0 0],'markersize',6);
% pl = plot(statEntry.coef.x, statEntry.coef.y, 'o','markersize',6);
pl = scatter(statEntry.coef.x, statEntry.coef.y);
hold on;

%fill([xVec fliplr(xVec)],[p95(1,:) fliplr(p95(2,:))],[ 0.75 0.75 0.75 ],'facealpha',.5,'edgealpha',0);
fill([xVec fliplr(xVec)],[p95(1,:) fliplr(p95(2,:))],[ 0.75 0.75 0.75 ],'facealpha',.7,'edgealpha',0);

%plot(xVec,yhat, 'color',[ 0.1 0.1 0.1 ],'linewidth',2);
if isempty(colorVec)
    ln = plot(xVec,yhat,'linewidth',3);
else
    ln = plot(xVec,yhat,'linewidth',3,'Color',colorVec);
end

if vertexBool > 0
    
    if ~strcmp(fitString,'linear')
        % vertex pos
        % put it at the bottom essentially 
        ylims = ylim() ;
        
        yRange = ylims(2) - ylims(1);

        if isempty(colorVec)
            plot(v95,repmat(ylims(1) + (yRange .* 0.02) + vertexBool,1,2),'linewidth',7);
            % set ythe ylim back to before
            %ylim([ylims(1)*1.02 ylims(2)]);
        else
            plot(v95,repmat(ylims(1) + (yRange .* 0.02) + vertexBool,1,2),'linewidth',7,...
                'Color',colorVec);
        end
   end
end

hold off

end









