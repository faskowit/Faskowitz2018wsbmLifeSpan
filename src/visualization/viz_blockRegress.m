function [] = viz_blockRegress(statEntry,vertexBool)

% take care of possible empty vars
if ~exist('vertexBool','var') || isempty(vertexBool)
    vertexBool = 0;
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
figure;

%plot(statEntry.coef.x, statEntry.coef.y, 'o','color',[0 0 0],'markersize',6);
plot(statEntry.coef.x, statEntry.coef.y, 'o','markersize',6);
hold on;

fill([xVec fliplr(xVec)],[p95(1,:) fliplr(p95(2,:))],[ 0.75 0.75 0.75 ],'facealpha',.5,'edgealpha',0);

plot(xVec,yhat, 'color',[ 0.1 0.1 0.1 ],'linewidth',2);
%plot(xVec,yhat,'linewidth',2);

if vertexBool
    
    if ~strcmp(fitString,'linear')
        % vertex pos
        % put it at the bottom essentially 
        ylims = ylim() ;
        plot(v95,repmat(ylims(1) .* 1.05,1,2),'linewidth',3,'color',[ 0.9 0.9 0.9 ]);
        % set ythe ylim back to befor
        ylim([ylims(1)*1.12 ylims(2)]);
    end
end
end
