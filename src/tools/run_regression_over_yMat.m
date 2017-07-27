function statMat = run_regression_over_yMat(xVec,yMat,fits,funcArgs)
% making the code more modular

if ~exist('funcArgs','var') || isempty(funcArgs)
    funcArgs = {1 500 [] 5000 } ;
end

% how many blocks we iterating over
numBlocks1 = size(yMat,1) ;
numBlocks2 = size(yMat,2) ;
numFits = length(fits) ;

% initialize cell mat, in which we will hold lots of info
statMat = cell([numBlocks1 numBlocks2 numFits]);

% regression work done with this function
% function [xvalR2, xvalsqErr, yhatLOOCV, coefStruct , lsFitStruct , permStruct ] = ... 
% nc_FitAndEvaluateModels(y, x, model, crossvalidate, bootIter, params , permIter)

% loop over blocks in upper triangle

% parallel?
parallel_pool = gcp ; 
parfor idx=1:numBlocks1
    for jdx=1:numBlocks2
        
        if jdx < idx
           continue 
        end
     
        % just for visual feedback
        disp(strcat(num2str(idx),'-',num2str(jdx)))
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the Y vector we want out of the statMap
        Y = squeeze(yMat(idx,jdx,:)) ;        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for fdx=1:numFits

            % initalize struct 
            %rdens_statMat{idx,jdx,fdx} = struct();
            statMatEntry = struct() ;
            [ statMatEntry.xvalR2 , ...
                statMatEntry.xvalsqErr, ... 
                statMatEntry.xvalYhat, ...
                statMatEntry.coef, ...
                statMatEntry.lsFitStruct,...
                statMatEntry.permStruct ] ...
                = nc_FitAndEvaluateModels(Y,xVec,fits{fdx},funcArgs{:}) ;
  
            % get meanSqErr
            statMatEntry.xvalRMSE = sqrt(mean(statMatEntry.xvalsqErr));
            % get medianSqErr
            statMatEntry.xvalRMedSE = sqrt(median(statMatEntry.xvalsqErr));
            % median absolute deviation    
            statMatEntry.xvalMADE = ...
                median( sqrt(statMatEntry.xvalsqErr) ); 

            statMat(idx,jdx,fdx) = {statMatEntry} ;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   end;
end;