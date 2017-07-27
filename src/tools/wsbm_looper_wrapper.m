function [ loopResults , allFitModels ] = wsbm_looper_wrapper( inputData , modelInputs , loopIters , scoreFunc )
% wrapper for the looper script already in the wsbm code

% initialize the variables we will write into
% so it parallels easily
numModels = numel(modelInputs) ;
loopResults = zeros(numModels, loopIters + 1) ; 
loopResults(1:numModels,1) = 1:numModels ;
%loopResults(numModels + 1, 1) = 999 ;

allFitModels = cell(numModels,loopIters) ;

if nargin < 4
  scoreFunc = @(model) model.Para.LogEvidence ;
  disp('using log evidence as score func')
else
  disp('using provided score func')
end

parallel_pool = gcp ; 
parfor idx = 1:loopIters

    disp('iteration:')
    disp(idx)
    
    % Fit
    [~, tempSores, tempModels] = wsbmLooper(inputData, ...
        modelInputs, ...
        scoreFunc);
    
    allFitModels(:,idx) = tempModels(:) ;
    
    %save the results 
    loopResults(:,idx + 1) = tempSores ;
    
end
