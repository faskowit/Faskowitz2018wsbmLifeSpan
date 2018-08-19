function [ loopResults , allFitModels ] = wsbm_looper_wrapper( inputData , modelInputs , loopIters , scoreFunc, numTrialPttrn, priorWeightPttrn)
% wrapper for the looper script already in the wsbm code

% initialize the variables we will write into
% so it parallels easily
numModels = numel(modelInputs) ;
loopResults = zeros(numModels, loopIters + 1) ; 
loopResults(1:numModels,1) = 1:numModels ;
%loopResults(numModels + 1, 1) = 999 ;

allFitModels = cell(numModels,loopIters) ;

if ~exist('scoreFunc','var') || isempty(scoreFunc)
  scoreFunc = @(model) model.Para.LogEvidence ;
  disp('using log evidence as score func')
else
  disp('using provided score func')
end

if ~exist('numTrialPttrn','var') || isempty(numTrialPttrn)
    numTrialPttrn = [] ;
end

if ~exist('priorWeightPttrn','var') || isempty(priorWeightPttrn)
    priorWeightPttrn = [] ;
end

parallel_pool = gcp ; 

%TMP
ppm1 = ParforProgMon('looper',loopIters,1) ;

parfor idx = 1:loopIters

    %disp('iteration:')
    %disp(idx)
    
    % Fit
    % function [Best_Models,Scores,Models] = wsbmLooper_2(E,ModelInputs,scorefuncs,numTrialPttrn,priorWeightPttrn)
    [~, tempSores, tempModels] = wsbmLooper_2(inputData, ...
        modelInputs, ...
        scoreFunc,...
        numTrialPttrn,...
        priorWeightPttrn);
    
    allFitModels(:,idx) = tempModels(:) ;
    
    %save the results 
    loopResults(:,idx + 1) = tempSores ;
    
    %TMP
    ppm1.increment();
    
end
