function [ centralModel , allModels ] = wsbmCentralFit( adjMat, rStruct , modelInputs , numFits , priorMu)
% will fit wsbm specificed number iterations
% and will return 'most central' based on variation of information between
% community labels
% 
% user has the option to also provide prior labels to additionally multiply
% by when computing the central mode. this will be useful if you want to
% keep the central model 'near' a prior 

if nargin < 4
    disp('need more args')
    return
end

if nargin < 5
    disp('no prior read, this is cool')
    priorMu = '';
end

tempModelStruct = struct() ; 

%fit numFits times at each iter
for idx=1:numFits % model fits per iteration 

    [~,tempModel] = wsbm(adjMat, ...
        rStruct, ...
        modelInputs{:} ) ;

    % to save space
    tempModel = rmfield(tempModel,'Data') ;
    tempModelStruct(idx).Model = tempModel ;

end %iterating over number of model fits

centralModel = central_model( tempModelStruct , priorMu) ;
allModels = tempModelStruct ;
