function [ MSE_e , MSE_w ] = crossValWSBM( trainedModel, testMat )
% lets wrap the cross-val code part of the WSBM package

testEdges = Adj2Edg(testMat,1) ;

Error = @(x,y) (x-y).^2; 

MSE_e = predictE_Error(trainedModel,testEdges,Error);
MSE_w = predictW_Error(trainedModel,testEdges,Error);



