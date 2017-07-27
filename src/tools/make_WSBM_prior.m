function [ priorMat , harshPriorMat ] = make_WSBM_prior(inputMu,wieghtMore)
% takes in a WSBM model, and ouputs at prior matrix based of the mu_0
% of the WSBM model provided, and a harsh prior with just the 
% 'winners' of the mu_0 mat

% lets make the prior assignment matrix
% initialize the matrix
%priorMat = ones(wsbmModel.R_Struct.k , length(wsbmModel.Para.mu)) ;

if nargin < 2
    wieghtMore = 0 ; 
end

% if inputMu is a struct, extract relevant variables
if isstruct(inputMu)
    k = size(inputMu.Para.mu,1);
    mu = inputMu.Para.mu ;
else
    k = size(inputMu,1) ;
    mu = inputMu ; 
end

%twice as likely (1 would be 100%, double more likely)

weight_label = (1 + wieghtMore ) / (k + wieghtMore) ;
weights_others = 1 / ( k + wieghtMore ) ; 

% and then fill the mu_prior initially with the weights_others
priorMat = (ones( k , length(mu))) * weights_others ;

%sanity check 
if (weight_label < weights_others) 
    disp('prior weight is too low dude')
    return
end

%loop through the best_model.para.mu 
for i = 1:length(mu),
    [~,pos] = max(mu(:,i)) ; 
    priorMat(pos,i) = weight_label ;
end

%lets alos make the harsh prior
harshPriorMat = zeros( k , length(mu)) ;
for i = 1:length(mu),
    [~,pos] = max(mu(:,i)) ; 
    harshPriorMat(pos,i) = 1 ;
end

