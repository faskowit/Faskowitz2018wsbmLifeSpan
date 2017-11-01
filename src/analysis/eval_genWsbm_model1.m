function [B,E,K] = eval_genWsbm_model(wsbmModel,D,numSims)
% EVAL_GENWSBM_MODEL     Generation and evaluation of synthetic networks
%
%   [B,E,K] = EVALUATE_GENERATIVE_MODEL(A,Atgt,D,m,modeltype,modelvar,params) 
%
%   Generates synthetic networks and evaluates their energy function (see
%   below) using the models described in the study by Betzel et al (2016)
%   in Neuroimage.
%
%   Inputs:
%           wsbmModel,  WSBM model struct output from wsbm fitting
%           D,          Euclidean distance/fiber length matrix
%           numSims,    Number of simulated networks to generate (1000)
%
%   Outputs:
%           B,          n x n x numSims matrix of synthetic networks
%           E,          energy for each synthetic network
%           K,          Kolmogorov-Smirnov statistics for each synthetic
%                       network.
%
%
%   Note: Energy is calculated in exactly the same way as in Betzel et
%   al (2016). There are four components to the energy are KS statistics
%   comparing degree, clustering coefficient, betweenness centrality, and 
%   edge length distributions. Energy is calculated as the maximum across
%   all four statistics.
%
%   Reference: Betzel et al (2016) Neuroimage 124:1054-64.
%
%   Richard Betzel, Indiana University/University of Pennsylvania, 2015
%   Josh Faskowitz yo

if nargin < 3
    numSims = 100 ;
end

Atgt = wsbmModel.Data.Raw_Data ;
%replace NaN's with 0
Atgt(~~isnan(Atgt)) = 0 ;
%remove diagonal to be safe
nNodes = size(Atgt,1) ;
Atgt(1:nNodes+1:end)=0; 

% emperical stats
x = cell(4,1);
x{1} = sum(single(triu(Atgt,1) > 0),2);
x{2} = local_assortativity_wu_sign(Atgt);
%x{3} = eigenvector_centrality_und(Atgt);
x{3} = betweenness_bin(single(triu(Atgt,1) > 0))';
x{4} = D(triu(Atgt,1) > 0);

% record K-S stats
K = zeros(numSims,4);

% record simulated networks
B = zeros([ nNodes nNodes numSims]);

for idx = 1:numSims
    
    % recover the no-NaN output
    [~,b] = genAdj_wsbm(wsbmModel) ;
    B(:,:,idx) = b ;
    
    y = cell(4,1);
    y{1} = sum(b,2);
    y{2} = local_assortativity_wu_sign(b); 
    y{3} = betweenness_bin(single(triu(b,1) > 0))';
    %y{3} = eigenvector_centrality_und(b); 
    y{4} = D(triu(b,1) > 0);
    for j = 1:4
        K(idx,j) = fcn_ks(x{j},y{j});
    end
    
    disp(idx)
end
E = max(K,[],2);

function kstat = fcn_ks(x1,x2)
binEdges    =  [-inf ; sort([x1;x2]) ; inf];

binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
kstat = max(deltaCDF);

