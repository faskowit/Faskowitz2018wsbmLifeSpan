function [B,E,K] = eval_genWsbm_model2(wsbmModel,D,numSims)
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
%           D x n,      Empirical measurement matricies that you want to 
%                       evaluate energy relative to. Examples could be 
%                       euclidean distance/fiber length matrix. Stack them.
%                       We will automatically test agains weighted degree
%                       distribution and assortativity too, therefore the
%                       number of distributions will be n + 2
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
%   Josh Yo

if nargin < 3
    numSims = 100 ;
end

Atgt = wsbmModel.Data.Raw_Data ;
%replace NaN's with 0
Atgt(~~isnan(Atgt)) = 0 ;
%remove diagonal to be safe
nNodes = size(Atgt,1) ;
Atgt(1:nNodes+1:end)=0; 

% get the matricies to test against
numMeas = size(D,3) + 2 ;

% emperical stats
x = cell(4,1);
x{1} = sum(Atgt,2);
x{2} = local_assortativity_wu_sign(Atgt);

for idx=3:(numMeas)
    
    tmpMat = D(:,:,idx) ;
    x{idx} = tmpMat(triu(Atgt,1) > 0);
    
end

%x{2} = clustering_coef_wu(Atgt);
%x{3} = betweenness_wei(Atgt);
%x{4} = D(triu(Atgt,1) > 0);

% record K-S stats
K = zeros(numSims,numMeas);

% record simulated networks
B = zeros([ nNodes nNodes numSims]);

for idx = 1:numSims
    
    % recover the no-NaN output
    [~,b] = genAdj_wsbm(wsbmModel) ;
    B(:,:,idx) = b ;
    
    y = cell(numMeas,1);
    y{1} = sum(b,2);
    y{2} = local_assortativity_wu_sign(b);
    
    for kdx=3:(numMeas)  
        tmpMat = D(:,:,kdx) ;
        y{kdx} = tmpMat(triu(b,1) > 0);  
    end
    
    %y{3} = betweenness_wei(b); % replaced with weighted
    %y{4} = D(triu(b,1) > 0);
    
    for jdx = 1:numMeas
        K(idx,jdx) = fcn_ks(x{jdx},y{jdx});
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

