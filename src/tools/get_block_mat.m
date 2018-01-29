function [weiBM,avgWeiBM,binBM,avgBinBM,stdWeiBM] = get_block_mat(CIJ,ca)
% given an adjacency matrix + community affiliations, return a block matrix
% returns: 
%           weighted sum block matrix
%           average weighted block matrix
%           binary sum block matrix
%           average binary block matrix
%           std weights block matrix 
%
% Josh Faskowtiz IU

% make ca column vec
if ~iscolumn(ca)
   ca = ca'; 
end

% number coms
nBlocks = length(unique(ca));

% do this in case the user has provided ca with numeric gaps, eg ca=[1,2,4] 
orderedBlocks = sort(unique(ca));

% number nodes per block
blockSizes = histcounts(sort(ca));

sizesMat = bsxfun(@times,...
    (bsxfun(@times,ones(nBlocks),blockSizes)),...
    blockSizes');

% initialize outputs
weiBM = zeros([nBlocks nBlocks]);
avgWeiBM = zeros([nBlocks nBlocks]);
binBM = zeros([nBlocks nBlocks]);
avgBinBM = zeros([nBlocks nBlocks]);
stdWeiBM = zeros([nBlocks nBlocks]);

% now loop it
for idx = 1:nBlocks % rows
    for jdx = 1:nBlocks %columns

        tmp = CIJ(ca == orderedBlocks(idx),ca == orderedBlocks(jdx));
         
        weiBM(idx,jdx) = nansum(tmp(:));
        avgWeiBM(idx,jdx) = nansum(tmp(:)) ./ sizesMat(idx,jdx);
        binBM(idx,jdx) = nansum(tmp(:) > 0);
        avgBinBM(idx,jdx) = nansum(tmp(:) > 0) ./ sizesMat(idx,jdx);
        stdWeiBM(idx,jdx) = nanstd(tmp(:)) ;
         
    end
end
