function [sumBlockMat,avgBlockMat,binBlockMat] = get_block_mat(CIJ,ca)
% given an adjacency matrix, return a block matrix 

% make sure ca is column
ca = ca(:);

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
sumBlockMat = zeros([nBlocks nBlocks]);
avgBlockMat = zeros([nBlocks nBlocks]);
binBlockMat = zeros([nBlocks nBlocks]);

% now loop it
for idx = 1:nBlocks % rows
    for jdx = 1:nBlocks %columns

        tmp = CIJ(ca == orderedBlocks(idx),ca == orderedBlocks(jdx));
         
        sumBlockMat(idx,jdx) = nansum(tmp(:));
        avgBlockMat(idx,jdx) = nansum(tmp(:)) ./ sizesMat(idx,jdx);
        binBlockMat(idx,jdx) = nansum(tmp(:) > 0);
              
    end
end

















