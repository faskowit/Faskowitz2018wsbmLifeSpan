function [ templateMatBothHemi , meanLens ] = make_template_mat(dataStruct, leftNodes, rightNodes, initThresh)
%% make average mat, and return it
% we will use function that preserves mat lengths
% therefore this function needs the additional function: fcn_group_average

selectNodes = [ leftNodes rightNodes ] ;

%% make average distCoorMM
% will compute both for fun but will only use the lengths of the fibers
% for when we run it, but we can change this

% euclidean distance

% arrayCoorMM = zeros([size(dataStruct(1).distCoorMM) length(dataStruct)]) ;
% for i=1:length(dataStruct)
%     arrayCoorMM(:,:,i) = dataStruct(i).distCoorMM ; 
% end
% meanCoorMM = mean(arrayCoorMM,3) ;
% 
% clear arrayCoorMM

% tract distance on average

arrayLens = zeros([size(dataStruct(1).lensMat) length(dataStruct)]) ;
for i=1:length(dataStruct)
    arrayLens(:,:,i) = dataStruct(i).lensMat ;
end
meanLens = mean(arrayLens,3) ;

%% make average mat
% by looping over data given 

% pre-allocate the matrix to store subject-level mats
sizeAllNodes = length(selectNodes) ;
arrayWeightMats = zeros(sizeAllNodes, sizeAllNodes, length(dataStruct)) ;

for idx=1:length(dataStruct)
           
    temp = dataStruct(idx).countVolNormMat(selectNodes,selectNodes) ;
    
    %get size of square mat
    n=size(temp,1);
    %make zero across diagonal
    temp(1:n+1:end) = 0;
    
    %make a mask based on thresh
    temp_mask = dataStruct(idx).countMat(selectNodes,selectNodes) > initThresh ;    
    temp_mask(temp_mask > 0) = 1 ;
    
    % threshold temp matrix by mask
    temp = temp .* temp_mask ;
    arrayWeightMats(:,:,idx) = temp ;
    
end

% make the lh rh membership mat
hemi_id = ones(sizeAllNodes,1);
temp = length(leftNodes) + 1 ; 
hemi_id(temp:end) = 2 ;

% using mean lengths here
template_mask = fcn_group_average(arrayWeightMats, ...
    meanLens(selectNodes,selectNodes), ...
    hemi_id) ; 

numEdg = sum(arrayWeightMats>0, 3);
avrWei = sum(arrayWeightMats, 3)./ numEdg;
avrWei(numEdg == 0) = 0;  % remove NANs

%mask the avrWei mat by the threshold mat
templateMatBothHemi = avrWei .* template_mask ;



