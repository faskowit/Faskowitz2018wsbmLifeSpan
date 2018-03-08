function gatherStruct = ma_gather_data2(dataStruct , datasetDemo, communityLabels, selectNodesFrmRaw, derivBool, thr, matDataStr)
% put the data gathering steps into function to make the code more compact

% gather the stuff into a struct
gatherStruct = struct() ; 

if ~exist('derivBool','var') || isempty(derivBool)
    derivBool = true ;
end

if ~exist('thr','var') || isempty(thr)
    thr = 0 ;
end
   
if ~exist('matDataStr','var') || isempty(matDataStr)
    matDataStr = 'countVolNormMat' ;
end

%% preallocate some vectors and friends

% make a matrix of all the raw data
% each entry of fitModelsStruct will have this information 
nSubj = length(dataStruct) ;
adjDataSize = size(dataStruct(1).(matDataStr)(selectNodesFrmRaw,selectNodesFrmRaw )) ;

numBlocks = length(unique(communityLabels)) ;

% pre-allocate arrays because thats how I roll
subjMatsArray = zeros( [ adjDataSize nSubj] ) ;

% vectors
totDensity = zeros( [ nSubj 1 ] ) ;
totDensityBin = zeros( [ nSubj 1 ] ) ;
sparseness = zeros( [ nSubj 1 ] ) ;

assortaCoef = zeros( [ adjDataSize(1) nSubj] ) ;
eignVecCent = zeros( [ adjDataSize(1) nSubj] ) ;

partiCoeff = zeros( [ adjDataSize(1) nSubj] ) ;
degreeZscor = zeros( [ adjDataSize(1) nSubj] ) ;

% mats
edgeBtwnessCent = zeros( [ adjDataSize nSubj] ) ;

% other network measures?

%% gather the data, get some measures
% gather data here into matrix arrays
for idx=1:nSubj
        
    disp(num2str(idx))
    
    % raw data
    tmpSubjMat = ...
        dataStruct(idx).(matDataStr)(selectNodesFrmRaw, ...
        selectNodesFrmRaw ) ;
     
    % get rid of the on-diagonal connections, as we didnt use these anyways
    % for any of the SBM stuff
    %make zero across diagonal
    n = adjDataSize;
    tmpSubjMat(1:n+1:end) = 0;
    
    % also threshold out by threshold we used initially when fitting the
    % wsbm... should be 1 or 0 (aka no thresh)
    temp_mask = ...
        dataStruct(idx).countMat(selectNodesFrmRaw,...
        selectNodesFrmRaw) > thr ;    
    temp_mask(temp_mask > 0) = 1 ;
    tmpSubjMat = tmpSubjMat .* temp_mask ;
    
    % data derivatives
    % whole mat descriptors
    totDensity(idx) = nansum(tmpSubjMat(:));

    partiCoeff(:,idx) = participation_coef(tmpSubjMat, ...
        communityLabels) ;
    degreeZscor(:,idx) = module_degree_zscore(tmpSubjMat, ...
        communityLabels, 0) ;
    
    if derivBool == true
        totDensityBin(idx) = sum(tmpSubjMat(:)>0);
        sparseness(idx) = density_und(tmpSubjMat);
        assortaCoef(idx) = assortativity_wei(tmpSubjMat,0) ;

        % node desciptions nx1
        eignVecCent(:,idx) = eigenvector_centrality_und(tmpSubjMat) ;

        edgeBtwnessCent(:,:,idx) = edge_betweenness_wei(-log(tmpSubjMat));
    end
    
    % subj data it in a mat array
    subjMatsArray(:,:,idx) = tmpSubjMat ;
    
end

%% loop over the blocks

dens = zeros([numBlocks numBlocks nSubj]) ;
bdens = zeros([numBlocks numBlocks nSubj]) ;
denseb = zeros([numBlocks numBlocks nSubj]) ;
rdens = zeros([numBlocks numBlocks nSubj]) ;
rbdens = zeros([numBlocks numBlocks nSubj]) ;
rdenseb = zeros([numBlocks numBlocks nSubj]) ;
rmdens = zeros([numBlocks numBlocks nSubj]) ;
rmbdens = zeros([numBlocks numBlocks nSubj]) ;
rmdenseb = zeros([numBlocks numBlocks nSubj]) ;

for idx=1:numBlocks
    for jdx=1:numBlocks
            
        ffi = find(communityLabels==idx);
        ffj = find(communityLabels==jdx);

        for kdx=1:nSubj
            
            % get adjMat
            tmpSubjMat = squeeze(subjMatsArray(ffi,ffj,kdx));
            dens(idx,jdx,kdx) = nansum(tmpSubjMat(:));
            bdens(idx,jdx,kdx) = sum(tmpSubjMat(:)>0);
            
            if derivBool == true
                % get mat of edge-between
                tmpSubjMat = squeeze(edgeBtwnessCent(ffi,ffj,kdx));
                denseb(idx,jdx,kdx) = nansum(tmpSubjMat(:));
            end
            
        end;

        % density with covariates
        [ ~ , ~ , res] = regress(squeeze(dens(idx,jdx,:)),[ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity]);
        rdens(idx,jdx,:) = res;
        
        % binary density with covariates
        [ ~ , ~ , res] = regress(squeeze(bdens(idx,jdx,:)),[ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity]);
        rbdens(idx,jdx,:) = res;
        
        if derivBool == true
            % eb with covariates 
            [ ~ , ~ , res] = regress(squeeze(denseb(idx,jdx,:)),[ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity]);
            rdenseb(idx,jdx,:) = res;
        end
        
        % density with covariates + movement
        [ ~ , ~ , res] = regress(squeeze(dens(idx,jdx,:)), ...
            [ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity datasetDemo.movement]);
        rmdens(idx,jdx,:) = res;
        
        % binary density with covariates + movement
        [ ~ , ~ , res] = regress(squeeze(bdens(idx,jdx,:)), ...
            [ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity datasetDemo.movement]);
        rmbdens(idx,jdx,:) = res;
        
        if derivBool == true
            % eb with covariates + movement
            [ ~ , ~ , res] = regress(squeeze(denseb(idx,jdx,:)),...
                [ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity datasetDemo.movement]);
            rmdenseb(idx,jdx,:) = res;
        end
        
    end;
end;

%% gather baby gather

% wsmb
gatherStruct.dens = dens ;
gatherStruct.rdens = rdens ;
gatherStruct.rbdens = rbdens ;

gatherStruct.rmdens = rmdens ;
gatherStruct.rmbdens = rmbdens ;

% subjMats in a mat array
gatherStruct.subjMatsArray = subjMatsArray ;

gatherStruct.deriv.particiCoef = partiCoeff ;
gatherStruct.deriv.degreeZscor = degreeZscor ;

if derivBool == true
    
    gatherStruct.denseb = denseb ; 
    gatherStruct.rdenseb = rdenseb ;
    gatherStruct.rmdenseb = rmdenseb ;

    % other derivatives
    gatherStruct.deriv.densities = totDensity ;
    gatherStruct.deriv.densityBin = totDensityBin ;
    gatherStruct.deriv.sparseness = sparseness ;
    gatherStruct.deriv.eignVecCent = eignVecCent ;
end

