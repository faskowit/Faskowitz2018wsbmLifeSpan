function gatherStruct = ma_gather_data(dataStruct , datasetDemo, templateModel, selectNodesFrmRaw, getModularity)
% put the data gathering steps into function to make the code more compact

%TODO MAKE THIS BETTER! MODULARIZEE

% gather the stuff into a struct
gatherStruct = struct() ; 

%% get labels

communityLabels_wsbm = community_assign(templateModel.Para.mu) ; 

if ~exist('getModularity','var') || isempty(getModularity) || getModularity ~= false
    [ communityLabels_mod , modResults ] = compute_mod(templateModel,1.0) ;
end

% this is now done in main preprocessing script
% %% condition the data based on exclusion criteria
% 
% nSubj = length(dataStruct) ;
% sparse_cutoff = 0.25 ;
% 
% removeVec = zeros( [nSubj 1] );
% 
% for idx=1:nSubj
%    
%     tmpSubjMat = ...
%     dataStruct(idx).countVolNormMat(selectNodesFrmRaw, ...
%     selectNodesFrmRaw ) ;
%     
%     sparseness = density_und(tmpSubjMat);
%     
%     if sparseness < sparse_cutoff
%        
%         removeVec(idx) = 1 ;
%     end
% end
% 
% % take out bad subjects
% % by filters by subjects not marked for removal
% dataStruct = dataStruct(removeVec == 0);
% datasetDemo = datasetDemo(removeVec == 0, :);

%% preallocate some vectors and friends

% make a matrix of all the raw data
% each entry of fitModelsStruct will have this information 
nSubj = length(dataStruct) ;
adjDataSize = size(dataStruct(1).countVolNormMat(selectNodesFrmRaw,selectNodesFrmRaw )) ;

% pre-allocate arrays because thats how I roll
subjMatsArray = zeros( [ adjDataSize nSubj] ) ;
% vectors
density = zeros( [ nSubj 1 ] ) ;
densityBin = zeros( [ nSubj 1 ] ) ;
sparseness = zeros( [ nSubj 1 ] ) ;
assortaCoef = zeros( [ nSubj 1] ) ;

eignVecCent = zeros( [ adjDataSize(1) nSubj] ) ;
particiCoef_wsbm = zeros( [ adjDataSize(1) nSubj] ) ;
particiCoef_mod = zeros( [ adjDataSize(1) nSubj] ) ;

% mats
edgeBtwnessCent = zeros( [ adjDataSize nSubj] ) ;

% covariate
% icv = datasetDemo.ICV ;

% partcipation coeff?
% other stuff?

%% gather the data, get some measures
% gather data here into matrix arrays
for idx=1:nSubj
        
    disp(num2str(idx))
    
    % raw data
    subjMatsArray(:,:,idx) = ...
        dataStruct(idx).countVolNormMat(selectNodesFrmRaw, ...
        selectNodesFrmRaw ) ;
    
    tmpSubjMat = squeeze(subjMatsArray(:,:,idx)) ;
     
    % get rid of the on-diagonal connections, as we didnt use these anyways
    % for any of the SBM stuff
    %make zero across diagonal
    n = adjDataSize;
    tmpSubjMat(1:n+1:end) = 0;
    
    % data derivatives
    % whole mat descriptors
    density(idx) = nansum(tmpSubjMat(:));
    densityBin(idx) = sum(tmpSubjMat(:)>0);
    sparseness(idx) = density_und(tmpSubjMat); 
    assortaCoef(idx) = assortativity_wei(tmpSubjMat,0) ;
    
    % node descriptions, nxn
    edgeBtwnessCent(:,:,idx) = edge_betweenness_wei(-log(tmpSubjMat));
    
    % node desciptions nx1
    eignVecCent(:,idx) = eigenvector_centrality_und(tmpSubjMat) ;
    particiCoef_wsbm(:,idx) = participation_coef(tmpSubjMat, ...
        communityLabels_wsbm(:,2)) ;
    particiCoef_mod(:,idx) = participation_coef(tmpSubjMat, ...
        communityLabels_mod(:,2)) ;
    
end

%% loop over the blocks

% next we will loop over blocks
numBlocks = templateModel.R_Struct.k ;

for idx=1:numBlocks
    for jdx=1:numBlocks
            
        ffi = find(communityLabels_wsbm(:,2)==idx);
        ffj = find(communityLabels_wsbm(:,2)==jdx);

        for kdx=1:nSubj
            
            % get adjMat
            tmpSubjMat = squeeze(subjMatsArray(ffi,ffj,kdx));
            dens(idx,jdx,kdx) = nansum(tmpSubjMat(:));
            bdens(idx,jdx,kdx) = sum(tmpSubjMat(:)>0);
            
            % get mat of edge-between
            tmpSubjMat = squeeze(edgeBtwnessCent(ffi,ffj,kdx));
            denseb(idx,jdx,kdx) = nansum(tmpSubjMat(:));

        end;
       
        % density with covariates
        [ ~ , ~ , res] = regress(squeeze(dens(idx,jdx,:)),[ones(nSubj,1) datasetDemo.SupraTentorialNotVent density]);
        rdens(idx,jdx,:) = res;
        
        % eb with covariates 
        [ ~ , ~ , res] = regress(squeeze(denseb(idx,jdx,:)),[ones(nSubj,1) datasetDemo.SupraTentorialNotVent density]);
        rdenseb(idx,jdx,:) = res;
        
    end;
end;

%% modularity gather

numNodes = size(templateModel.Para.mu,2) ;

% pick out the partician that is same # blocks
idx = max(modResults(2:end,:)) == templateModel.R_Struct.k ;
modResultsSubS = modResults(:,idx) ;
[ ~ , idx ] = max(modResults(1,idx)) ; 
communityLabels_mod_matchK = [ (1:numNodes)' modResultsSubS(2:end,idx) ] ;

for idx=1:numBlocks
    for jdx=1:numBlocks
            
        ffi = find(communityLabels_mod_matchK(:,2)==idx);
        ffj = find(communityLabels_mod_matchK(:,2)==jdx);

        for kdx=1:nSubj
            
            % get adjMat
            tmpSubjMat = squeeze(subjMatsArray(ffi,ffj,kdx));
            dens(idx,jdx,kdx) = nansum(tmpSubjMat(:));
            %bdens(idx,jdx,kdx) = sum(tmpSubjMat(:)>0);
            
            % get mat of edge-between
            tmpSubjMat = squeeze(edgeBtwnessCent(ffi,ffj,kdx));
            denseb(idx,jdx,kdx) = nansum(tmpSubjMat(:));

        end;
       
        [ ~ , ~ , res ] = regress(squeeze(dens(idx,jdx,:)),[ones(nSubj,1) datasetDemo.ICV density]);
        mod_rdens(idx,jdx,:) = res;
      
        [ ~ , ~ , res ] = regress(squeeze(denseb(idx,jdx,:)),[ones(nSubj,1) datasetDemo.ICV density]);
        mod_rdenseb(idx,jdx,:) = res;
        
    end;
end;

%% gather baby gather

% wsmb
gatherStruct.wsbm.dens = dens ;
gatherStruct.wsbm.rdens = rdens ;
gatherStruct.wsbm.denseb = denseb ; 
gatherStruct.wsbm.rdenseb = rdenseb ;

% mod
gatherStruct.mod.modResults = modResults ;
gatherStruct.mod.rdens = mod_rdens;
gatherStruct.mod.rdenseb = mod_rdenseb ;

% other derivatives
gatherStruct.deriv.densities = density ;
gatherStruct.deriv.densityBin = densityBin ;
gatherStruct.deriv.sparseness = sparseness ;
gatherStruct.deriv.eignVecCent = eignVecCent ;
gatherStruct.deriv.particiCoef_wsbm = particiCoef_wsbm ;
gatherStruct.deriv.particiCoef_mod = particiCoef_mod ;
gatherStruct.deriv.assortaCoef = assortaCoef ;

% conditioned data and demo
gatherStruct.datasetDemo = datasetDemo ;
gatherStruct.dataStruct = dataStruct ; 

