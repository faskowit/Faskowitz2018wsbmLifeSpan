%% analysis script
% where I'll run some stuff

% add BCT path
addpath('~/JOSHSTUFF/scripts/BCT/2017_01_15_BCT/')

% read in the demographic files
%demoTable = datasetDemo ;

% make a matrix of all the raw data
% each entry of fitModelsStruct will have this information 
nSubj = length(dataStruct) ;
adjDataSize = size(dataStruct(1).countVolNormMat(selectNodesFrmRaw,selectNodesFrmRaw )) ;

% pre-allocate arrays because thats how I roll
subjMats = zeros( [ adjDataSize nSubj] ) ;
% vectors
%density = zeros( [ nSubj] ) ;
%densityBin = zeros( [  nSubj] ) ;
% mats
eb = zeros( [ adjDataSize nSubj] ) ;
ec = zeros( [ adjDataSize nSubj] ) ;
% partcipation coeff?
% other stuff?

%% gather the data, get some measures
% gather data here into matrix arrays
for idx=1:nSubj
        
    disp(num2str(idx))
    
    % raw data
    subjMats(:,:,idx) = ...
        dataStruct(idx).countVolNormMat(selectNodesFrmRaw, ...
        selectNodesFrmRaw ) ;
    
    subjMat = squeeze(subjMats(:,:,idx)) ;
    
    % data derivatives
    % whole mat descriptors
    density(idx) = nansum(subjMat(:));
    densityBin(idx) = sum(subjMat(:)>0);
    
    % node descriptors 
    eb(:,:,idx) = edge_betweenness_wei(-log(subjMat));
    
end

%%

% next we will loop over blocks
numBlocks = templateModel.R_Struct.k ;

for idx=1:numBlocks
    for jdx=1:numBlocks
            
        ffi = find(communityLabels(:,2)==idx);
        ffj = find(communityLabels(:,2)==jdx);
%        ffi = find(modular_community(:,2)==i);
%        ffj = find(modular_community(:,2)==j);

        for kdx=1:nSubj
            
            % get adjMat
            subjMat = squeeze(subjMats(ffi,ffj,kdx));
            dens(idx,jdx,kdx) = nansum(subjMat(:));
            bdens(idx,jdx,kdx) = sum(subjMat(:)>0);
            
            % get mat of edge-between
            subjMat = squeeze(eb(ffi,ffj,kdx));
            denseb(idx,jdx,kdx) = nansum(subjMat(:));

        end;
       
        [a b res] = regress(squeeze(dens(idx,jdx,:)),[density' ones(1,nSubj)']);
        rdens(idx,jdx,:) = res;
        
    end;
end;

%% make some R mats

age = datasetDemo.age;

for i=1:numBlocks
    for j=1:numBlocks
        
        LM = regstats(squeeze(dens(i,j,:)),age,'linear');
        Rl(i,j) = sqrt(LM.rsquare);
        LM = regstats(squeeze(dens(i,j,:)),age,'quadratic');
        Rq(i,j) = sqrt(LM.rsquare);
        
        LMb = regstats(squeeze(bdens(i,j,:)),age,'linear');
        Rlb(i,j) = sqrt(LMb.rsquare);
        LMb = regstats(squeeze(bdens(i,j,:)),age,'quadratic');
        Rqb(i,j) = sqrt(LMb.rsquare);
        
        LMr = regstats(squeeze(rdens(i,j,:)),age,'linear');
        Rlr(i,j) = sqrt(LMr.rsquare);
        LMr = regstats(squeeze(rdens(i,j,:)),age,'quadratic');
        Rqr(i,j) = sqrt(LMr.rsquare);
        
        LMe = regstats(squeeze(denseb(i,j,:)),age,'linear');
        Rle(i,j) = sqrt(LMe.rsquare);
        LMe = regstats(squeeze(denseb(i,j,:)),age,'quadratic');
        Rqe(i,j) = sqrt(LMe.rsquare);
   end;
end;






%% JUNK FROM BEFORE

uuu

%% lets display some reults on the brain 

%initialize mat
[temp1,temp2] = size(templateModel.Para.mu) ;
tempMuMats = zeros( temp1 , temp2 , length(fitWSBMAllStruct)) ; 

%lets quickly just plot binary mu
for i=1:length(fitWSBMAllStruct)
    
    disp(i)
    
    subjAdjMat = fitWSBMAllStruct(i).Model(1).Para.mu ;
    %threshold absolute
    subjAdjMat(subjAdjMat < 0.5) = 0 ; 
    %binarize
    subjAdjMat(subjAdjMat >= 0.5) = 1 ;
    size(subjAdjMat) ; 
    
    tempMuMats(:,:,i) = subjAdjMat ;
   
end

overall_muMean = mean(tempMuMats,3) ;

%% load the seven parcellation

load('aux_stuff/seven_network.mat')
load('aux_stuff/seventeen_network.mat')

% lets condition the networks if we only 
% looking at one hemisphere
if strcmp(HEMI_CHOICE,'left')
    seven_network = seven_network(:,1:57) ;
elseif strcmp(HEMI_CHOICE,'right')
    seven_network = seven_network(:,58:114) ;
elseif strcmp(HEMI_CHOICE,'both')
    disp('good to go')
end
    
%% make the results table
% lets get some stats

[numBlocks , ~ ] = size(templateModel.Para.mu) ;
logEvidences = zeros(length(fitWSBMAllStruct),1) ;

%adding more stuff
densities = zeros(length(fitWSBMAllStruct),1) ; 
% yeo7Dist
yeo7Dist = zeros(length(fitWSBMAllStruct),1) ;
% sim to prior
priorSim = zeros(length(fitWSBMAllStruct),1) ;
% dist to prior
priorDist = zeros(length(fitWSBMAllStruct),1) ;

% save a vector of the block existences
% predictE_Vecs = zeros(length(fitWSBMAllStruct), ...
%     ((numBlocks+1)*numBlocks/2 )) ;
predictE_Vecs = zeros(length(fitWSBMAllStruct), ...
    length(templateModel.Para.predict_e)) ;

modelIterIndx = 5 ;

for idx=1:length(fitWSBMAllStruct)
    
    logEvidences(idx) = fitWSBMAllStruct(idx).Model(modelIterIndx).Para.LogEvidence ;
    
    densities(idx) = density_und(fitWSBMAllStruct(idx).Raw_Data > 0) ;
        
    [ ~ , mu_temp ] = make_WSBM_prior(fitWSBMAllStruct(idx).Model(modelIterIndx)) ;
    
    yeo7Dist(idx) = varInfo(mu_temp,seven_network) ;
    priorSim(idx) = nmi(mu_temp , harsh_mu) ; 
    priorDist(idx) = varInfo(mu_temp, harsh_mu ) ;
    
    % get the predictE_Vecs
%     predictE_Vecs(idx,:) = nonzeros(triu(reshape(...
%         fitWSBMAllStruct(idx).Model(modelIterIndx).Para.theta_e,...
%         [numBlocks numBlocks])))' ; 
    predictE_Vecs(idx,:) = fitWSBMAllStruct(idx).Model(modelIterIndx).Para.predict_e ;
    
    % to recover the triu use the following:
    % b = triu(ones(numBlocks),0)
    % b(~~b) = predictE_Vecs(idx,:)'
    
end

%put it into a simple table
col_names = {'ids' 'log_evid'}; 
results_table = table({fitWSBMAllStruct.id}',logEvidences) ;%'RowNames',{fit_Models_Stuct.id}') ;  
results_table.Properties.VariableNames = col_names ;
results_table.Properties.RowNames = [fitWSBMAllStruct.id]' ;

%adding more stuff
results_table.densities = densities ; 
results_table.priorSim = priorSim ;
results_table.priorDist = priorDist ; 
results_table.yeo7Dist = yeo7Dist ;
results_table.predictE_Vecs = predictE_Vecs ;

results_table = join( results_table , datasetDemo , 'Keys' , 'RowNames') ; 

% remove density outliers
% watch out...
% will mess up corresponence to dataStruct
% removeVec = results_table.densities > 0.2 ;
% results_table = results_table(removeVec,:) ;

%% more cleaning to results table

results_table.sex = cellstr(results_table.sex) ;

%% lets do some binning

[ iterate , ~ ] = size(results_table) ;

ageBin = zeros(iterate,1) ;

numBins = 5 ;
[~,edge] = histcounts(results_table.age, numBins) ;
ageBin = discretize(results_table.age,edge) ;

results_table.ageBin1 = ageBin ;

% other type of age bin...

ageBin = zeros(iterate,1) ;

thresholds = 4 ;

low_quantile = [ 0 quantile(results_table.age,thresholds) ] ; 
high_quantile = [ quantile(results_table.age,thresholds) 100 ] ;

for idx=1:iterate
       
    for jdx=1:(thresholds+1)
    
       if ( results_table.age(idx) > low_quantile(jdx) ) && ( results_table.age(idx) <= high_quantile(jdx) )
           
           disp(idx)
           
           ageBin(idx) = jdx ;
           
       end     
        
    end    
end    

results_table.ageBin2 = ageBin ;

%% make some covariance at each age

numSubsets = 5 ;

predictEVec_crossAge = zeros( [ numBlocks numBlocks numSubsets ] ) ;

for subset=1:numSubsets

    pickoutVec = results_table.ageBin1 == subset ;

    subsetTable = results_table(pickoutVec,:) ;

    [ iterate , ~ ] = size(subsetTable) ;

    a = zeros( [ numBlocks numBlocks iterate ] ) ;

    for idx=1:1

        b = triu(ones(numBlocks),0) ;
        b(~~b) = subsetTable.predictE_Vecs(idx,:)' ;

        a(:,:,idx) = b ;

    end

    predictEVec_crossAge(:,:,subset) = mean(a,3) ;

end

%% overall muMean 

%make the overall mu. baed on overall_muMean
overall_mu = zeros(size(overall_muMean)) ; 
for i=1:length(overall_muMean)
    
    overall_mu(:,i) = overall_muMean(:,i) == max(overall_muMean(:,i)) ;  
 
end

%% plotting

subplot(1,3,1)
plotMu(overall_mu)
subplot(1,3,2)
plotMu(overall_muMean)
subplot(1,3,3)
plotMu(muPrior)

%% get comunity for each of the verticies 

community_assign2 = zeros(length(overall_mu),2) ; 

for i = 1:length(overall_mu),
    
    %make first column the label number
    community_assign2(i,1) = selectNodesFrmRaw(i) ; 
    
    %make second column the 
    [~,temp] = max(overall_mu(:,i)) ; 
    community_assign2(i,2) = temp ; 
    
end
    
disp('the block assisgnments of Init Model are:\n')
disp(community_assign2)

%remove temp vars
clear temp

%% what if from the looper, we read something like best 3 models...
% and then used x-val to pick the one with the lower error

% also can make cross val parition 
xValInd = crossvalind('Kfold',length(templateSubj_data),10) ;

trainData = templateSubj_data(xValInd ~= 10) ; 
testData = templateSubj_data(xValInd == 10) ;

% first lets get this working to test it out
trainTemplate = make_template_mat(trainData, ...
    LEFT_HEMI_NODES, ...
    RIGHT_HEMI_NODES, ...
    MASK_THR_INIT) ; 
trainTemplate(trainTemplate == 0) = NaN ;

testTemplate = make_template_mat(testData, ...
    LEFT_HEMI_NODES, ...
    RIGHT_HEMI_NODES, ...
    MASK_THR_INIT) ; 
testTemplate(testTemplate == 0) = NaN ;

% fit model to training data
% k=5 for test
[~,trainModel] = wsbm(trainTemplate, ...
    8, ... %diagWBlocks_RStruct(kBest), ... %sym_RStruct(kBest), ...
    'W_Distr', WEIGHT_DIST, ...
    'E_Distr',EDGE_DIST, ...
    'alpha', 0.5, ...
    'numTrials', 100 , ...
    'mainMaxIter', MODEL_WSBM_MAIN_ITER , ...
    'muMaxIter', MODEL_WSBM_MU_ITER );

[ MSEe , MSEw ] = crossValWSBM(trainModel,testTemplate) ;




