clc
clearvars

load('data/processed/yeo_both_normalpoisson_a0p5_fit_wsbm_script_v7p3.mat')
load('data/interim/yeo_both_normalpoisson_a0p5_comVecs.mat')

%% need to compute mod on everybody

fitModAllStruct = struct() ;

nSubj = length(fitWSBMAllStruct);

gammaRange = 0.5:0.001:4.0 ;

parallel_pool = gcp ; 
%for idx = 1:1
parfor idx=1:nSubj
    
    tmpAdj = fitWSBMAllStruct(idx).centModel.Data.Raw_Data ;
    tmpAdj(isnan(tmpAdj)) = 0 ;
   
    [fitModAllStruct(idx).sweepModCI] = compute_mod_sweepG(tmpAdj,gammaRange);
     
end

%% now get the best parition outta the sweep...

caMod_all_aligned = zeros([templateModel.Data.n nSubj]);

for idx=1:nSubj
    
    modSubjResults = fitModAllStruct(idx).sweepModCI ;
    
    matchKind = max(modSubjResults(2:end,:)) == templateModel.R_Struct.k ;
    
    if sum(matchKind) == 0
        
        fitModAllStruct(idx).caMod_subj = zeros([114 1]) ;
        fitModAllStruct(idx).caMod_gamma = 0 ;
        caMod_all_aligned(:,idx) = zeros([114 1]) ;
        
        continue
    end
    
    modResultsSubS = modSubjResults(:,matchKind) ;
    modResultsSubGammaVals = gammaRange(matchKind) ;
      
    % lets loop over the results to find the modular parition that is least
    % distant from the WSBM parition-- for a fair comparison. 
    modDist2Wsbm = zeros([ size(modResultsSubS,2) 1 ]);
    for jdx=1:size(modResultsSubS,2)

        modDist2Wsbm(jdx) = partition_distance(comVecs.wsbm, ...
            modResultsSubS(2:end,jdx));

    end

    [ ~ , minIdx ] = min(modDist2Wsbm) ; 

    % align labels
    subj_comVecMod = modResultsSubS(2:end,minIdx) ;
    subj_comVecMod = CBIG_HungarianClusterMatch(comVecs.wsbm,subj_comVecMod);
    subj_comVecModGamma = modResultsSubGammaVals(minIdx);
       
    fitModAllStruct(idx).caMod_subj = subj_comVecMod ;
    fitModAllStruct(idx).caMod_gamma = subj_comVecModGamma ;
    caMod_all_aligned(:,idx) = subj_comVecMod ;
    
end

%% save results 

outName = '/home/jfaskowi/JOSHSTUFF/projects/sbm3/data/processed/yeo_both_normalpoisson_a0p5_fit_mod_v7p3.mat' ;
save(outName,...
    'fitModAllStruct',...
    'caMod_all_aligned',...
    '-v7.3')









