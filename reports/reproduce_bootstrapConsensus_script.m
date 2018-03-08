%% load the data
clc 
clearvars

load('/home/jfaskowi/JOSHSTUFF/projects/sbm3/reports/reproduce_bootstrapConsensus_script_data.mat')

%% test consensus strategy
% lets do a test where we randomly pick x proportion of LOOPER_ITER and run
% the consensus method, and then pick consensus based on that proportion
% kiter initialization. this way we can test how robust our consensus
% method is.

%init a cell array to hold all boot iters results...
bt_kiter_prior_results = cell([50 1]) ;
bt_consensus_model = cell([50 1]) ;
bt_C = cell([50 1]) ;

%% run it

for idx=1:100
   
    btIdx(idx,:) = datasample(1:100,100) ;
    
end

parallel_pool = gcp ; 
ppm1 = ParforProgMon('indivFits',100,1) ;

for btiter=1:100

    %btIdx = datasample(1:100,100) ;

    btkBestModels_ca_algn_sort = kBestModels_ca_algn_sort(:,btIdx(btiter,:)) ;

    % iterate over each node
    for idx=1:(kCentralModel.Data.n)

        btkiter_prior(:,idx,1) = ...
            sum(bsxfun(@eq,btkBestModels_ca_algn_sort(idx,1:cutoff), ...
            [1:(kCentralModel.R_Struct.k)]'),2) ./ cutoff ;

    end  

    CONSENSUS_ITER = 10 ;
    btC = zeros([CONSENSUS_ITER 1]) ;

    for idx=1:CONSENSUS_ITER

        runNTimes = 100 ;

        rr = sym_RStruct(kBest) ;
        modIn = { ... 
            'W_Distr', WEIGHT_DIST, ...
            'E_Distr', EDGE_DIST, ...
            'alpha', INIT_ALPHA, ...
            'mainMaxIter', LOOPER_MAIN_ITER , ...
            'muMaxIter' , LOOPER_MU_ITER,...
            'mainTol',0.01, ...
            'muTol', 0.01 ,...
            'verbosity', 0, ...
            'numTrials', 50 ,...
            'mu_0', btkiter_prior(:,:,idx)};

        % function [ allModels ] = wsbmFitNTimes( adjMat, rStruct , modelInputs , numFits , numCores)
        btcnsnsusModels = wsbmFitNTimes(kCentralModel.Data.Raw_Data,...
            rr,...
            modIn,...
            runNTimes, 16) ;

        tmpCnsnsusCa = zeros([ kCentralModel.Data.n runNTimes ]) ;

        for jdx=1:runNTimes
            [~,tmpCnsnsusCa(:,jdx)] = community_assign(btcnsnsusModels(jdx).Model) ;
        end

        tmpAgreeMat = agreement(tmpCnsnsusCa) ./ runNTimes ;

        % get the consensus consitency
        btC(idx) = consensus_consistency(tmpAgreeMat) ;

        % make new kiter_prior for new loop
        for kdx=1:(kCentralModel.Data.n)

            btkiter_prior(:,kdx,idx+1) = ...
                sum(bsxfun(@eq,tmpCnsnsusCa(kdx,:), ...
                [1:(kCentralModel.R_Struct.k)]'),2) ./ cutoff ;

        end    

        % have we converged? or are we at the end of loop?
        if btC(idx) >= 0.95 || idx == CONSENSUS_ITER
            btconsensus_kiter_prior = btkiter_prior(:,:,idx+1) ;
            %[~,consensus_ca] = community_assign(consensus_kiter_prior) ; 
            btconsensus_kCentralModel = central_model(btcnsnsusModels) ;
            % also add the data back into consensus_kCentral...
            btconsensus_kCentralModel.Data = kCentralModel.Data ;
            break 
        end

    end
    
    %save results
    bt_C{btiter} = btC ;
    bt_kiter_prior_results{btiter} =  btconsensus_kiter_prior ;
    bt_consensus_model{btiter} = btconsensus_kCentralModel ;
                
    ppm1.increment() 
    
end

%% evaluate how similar all the bt_consensus_models are...

% function [ centralModel , simVec , simTriu ] = central_model( modelsStruct , priorMu , priorWeightTune)
[ bt_conCent , bt_conCent_simVec, bt_conCent_simTriu ] = central_model( bt_consensus_model );

%% first save results 

save('reproduce_bootstrapConsensus_script_results.mat',...
    'bt_conCent','bt_conCent_simVec','bt_conCent_simTriu',...
    'bt_C','bt_kiter_prior_results','bt_consensus_model')

%% evaluate how these all differ from the original model we got?

mi_btConMods_2_templateModel = zeros([50 1]) ;
vi_btConMods_2_templateModel = zeros([50 1]) ; 

% load consensus model
templateModel = load('/home/jfaskowi/JOSHSTUFF/projects/sbm3/data/interim/yeo_both_normalpoisson_a0p5_templateModel_1.mat') ;
templateModel = templateModel.templateModel ;
[~,tm_ca] = community_assign(templateModel) ;


TODODODODO

for idx=1:50
    
    [~,tmp_ca] = community_assign( bt_consensus_model{idx} );
    
   % just compute nmi from each bt_consensus model to the actual consensus
   % model...
    [vi_btConMods_2_templateModel(idx),mi_btConMods_2_templateModel(idx)] = ...
        partition_distance(tm_ca,tmp_ca) ;
    
end



