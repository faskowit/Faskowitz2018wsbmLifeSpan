%% load data

load('/home/jfaskowi/JOSHSTUFF/projects/sbm3/reports/reproduce_consistencyConsensus_script_data.mat')

%% test conistency of consensus method 

%init a cell array to hold all boot iters results...
cc_kiter_prior_results = cell([50 1]) ;
cc_consensus_model = cell([50 1]) ;
cc_C = cell([50 1]) ;

for cciter=1:50

    % iterate over each node
    for idx=1:(kCentralModel.Data.n)

        cckiter_prior(:,idx,1) = ...
            sum(bsxfun(@eq,kBestModels_ca_algn_sort(idx,1:cutoff), ...
            [1:(kCentralModel.R_Struct.k)]'),2) ./ cutoff ;

    end  

    CONSENSUS_ITER = 10 ;
    ccC = zeros([CONSENSUS_ITER 1]) ;

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
            'mu_0', cckiter_prior(:,:,idx)};

        % function [ allModels ] = wsbmFitNTimes( adjMat, rStruct , modelInputs , numFits , numCores)
        cccnsnsusModels = wsbmFitNTimes(kCentralModel.Data.Raw_Data,...
            rr,...
            modIn,...
            runNTimes, 16) ;

        tmpCnsnsusCa = zeros([ kCentralModel.Data.n runNTimes ]) ;

        for jdx=1:runNTimes
            [~,tmpCnsnsusCa(:,jdx)] = community_assign(cccnsnsusModels(jdx).Model) ;
        end

        tmpAgreeMat = agreement(tmpCnsnsusCa) ./ runNTimes ;

        % get the consensus consitency
        ccC(idx) = consensus_consistency(tmpAgreeMat) ;

        % make new kiter_prior for new loop
        for kdx=1:(kCentralModel.Data.n)

            cckiter_prior(:,kdx,idx+1) = ...
                sum(bsxfun(@eq,tmpCnsnsusCa(kdx,:), ...
                [1:(kCentralModel.R_Struct.k)]'),2) ./ cutoff ;

        end    

        % have we converged? or are we at the end of loop?
        if ccC(idx) >= 0.95 || idx == CONSENSUS_ITER
            ccconsensus_kiter_prior = cckiter_prior(:,:,idx+1) ;
            %[~,consensus_ca] = community_assign(consensus_kiter_prior) ; 
            ccconsensus_kCentralModel = central_model(cccnsnsusModels) ;
            % also add the data back into consensus_kCentral...
            ccconsensus_kCentralModel.Data = kCentralModel.Data ;
            break 
        end

    end
    
    %save results
    cc_C{cciter} = ccC ;
    cc_kiter_prior_results{cciter} =  ccconsensus_kiter_prior ;
    cc_consensus_model{cciter} = ccconsensus_kCentralModel ;
    
end

%% evaluate how similar all the bt_consensus_models are...

% function [ centralModel , simVec , simTriu ] = central_model( modelsStruct , priorMu , priorWeightTune)
[ cc_conCent , cc_conCent_simVec, cc_conCent_simTriu ] = central_model( cc_consensus_model(1:50) );

%% save results

save('reports/reproduce_consistencyConsensus_script_results.mat',...
    'cc_conCent','cc_conCent_simVec','cc_conCent_simTriu','cc_C',...
    'cc_kiter_prior_results','cc_consensus_model')