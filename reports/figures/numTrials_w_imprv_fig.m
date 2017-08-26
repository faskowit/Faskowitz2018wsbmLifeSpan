% load the results 
addpath('/home/jfaskowi/JOSHSTUFF/projects/sbm3/src/external/WSBM_v1.2/')
load('/home/jfaskowi/JOSHSTUFF/projects/SBM2/analysis/numTrials_results_poisson7.mat')
load('/home/jfaskowi/JOSHSTUFF/projects/SBM2/workspaces/boosting_test_15iters_wProgPrior_0p5.mat')
load('/home/jfaskowi/JOSHSTUFF/projects/SBM2/workspaces/boosting_test_15iters_wProgPrior_0p5_rH.mat')
data = numTrials_results{1,1}.Data.Raw_Data ;

%% plot what we just ran above

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first plot baseline results %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
cm = parula(length(prirPtrn)-1);

logE_numTrials = cellfun(@(a)a.Para.LogEvidence,numTrials_results) ;

%nT = log(numTrials_to_test);
nT = numTrials_to_test;

for idx=1:size(logE_numTrials,2)

    scatter(nT,logE_numTrials(:,idx),[],[0 0.5 1])
    hold on
    
end

scatter(nT,median(logE_numTrials,2),135,[0.1 0.5 1],'filled')
ylabel('log evidence, dots: median, line: mean')
xlabel('num trials')
% hold on
% plot(numTrials_to_test,mean(logE_numTrials,2))
%hold off
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now plot the boosted iters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx=1:length(nT)
    
    tmp_results = results_here{idx} ;
    tmp_results = squeeze(struct2cell(tmp_results));
    
    %loop through the 5 iterations 
    %for iidx=1:20
    for iidx=1:15
       
        tmp_vals = zeros([25 1]) ;
        
        for iiidx=1:25
            
            tmp_vals(iiidx) = tmp_results{iiidx}(iidx).Para.LogEvidence ;
         
        end 
        
        % plot it
%         scatter(numTrials_to_test(idx) .* (ones([length(tmp_vals) 1])),...
%             tmp_vals(:),...
%             [],[ 1 (iidx)*0.05 (iidx)*0.1 ]) ;
%         hold on
%         scatter(numTrials_to_test(idx), median(tmp_vals),135,...
%             [ 1 (iidx)*0.05 (iidx)*0.1 ],'filled')

        % plotting the datapoints
        scatter(nT(idx) .* (ones([length(tmp_vals) 1])),...
            tmp_vals(:),...
            [], cm(iidx,:) ) ;
        hold on
        
        % plotting also the center 
        scatter(nT(idx), median(tmp_vals),135,...
            cm(iidx,:),'filled')
        
        hold on
        
    end
end

hold off

%% also get imporvement

resultsLooking = results_here{8} ;

for idx=1:25

    tmp_results = resultsLooking(idx) ;

    for iidx=1:(length(prirPtrn)-1)

       log_ediv_runs(idx,iidx) = tmp_results.ModelRun(iidx).Para.LogEvidence ;

       % compute improvement
       if iidx == 1
            continue
       else
            improvemnt(idx,iidx) =  tmp_results.ModelRun(iidx).Para.LogEvidence ...
                - tmp_results.ModelRun(iidx-1).Para.LogEvidence ;
       end
       
    end

end
