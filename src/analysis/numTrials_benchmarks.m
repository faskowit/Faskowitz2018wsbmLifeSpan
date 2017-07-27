%% test how many iterations we need

data = avgTemp ;

% different options we should consider playing with
% --Inference Options:
%   'numTrials' - number of trials with different random initial conditions
%   'algType' - 'vb'* naive bayes, 'bp' belief propagation
%   'networkType' - 'dir'* directed,'sym' = symmetric,'asym' = asymmetric
%   'nanType' - 'missing'* nans are missing, 'nonedge' = nans are nonedge 
%   'mainMaxIter' - Maximum number of iterations in main_loop
%   'mainTol' - Minimum (Max Norm) convergence tolerance in main_loop
%   'muMaxIter' - Maximum number of iterations in mu_loop
%   'muTol' - Minimum (Max Norm) convergence tolerance in mu_loop

%     % Setup Default Options Struct
%     Options = struct('algType','vb',...
%                      'alpha',0.5,...
%                      'networkType','directed',...
%                      'nanType','missing',...
%                      'verbosity',1,...
%                      'numTrials',50,...
%                      'mainMaxIter',80,...
%                      'mainTol',0.001,...
%                      'muMaxIter',50,...
%                      'muTol',0.001,...
%                      'mu_0',[]);
        
%% first lets just see effect of numTrials?
% test 1

% levels = 10 ;
% repeat_iters = 25 ; 
% 
% numTrials_to_test = floor(logspace(1.6990,3.6990,levels)) ; 
% numTrials_results = cell([levels repeat_iters]);
% numTrials_time_results = zeros([levels repeat_iters]) ;

%%
            %'E_Distr', setup_distr('poisson',[0.1,0.01]), ...

modelInputs = {'W_Distr', 'Normal', ...
            'E_Distr', setup_distr('poisson',[0.1,0.01]), ...
            'alpha', 0.5 , ...
            'mainMaxIter', 100 , ...
            'mainTol',0.1, ...
            'muMaxIter', 50,  ...
            'muTol', 0.1, ...
            'verbosity' , 0};

%            'numTrials', 150 , ...
        
r_struct = sym_RStruct(7) ;

%% running just once

% tic
% [~,oneRun_results] = wsbm(data, ...
%     r_struct, ...
%     modelInputs{:},...
%     'numTrials', 50 ,...
%     'verbosity' , 1) ;  
% toc

%%
parallel_pool = gcp ; 
for idx=1:levels
    
    parfor iidx=1:repeat_iters

            t1=tic
            
            [~,numTrials_results{idx,iidx}] = wsbm(data, ...
                r_struct, ...
                modelInputs{:},...
                'numTrials', numTrials_to_test(idx) ) ;  
            
            numTrials_time_results(idx,iidx) = toc(t1) ;
            
    end

end

save('analysis/numTrials_results2.mat','numTrials_results') ;

numTrialsTest_logE = cellfun(@(a)a.Para.LogEvidence,numTrials_results) ;

%% now run more from the starting point of already high logEvvid runs...
% use results of 2000 as example 

nt_2000 = numTrials_results(8,:) ;
nt_2000_logE = cellfun(@(a)a.Para.LogEvidence,nt_2000) ;
[~,sortI] = sort(nt_2000_logE,'descend') ;
top5_models_idx = sortI(1:5) ;

modelInputs = {'W_Distr', 'Normal', ...
            'E_Distr', setup_distr('poisson',[0.1,0.01]), ...
            'alpha', 0.5 , ...
            'mainMaxIter', 100 , ...
            'mainTol',0.01, ...
            'muMaxIter', 50,  ...
            'muTol', 0.01, ...
            'verbosity' , 0};

% tic
% [~,oneRun_results] = wsbm(data, ...
%     r_struct, ...
%     modelInputs{:},...
%     'numTrials', 50 ,...
%     'verbosity' , 2,...
%     'mu_0', muPrior) ;  
% toc

numTrials_results_additional = cell([5 repeat_iters]) ;
numTrials_time_results2 = zeros([5 repeat_iters]) ;

for idx=1:5

    muPrior = make_WSBM_prior(nt_2000{top5_models_idx(idx)},2) ;
    
    parfor iidx=1:repeat_iters

            t1=tic
            
            [~,numTrials_results_additional{idx,iidx}] = wsbm(data, ...
                r_struct, ...
                modelInputs{:},...
                'numTrials', 500,...
                'mu_0', muPrior) ;  
            
            numTrials_time_results2(idx,iidx) = toc(t1) ;
            
    end

end

additional_logE = cellfun(@(a)a.Para.LogEvidence,numTrials_results_additional) ;

%% test 2

levels = 10 ;
repeat_iters = 25 ; 

numTrials_to_test = floor(logspace(1.6990,3.6990,levels)) ; 
numTrials_results = cell([levels repeat_iters]);
numTrials_time_results = zeros([levels repeat_iters]) ;

modelInputs = {'W_Distr', 'Normal', ...
            'E_Distr', setup_distr('poisson',[0.1,0.01]), ...
            'alpha', 0.5 , ...
            'mainMaxIter', 100 , ...
            'mainTol',0.1, ...
            'muMaxIter', 50,  ...
            'muTol', 0.1, ...
            'verbosity' , 0};

%            'numTrials', 150 , ...
        
r_struct = sym_RStruct(9) ;

%% running just once

% tic
% [~,oneRun_results] = wsbm(data, ...
%     r_struct, ...
%     modelInputs{:},...
%     'numTrials', 100 ) ;  
% toc

%%
parallel_pool = gcp ; 
for idx=1:levels
    
    parfor iidx=1:repeat_iters

            t1=tic
            
            [~,numTrials_results{idx,iidx}] = wsbm(data, ...
                r_struct, ...
                modelInputs{:},...
                'numTrials', numTrials_to_test(idx) ) ;  
            
            numTrials_time_results(idx,iidx) = toc(t1) ;
            
    end

end

save('analysis/numTrials_results2.mat','numTrials_results') ;

numTrialsTest_logE2 = cellfun(@(a)a.Para.LogEvidence,numTrials_results) ;

%% add more trials

numTrials_to_test = [ numTrials_to_test 6000 7000 ] ;

parallel_pool = gcp ; 
for idx=11:12
    
    parfor iidx=1:repeat_iters

            t1=tic
            
            [~,numTrials_results{idx,iidx}] = wsbm(data, ...
                r_struct, ...
                modelInputs{:},...
                'numTrials', numTrials_to_test(idx) ) ;  
            
            numTrials_time_results(idx,iidx) = toc(t1) ;
            
    end

end

%% make some plots

subplot(1,3,1)
scatter(numTrials_to_test,median(numTrialsTest_logE2,2))
ylabel('log evidence, dots: median, line: mean')
xlabel('num trials')
hold on
plot(numTrials_to_test,mean(numTrialsTest_logE2,2))
hold off

subplot(1,3,2)
scatter(numTrials_to_test,std(numTrialsTest_logE2,0,2))
ylabel('log evidence std')
xlabel('num trials')
hold on
plot(numTrials_to_test,std(numTrialsTest_logE2,0,2))
hold off

subplot(1,3,3)
scatter(numTrials_to_test,median(numTrials_time_results,2)./3600)
ylabel('hours runtime, dots: median, line: mean')
xlabel('num trials')
hold on
plot(numTrials_to_test,mean(numTrials_time_results,2)./3600)
hold off

%%
% 
% mmm = make_WSBM_prior(oneRun_results.Para.mu,5) ;
% 
% mu_seed = rand(9,114).* make_WSBM_prior(oneRun_results.Para.mu,8);
% mu_seed = (mu_seed./(ones(9,1)*max(mu_seed,[],1))).^ 5;
% theSeed = mu_seed./(ones(9,1)*sum(mu_seed,1))


%% test 3

levels = 10 ;
repeat_iters = 25 ; 

numTrials_to_test = floor(logspace(1.6990,3.6990,levels)) ; 
numTrials_results = cell([levels repeat_iters]);
numTrials_time_results = zeros([levels repeat_iters]) ;

%%

modelInputs = {'W_Distr', 'Normal', ...
            'E_Distr', setup_distr('dcpoisson',[0.1,0.01]), ...
            'alpha', 0.5 , ...
            'mainMaxIter', 100 , ...
            'mainTol',0.1, ...
            'muMaxIter', 50,  ...
            'muTol', 0.1, ...
            'verbosity' , 0};

%            'numTrials', 150 , ...
        
r_struct = sym_RStruct(9) ;

%% running just once

tic
[~,oneRun_results] = wsbm(data, ...
    r_struct, ...
    modelInputs{:},...
    'numTrials', 50 ,...
    'verbosity', 1) ;  
toc

%%
parallel_pool = gcp ; 
for idx=1:levels
    
    parfor iidx=1:repeat_iters

            t1=tic
            
            [~,numTrials_results{idx,iidx}] = wsbm(data, ...
                r_struct, ...
                modelInputs{:},...
                'numTrials', numTrials_to_test(idx) ) ;  
            
            numTrials_time_results(idx,iidx) = toc(t1) ;
            
    end

end

save('analysis/numTrials_results2.mat','numTrials_results') ;

numTrialsTest_logE2 = cellfun(@(a)a.Para.LogEvidence,numTrials_results) ;