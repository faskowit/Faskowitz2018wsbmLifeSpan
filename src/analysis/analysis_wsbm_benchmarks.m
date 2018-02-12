%% load necessary data

clc
clearvars

% need to edit config file string to match what you want!!
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

load('data/interim/yeo_both_normalpoisson_a0p5_templateModel_1.mat')

%% test how many iterations we need

dataMat = templateModel.Data.Raw_Data ;
   
%% first lets just see effect of numTrials?
% test 1

levels = 10 ;
repeat_iters = 10 ; 

numTrials_to_test = floor(logspace(1.6990,3.39795,levels)) ; 
numTrials_results = cell([levels repeat_iters]);
numTrials_time_results = zeros([levels repeat_iters]) ;

%% setup input

modelInputs = {'W_Distr', 'Normal', ...
            'E_Distr', setup_distr('poisson',[0.1,0.01]), ...
            'alpha', 0.5 , ...
            'mainMaxIter', 100 , ...
            'muMaxIter', 50,  ...
            'verbosity' , 0};
        
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
ppm1 = ParforProgMon('indivFits',(repeat_iters * levels),1) ;

for idx=1:levels
    
    parfor iidx=1:repeat_iters

            t1=tic
            
            [~,numTrials_results{idx,iidx}] = wsbm(dataMat, ...
                r_struct, ...
                modelInputs{:},...
                'numTrials', numTrials_to_test(idx) ) ;  
            
            numTrials_time_results(idx,iidx) = toc(t1) ; 
    
            ppm1.increment() 
    end
    
end

%% now run more from the starting point of already high logEvvid runs...
% use results of 2000 as example 

boostBegin = numTrials_results(4,:) ;
boostBeginLogE = cellfun(@(a)a.Para.LogEvidence,boostBegin) ;
[~,sortI] = sort(boostBeginLogE,'descend') ;
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
numTrials_time_additional = zeros([5 repeat_iters]) ;

ppm2 = ParforProgMon('boosting',(repeat_iters * 5),1) ;

for idx=1:5

    muPrior = make_WSBM_prior(boostBegin{top5_models_idx(idx)},2) ;
    
    parfor iidx=1:repeat_iters

            t1=tic
            
            [~,numTrials_results_additional{idx,iidx}] = wsbm(dataMat, ...
                r_struct, ...
                modelInputs{:},...
                'numTrials', 100,...
                'mu_0', muPrior) ;  
            
            numTrials_time_additional(idx,iidx) = toc(t1) ;
            
            ppm2.increment() 
    end

end

% TODO FIX THIS MEASUREMENT
additional_logE = cellfun(@(a)a.Para.LogEvidence,numTrials_results_additional) ;

%% plot stuff

numTrials_to_test_unrollable = repmat(numTrials_to_test',1,10) ;

% time it takes
scatter(numTrials_to_test_unrollable(:),numTrials_time_results(:))

% get the absolute logEvidence from the results of first test
numTrials_results_abs = cell(size(numTrials_results)) ;

for idx = 1:size(numTrials_results,1)
    for jdx = 1:size(numTrials_results,2)

        [~,tmp_harsh] = make_WSBM_prior(numTrials_results{idx,jdx});
        
        modelInputs = {'W_Distr', numTrials_results{1,1}.W_Distr, ...
            'E_Distr', numTrials_results{1,1}.E_Distr, ...
            'alpha', numTrials_results{1,1}.Options.alpha, ...
            'mainMaxIter', 1 , ...
            'muMaxIter', 50,  ...
            'mu_0', tmp_harsh, ...
            'verbosity' , 0};

         [~,tmp_model] = wsbm(dataMat, ...
            numTrials_results{1,1}.R_Struct.R, ...
            modelInputs{:},...
            'numTrials', 1) ;  
        
        numTrials_results_abs{idx,jdx} = tmp_model ;
        
    end
end
    
tmp1 = cellfun(@(a)a.Para.LogEvidence,numTrials_results)
tmp2 = cellfun(@(a)a.Para.LogEvidence,numTrials_results_abs)

% %% test 2
% 
% levels = 10 ;
% repeat_iters = 10 ; 
% 
% numTrials_to_test = floor(logspace(1.6990,3.6990,levels)) ; 
% numTrials_results = cell([levels repeat_iters]);
% numTrials_time_results = zeros([levels repeat_iters]) ;
% 
% modelInputs = {'W_Distr', 'Normal', ...
%             'E_Distr', setup_distr('poisson',[0.1,0.01]), ...
%             'alpha', 0.5 , ...
%             'mainMaxIter', 100 , ...
%             'mainTol',0.1, ...
%             'muMaxIter', 50,  ...
%             'muTol', 0.1, ...
%             'verbosity' , 0};
% 
% %            'numTrials', 150 , ...
%         
% r_struct = sym_RStruct(9) ;
% 
% %% running just once
% 
% % tic
% % [~,oneRun_results] = wsbm(data, ...
% %     r_struct, ...
% %     modelInputs{:},...
% %     'numTrials', 100 ) ;  
% % toc
% 
% %%
% parallel_pool = gcp ; 
% for idx=1:levels
%     
%     parfor iidx=1:repeat_iters
% 
%             t1=tic
%             
%             [~,numTrials_results{idx,iidx}] = wsbm(dataMat, ...
%                 r_struct, ...
%                 modelInputs{:},...
%                 'numTrials', numTrials_to_test(idx) ) ;  
%             
%             numTrials_time_results(idx,iidx) = toc(t1) ;
%             
%     end
% 
% end
% 
% save('analysis/numTrials_results2.mat','numTrials_results') ;
% 
% numTrialsTest_logE2 = cellfun(@(a)a.Para.LogEvidence,numTrials_results) ;
% 
% %% add more trials
% 
% numTrials_to_test = [ numTrials_to_test 6000 7000 ] ;
% 
% parallel_pool = gcp ; 
% for idx=11:12
%     
%     parfor iidx=1:repeat_iters
% 
%             t1=tic
%             
%             [~,numTrials_results{idx,iidx}] = wsbm(dataMat, ...
%                 r_struct, ...
%                 modelInputs{:},...
%                 'numTrials', numTrials_to_test(idx) ) ;  
%             
%             numTrials_time_results(idx,iidx) = toc(t1) ;
%             
%     end
% 
% end
% 
% %% make some plots
% 
% subplot(1,3,1)
% scatter(numTrials_to_test,median(numTrialsTest_logE2,2))
% ylabel('log evidence, dots: median, line: mean')
% xlabel('num trials')
% hold on
% plot(numTrials_to_test,mean(numTrialsTest_logE2,2))
% hold off
% 
% subplot(1,3,2)
% scatter(numTrials_to_test,std(numTrialsTest_logE2,0,2))
% ylabel('log evidence std')
% xlabel('num trials')
% hold on
% plot(numTrials_to_test,std(numTrialsTest_logE2,0,2))
% hold off
% 
% subplot(1,3,3)
% scatter(numTrials_to_test,median(numTrials_time_results,2)./3600)
% ylabel('hours runtime, dots: median, line: mean')
% xlabel('num trials')
% hold on
% plot(numTrials_to_test,mean(numTrials_time_results,2)./3600)
% hold off
% 
% %%
% % 
% % mmm = make_WSBM_prior(oneRun_results.Para.mu,5) ;
% % 
% % mu_seed = rand(9,114).* make_WSBM_prior(oneRun_results.Para.mu,8);
% % mu_seed = (mu_seed./(ones(9,1)*max(mu_seed,[],1))).^ 5;
% % theSeed = mu_seed./(ones(9,1)*sum(mu_seed,1))


% for idx = 1:100 
% 
%     tmp_mu = numTrials_results{10,1}.Para.mu    ;
%     tmp_mu = tmp_mu(:,randperm(114)) ;
% 
%     [~,tmp_mu] = make_WSBM_prior(tmp_mu,0) ;
%     
%     modelInputs = {'W_Distr', numTrials_results{1,1}.W_Distr, ...
%         'E_Distr', numTrials_results{1,1}.E_Distr, ...
%         'alpha', numTrials_results{1,1}.Options.alpha, ...
%         'mainMaxIter', 1 , ...
%         'muMaxIter', 50,  ...
%         'mu_0', tmp_mu, ...
%         'verbosity' , 0};
%     [~,tmp_model] = wsbm(dataMat, ...
%         numTrials_results{1,1}.R_Struct.R, ...
%         modelInputs{:},...
%         'numTrials', 1) ; 
% 
%     wwwwww(idx) = tmp_model.Para.LogEvidence ;
%     
% end

for Rval = 8:13     
        
    r_struct = sym_RStruct(Rval) ;

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
    ppm1 = ParforProgMon('indivFits',(repeat_iters * levels),1) ;

    for idx=1:levels

        parfor iidx=1:repeat_iters

                t1=tic

                [~,numTrials_results{idx,iidx}] = wsbm(dataMat, ...
                    r_struct, ...
                    modelInputs{:},...
                    'numTrials', numTrials_to_test(idx) ) ;  

                numTrials_time_results(idx,iidx) = toc(t1) ; 

                ppm1.increment() 
        end

    end

end