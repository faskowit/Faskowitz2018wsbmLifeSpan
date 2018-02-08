clc 
clearvars

config_file='config_template.m';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath(strcat(pwd,'/config'))
run(config_file);

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_templateModel_1.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/interim/', OUTPUT_STR, '_comVecs.mat');
load(loadName) ;

% the stats maps for wsbm, mod
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_comVec_dynamic_results.mat');
load(loadName) ;

loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
load(loadName) ;

%% before we run the regressions, have to deal with module exclude
 
wsbm_weiVec_cos2 = wsbm_weiVec_cos(~modExclude);
mod_weiVec_cos2 = mod_weiVec_cos(~modExclude) ;
wsbm_weiVec_cb2 = wsbm_weiVec_cb(~modExclude) ;
mod_weiVec_cb2 = mod_weiVec_cb(~modExclude) ;

age2 = datasetDemo.age(~modExclude) ;

%% run the regression for cosine similiarity
funcArgs = {1 500 [] 0 } ;
%funcArgs = {0 2 [] 0 } ;


fits = {'linear' 'quadratic' 'poisson'} ;

xVec = age2;
nSubj = length(age2) ;

regResultsCov_cos = cell([2 1]) ;
regResultsCovWMov_cos = cell([2 1]) ;
regResultsNoCov_cos = cell([2 1]) ;

regInputs = cell([2 1]);
% regInputs{1} = wsbm_weiVec_cos;
% regInputs{2} = mod_weiVec_cos;
regInputs{1} = wsbm_weiVec_cos2;
regInputs{2} = mod_weiVec_cos2;

for idx = 1:length(regInputs)
      
    disp(idx)
    
    % setup results for wsbm, mod... NO yeo....
    regResultsCov_cos{idx} = cell([length(fits) 1]);
    regResultsCovWMov_cos{idx} = cell([length(fits) 1]);
    regResultsNoCov_cos{idx} = cell([length(fits) 1]);

    Y = regInputs{idx}';
        
    sexCovar = (datasetDemo.sex(:,1)=='M') ;
    sexCovar = sexCovar(~modExclude) ;
    denCovar = totDensity(~modExclude) ;
    movCovar = datasetDemo.movement(~modExclude);
    
    % should we regress out age and totDensity, and then look at the similarity
    % or distance metrics? lets check this out!
    [ ~ , ~ , Y_nusReg ] = regress(Y,[ones(nSubj,1) sexCovar denCovar]);
    [ ~ , ~ , Y_nusReg_wMov ] = regress(Y,[ones(nSubj,1) sexCovar denCovar movCovar]);
    
    for jdx = 1:length(fits)
    
        disp(jdx)
        
        regressResults = struct() ; 
        [ regressResults.xvalR2 , ...
            regressResults.xvalsqErr, ... 
            regressResults.xvalYhat, ...
            regressResults.coef, ...
            regressResults.lsFitStruct,...
            regressResults.permStruct ] ...
            = nc_FitAndEvaluateModels(Y_nusReg,xVec,fits{jdx},funcArgs{:}) ;
        regResultsCov_cos{idx}{jdx} = regressResults ;
    
        regressResults = struct() ; 
        [ regressResults.xvalR2 , ...
            regressResults.xvalsqErr, ... 
            regressResults.xvalYhat, ...
            regressResults.coef, ...
            regressResults.lsFitStruct,...
            regressResults.permStruct ] ...
            = nc_FitAndEvaluateModels(Y_nusReg_wMov,xVec,fits{jdx},funcArgs{:}) ;
        regResultsCovWMov_cos{idx}{jdx} = regressResults ;
        
        regressResults = struct() ; 
        [ regressResults.xvalR2 , ...
            regressResults.xvalsqErr, ... 
            regressResults.xvalYhat, ...
            regressResults.coef, ...
            regressResults.lsFitStruct,...
            regressResults.permStruct ] ...
            = nc_FitAndEvaluateModels(Y,xVec,fits{jdx},funcArgs{:}) ;
        regResultsNoCov_cos{idx}{jdx} = regressResults ;
        
    end
end

%% run the regression for cb dist
funcArgs = {1 500 [] 0 } ;
%funcArgs = {0 2 [] 0 } ;

fits = {'linear' 'quadratic' 'poisson'} ;

xVec = age2;
nSubj = length(age2) ;
 
regResultsCov_cb = cell([2 1]) ;
regResultsCovWMov_cb = cell([2 1]) ;
regResultsNoCov_cb = cell([2 1]) ;

regInputs = cell([2 1]);
% regInputs{1} = wsbm_weiVec_cb;
% regInputs{2} = mod_weiVec_cb;
regInputs{1} = wsbm_weiVec_cb2;
regInputs{2} = mod_weiVec_cb2;
%regInputs{3} = yeo_weiVec_cb;

for idx = 1:length(regInputs)
        
    % setup results for wsbm, mod, 
    regResultsCov_cb{idx} = cell([length(fits) 1]);
    regResultsCovWMov_cb{idx} = cell([ length(fits) 1]);
    regResultsNoCov_cb{idx} = cell([length(fits) 1]);

    Y = regInputs{idx}';
    
    sexCovar = (datasetDemo.sex(:,1)=='M') ;
    sexCovar = sexCovar(~modExclude) ;
    denCovar = totDensity(~modExclude) ;
    movCovar = datasetDemo.movement(~modExclude);
    
    % should we regress out age and totDensity, and then look at the similarity
    % or distance metrics? lets check this out!
    [ ~ , ~ , Y_nusReg ] = regress(Y,[ones(nSubj,1) sexCovar denCovar]);
    [ ~ , ~ , Y_nusReg_wMov ] = regress(Y,[ones(nSubj,1) sexCovar denCovar movCovar]);
    
    for jdx = 1:length(fits)
    
         disp(jdx)
        
        regressResults = struct() ; 
        [ regressResults.xvalR2 , ...
            regressResults.xvalsqErr, ... 
            regressResults.xvalYhat, ...
            regressResults.coef, ...
            regressResults.lsFitStruct,...
            regressResults.permStruct ] ...
            = nc_FitAndEvaluateModels(Y_nusReg,xVec,fits{jdx},funcArgs{:}) ;
        regResultsCov_cb{idx}{jdx} = regressResults ;
    
        regressResults = struct() ; 
        [ regressResults.xvalR2 , ...
            regressResults.xvalsqErr, ... 
            regressResults.xvalYhat, ...
            regressResults.coef, ...
            regressResults.lsFitStruct,...
            regressResults.permStruct ] ...
            = nc_FitAndEvaluateModels(Y_nusReg_wMov,xVec,fits{jdx},funcArgs{:}) ;
        regResultsCovWMov_cb{idx}{jdx} = regressResults ;
        
        regressResults = struct() ; 
        [ regressResults.xvalR2 , ...
            regressResults.xvalsqErr, ... 
            regressResults.xvalYhat, ...
            regressResults.coef, ...
            regressResults.lsFitStruct,...
            regressResults.permStruct ] ...
            = nc_FitAndEvaluateModels(Y,xVec,fits{jdx},funcArgs{:}) ;
        regResultsNoCov_cb{idx}{jdx} = regressResults ;
        
    end
end

%% save

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_comVec_dynamic_regResults.mat');
save(outName,...
    'regResults*',...
    ...
    '-v7.3')
