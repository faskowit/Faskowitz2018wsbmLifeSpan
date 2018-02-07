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

% the stats maps for wsbm, mod, yeo
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_comVec_static_results.mat');
load(loadName) ;

% the stats maps for wsbm, mod, yeo
loadName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_basicData_v7p3.mat');
load(loadName) ;

FIGURE_NAME = 'figE' ;

outputdir = strcat(PROJECT_DIR,'/reports/figures/',FIGURE_NAME,'/');
mkdir(outputdir)

%% run the regression for cosine similiarity
funcArgs = {1 500 [] 0 } ;
%funcArgs = {0 2 [] 0 } ;

nSubj = length(datasetDemo.age) ;

fits = {'linear' 'quadratic' 'poisson'} ;

xVec = datasetDemo.age;

regResultsCov_cos = cell([3 1]) ;
regResultsCovWMov_cos = cell([3 1]) ;
regResultsNoCov_cos = cell([3 1]) ;

regInputs = cell([3 1]);
regInputs{1} = wsbm_weiVec_cos;
regInputs{2} = mod_weiVec_cos;
regInputs{3} = yeo_weiVec_cos;

for idx = 1:length(regInputs)
        
    % setup results for wsbm, mod, yeo....
    regResultsCov_cos{idx} = cell([length(fits) 1]);
    regResultsCovWMov_cos{idx} = cell([length(fits) 1]);
    regResultsNoCov_cos{idx} = cell([length(fits) 1]);

    Y = regInputs{idx}';
    
    % should we regress out age and totDensity, and then look at the similarity
    % or distance metrics? lets check this out!
    [ ~ , ~ , Y_nusReg ] = regress(Y,[ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity]);
    [ ~ , ~ , Y_nusReg_wMov ] = regress(Y,[ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity datasetDemo.movement]);
    
    for jdx = 1:length(fits)
    
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

xVec = datasetDemo.age;
% 
regResultsCov_cb = cell([3 1]) ;
regResultsCovWMov_cb = cell([3 1]) ;
regResultsNoCov_cb = cell([3 1]) ;

regInputs = cell([3 1]);
regInputs{1} = wsbm_weiVec_cb;
regInputs{2} = mod_weiVec_cb;
regInputs{3} = yeo_weiVec_cb;

for idx = 1:length(regInputs)
        
    % setup results for wsbm, mod, yeo....
    regResultsCov_cb{idx} = cell([length(fits) 1]);
    regResultsCovWMov_cb{idx} = cell([ length(fits) 1]);
    regResultsNoCov_cb{idx} = cell([length(fits) 1]);

    Y = regInputs{idx}';
    
    % should we regress out age and totDensity, and then look at the similarity
    % or distance metrics? lets check this out!
    [ ~ , ~ , Y_nusReg ] = regress(Y,[ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity]);
    [ ~ , ~ , Y_nusReg_wMov ] = regress(Y,[ones(nSubj,1) (datasetDemo.sex(:,1)=='M') totDensity datasetDemo.movement]);
    
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

outName = strcat(OUTPUT_DIR, '/processed/', OUTPUT_STR, '_comVec_static_regResults.mat');
save(outName,...
    'regResults*',...
    ...
    '-v7.3')

%% figure out the fits that are best

all_reg_results = cell([6 1]) ;
all_reg_results{1} = regResultsCov_cos ;
all_reg_results{2} = regResultsCovWMov_cos ;
all_reg_results{3} = regResultsNoCov_cos ;
all_reg_results{4} = regResultsCov_cb ;
all_reg_results{5} = regResultsCovWMov_cb ;
all_reg_results{6} = regResultsNoCov_cb ;

reg_result_names = {'Cos sim. w/ covar.' 'Cos sim. w/ covar+mov' 'Cos sim' ...
    'CB dist. w/ covar' 'CB dist. w/ covar+mov' 'CB dist'};

reg_ylabel_names = { 'Cos similarity adj. for covar' 'Cos similarity adj. for covar, mov' 'Cos similarity' ...
    'CB distance adj. covar' 'CB distance adj. for covar, mov' 'CB distance'} ;

default_cmap = [0    0.4470    0.7410 ;
                0.8500    0.3250    0.0980;
                0.9290    0.6940    0.1250 ] ;

subp = tight_subplot(2,3,[.10 .05],[.1 .05],[.1 .05]) ;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);      

for idx = 1:length(all_reg_results)
    
    currentRegResults = all_reg_results{idx} ;

    [~,wsbm_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults{1}));
    [~,mod_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults{2}));
    [~,yeo_minIdx] = min(cellfun(@(x) sqrt(mean(x.xvalsqErr)), currentRegResults{3}));

    wsbm_trend = currentRegResults{1}{wsbm_minIdx};
    mod_trend = currentRegResults{2}{mod_minIdx};
    yeo_trend = currentRegResults{3}{yeo_minIdx};
    
%% plot it

    %fig1 = figure ;
    axes(subp(idx)) 
    
    % yeo
    pl = viz_blockRegress(yeo_trend,0,default_cmap(3,:));
    set(pl,'MarkerFaceColor', default_cmap(3,:),...
        'MarkerEdgeColor',default_cmap(3,:),...
        'MarkerFaceAlpha',.15,...
        'MarkerEdgeAlpha',.1) 
    hold

    % mod
    pl = viz_blockRegress(mod_trend,0,default_cmap(2,:));
    set(pl,'MarkerFaceColor', default_cmap(2,:),...
        'MarkerEdgeColor',default_cmap(2,:),...
        'MarkerFaceAlpha',.15,...
        'MarkerEdgeAlpha',.1) 
    hold

    % wsbm
    pl = viz_blockRegress(wsbm_trend,0, default_cmap(1,:)) ;
    set(pl,'MarkerFaceColor',default_cmap(1,:),...
        'MarkerFaceColor',default_cmap(1,:),...
        'MarkerFaceAlpha',.15,...
        'MarkerEdgeAlpha',.1)
    hold
   
    xlim([min(pl.XData)-3 max(pl.XData)+3])
    
    ylabel(reg_ylabel_names{idx})
    xlabel('Age')
    
    title(reg_result_names{idx})
    
    yrange = ylim ;
    yrangeAbs = yrange(2) - yrange(1) ;
    % add the R2!! 
    if idx < 4
        ypos = yrangeAbs * 0.25;
        ypos = yrange(1) + ypos ;
    else
        ypos = yrangeAbs * 0.95;
        ypos = yrange(1) + ypos ;
    end

    annotText = { strcat('wsbm R^2:',32,num2str(round(wsbm_trend.xvalR2 / 100,3))) ...
        strcat('mod R^2:',32,num2str(round(mod_trend.xvalR2 / 100,3))) ...
        strcat('yeo R^2:',32,num2str(round(yeo_trend.xvalR2 / 100,3))) ...
        } ;

    text((min(pl.XData)+2),ypos, annotText,'FontSize',10,'VerticalAlignment','cap')    
    
end

%% but also just visualize the similarity and distance...





%% extra viz for communities...

y_lim = [0.5 1] ;

for idx=1:10
   subplot(2,5,idx) 
   plot(fitlm(datasetDemo.age,wsbm_comms_weiVec_cos(idx,:),'quadratic'));
   %ylim([0 0.5])
   ylim(y_lim)
end
suptitle('wsbm comm cos')

figure
for idx=1:10
   subplot(2,5,idx) 
   plot(fitlm(datasetDemo.age,mod_comms_weiVec_cos(idx,:),'quadratic'));
   %ylim([0 0.5])
   ylim(y_lim)
end
suptitle('mod comm cos')

% figure;
% for idx=1:10
%    subplot(2,5,idx) 
%    plot(fitlm(datasetDemo.age,wsbm_comms_weiVec_eud(idx,:),'quadratic'));
%    
% end

%%

y_lim = [0 0.75] ;

for idx=1:10
   subplot(2,5,idx) 
   plot(fitlm(datasetDemo.age,wsbm_comms_weiVec_cb(idx,:),'quadratic'));
   ylim(y_lim)
end
suptitle('wsbm comm eud')

figure
for idx=1:10
   subplot(2,5,idx) 
   plot(fitlm(datasetDemo.age,mod_comms_weiVec_cb(idx,:),'quadratic'));
   ylim(y_lim)
end
suptitle('mod comm eud')

