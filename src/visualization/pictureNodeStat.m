function [f1 , f2] = pictureNodeStat( nodeStat, plotWhat , colorStr, atlas, cb_lim) 
% inputs
% 1) node-wise statistics
% 2) what are they describing? (left, right, both?) only both rn
% 3) your only choice right now is yeo

if nargin < 2
    disp('need at least 2 args')
    return 
end

if ~exist('colorStr','var') || isempty(colorStr)
    colorStr = 'YlOrRd';
end

if ~exist('cb_lim','var') || isempty(cb_lim)
    cb_lim = '';
end

if nargin < 4
    atlas = 'yeo' ;
end

if ~strcmp(atlas,'yeo') 
    disp('your only choice is yeo right now')
    return
end

%% continue

if strcmp(plotWhat,'both')
    disp('plot both')
    numParcels = length(nodeStat) / 2; 
    
% elseif strcmp(plotWhat,'left')
%     disp('plot left')
%     numParcels = length(nodeStat) ;
%     
% elseif strcmp(plotWhat,'right')
%     disp('plot right')
%     numParcels = length(nodeStat) ;
    
else
    disp('invalid choice')
    return  
end

%% Mapping Networks or Modules
% load the vars I need
% load('4_josh_thanks.mat')
% load('4_josh_thanks2.mat')
load('ccs_surfstat_yeo_info.mat')

modules_lh = zeros(nVertices_lh,1);
modules_rh = zeros(nVertices_rh,1);

% addpath('figures/DrosteEffect-BrewerMap-04533de/')
cmap = brewermap(numParcels+2,colorStr);
%force 0 to be grey
%cmap(1,:) = [0.6 0.6 0.6];                        

for parcelID=1:numParcels

    disp(parcelID)

    if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'left'))
        tmpidx = parcel_vertexid_lh{parcelID};
        modules_lh(tmpidx) = nodeStat(parcelID);
    end
    
    if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'right'))
        tmpidx = parcel_vertexid_rh{parcelID};
        modules_rh(tmpidx) = nodeStat((parcelID)+numParcels);
    end
    
end
   
disp(plotWhat)

statRange = max(nodeStat) - min(nodeStat) ;
colorbarMod = 0.025 * statRange ;

if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'left'))
    
    disp('rendering left')
    
    medialwall_lh = medialwall_lh .* (min(nodeStat)-0.2);
    
    %render lh surfaces
    modules_lh = modules_lh + medialwall_lh;
    f1 = figure('Units', 'pixel', 'Position', [100 100 800 800]); 
    axis off
    SurfStatView(modules_lh, fs_lh, ' ', 'white', 'true'); 
    %colormap([0 0 0; 0.5 0.5 0.5; cmap]); 
    colormap([0 0 0; cmap]);
    
    if isempty(cb_lim)
        SurfStatColLim([min(nodeStat)-colorbarMod max(nodeStat)+colorbarMod]);
    else
        SurfStatColLim(cb_lim);
    end
    
    set(gcf, 'PaperPositionMode', 'auto');

end

%figout = [ana_dir '/figures/' group{i} '.modules.lh.png'];
%print('-dpng', '-r300', figout); close

if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'right'))
    
    disp('rendering right')
    
    %render rh surfaces
    modules_rh = modules_rh + medialwall_rh;
    f2 = figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(modules_rh, fs_rh, ' ', 'white', 'true'); 
    %colormap([0 0 0; 0.5 0.5 0.5; cmap]); 
    colormap([0 0 0; cmap]);
    
    if isempty(cb_lim)
        SurfStatColLim([min(nodeStat)-colorbarMod max(nodeStat)+colorbarMod]);
    else
        SurfStatColLim(cb_lim);
    end
    
    set(gcf, 'PaperPositionMode', 'auto');
end

%figout = [ana_dir '/figures/' group{i} '.modules.rh.png'];
%print('-dpng', '-r300', figout); close

