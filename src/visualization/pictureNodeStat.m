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

% addpath('figures/DrosteEffect-BrewerMap-04533de/')
cmap = brewermap(100,colorStr);
%force 0 to be grey
%cmap(1,:) = [0.6 0.6 0.6];                        

statRange = max(nodeStat) - min(nodeStat) ;
colorbarMod = 0.01 * statRange ;

modules_lh = ones(nVertices_lh,1);
modules_rh = ones(nVertices_rh,1);

if min(nodeStat) < 0
    mm = -( max(abs([min(nodeStat) max(nodeStat)]))) ; 
    mx = ( max(abs([min(nodeStat) max(nodeStat)]))) ; 
    disp('zero in nodeStat')
else
    mm = min(nodeStat) ;
    mx = max(nodeStat) ;
end

% modules_lh = (min(nodeStat) - (2*colorbarMod)) .* modules_lh;
% modules_rh = (min(nodeStat) - (2*colorbarMod)) .* modules_rh;
modules_lh = (mm - (2*colorbarMod)) .* modules_lh;
modules_rh = (mm - (2*colorbarMod)) .* modules_rh;

% if there are NaNs, make it border
nodeStat(isnan(nodeStat)) = (mm-(1*colorbarMod)) ;

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
   
%disp(plotWhat)

if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'left'))
    
    disp('rendering left')
    
    % clear any vars on medial wall
    modules_lh = modules_lh .* ~(medialwall_lh) ;
    
    %medialwall_lh = medialwall_lh .* (min(nodeStat)-(1*colorbarMod));
    medialwall_lh = medialwall_lh .* (mm-(1*colorbarMod));
    
    %render lh surfaces
    modules_lh = modules_lh + medialwall_lh;
    f1 = figure('Units', 'pixel', 'Position', [100 100 800 800]); 
    axis off
    SurfStatView(modules_lh, fs_lh, ' ', 'white', 'true'); 
    colormap([0 0 0; 0.5 0.5 0.5;  cmap]); 
    %colormap([0 0 0; cmap]);
    
    if isempty(cb_lim)
        %SurfStatColLim([min(nodeStat)-(3*colorbarMod) max(nodeStat)+(3*colorbarMod)]);
        %SurfStatColLim([min(nodeStat)-(3*colorbarMod) max(nodeStat)+(3*colorbarMod)]);
        SurfStatColLim([mm-(3*colorbarMod) mx]);
    else
        SurfStatColLim(cb_lim);
    end
    
    set(gcf, 'PaperPositionMode', 'auto');

end

%figout = [ana_dir '/figures/' group{i} '.modules.lh.png'];
%print('-dpng', '-r300', figout); close

if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'right'))
    
    disp('rendering right')
    
    % clear any vars on medial wall
    modules_rh = modules_rh .* ~(medialwall_rh) ;
    
%     medialwall_rh = medialwall_rh .* (min(nodeStat)-(1*colorbarMod));
    medialwall_rh = medialwall_rh .* (mm - (1*colorbarMod));
    
    %render rh surfaces
    modules_rh = modules_rh + medialwall_rh;
    f2 = figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(modules_rh, fs_rh, ' ', 'white', 'true'); 
    colormap([0 0 0; 0.5 0.5 0.5;   cmap]); 
    %colormap([0 0 0; cmap]);
    
    if isempty(cb_lim)
%         SurfStatColLim([min(nodeStat)-(3*colorbarMod) max(nodeStat)+(3*colorbarMod)]);
        SurfStatColLim([mm-(3*colorbarMod) mx]);
    else
        SurfStatColLim(cb_lim);
    end
    
    set(gcf, 'PaperPositionMode', 'auto');
end

%figout = [ana_dir '/figures/' group{i} '.modules.rh.png'];
%print('-dpng', '-r300', figout); close

