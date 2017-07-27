function [] = pictureWSBM( communityAssign, plotWhat , atlas) 
% inputs
% 1) community assigments
% 2) what are they describing? (left, right, both?)
% 3) your only choice right now is yeo

if nargin < 2
    disp('need at least 2 args')
    return 
end

if nargin < 3
    atlas = 'yeo' ;
end

if ~strcmp(atlas,'yeo')
    disp('your only choice is yeo right now')
    return
end

%% if communityAssign is a struct, extract communities
if isstruct(communityAssign)
   
   communityAssign = community_assign(communityAssign.Para.mu);
end

%% continue

if strcmp(plotWhat,'both')
    disp('plot both')
    numParcels = length(communityAssign) / 2; 
    
elseif strcmp(plotWhat,'left')
    disp('plot left')
    numParcels = length(communityAssign) ;
    
elseif strcmp(plotWhat,'right')
    disp('plot right')
    numParcels = length(communityAssign) ;
    
else
    disp('invalid choice')
    return  
end

%% Mapping Networks or Modules
%cmap_spec = ccs_mkcolormap([ccs_dir '/vistool/fimages/spectrum_afni.tif']);

% load the vars I need
% addpath('/home/jfaskowi/JOSHSTUFF/projects/SBM')
load('4_josh_thanks.mat')
load('4_josh_thanks2.mat')

if size(communityAssign,2) == 2
    numModules = max(max(communityAssign(:,2)));
    colIdx = 2 ;
else
    numModules = max(max(communityAssign(:,1)));
    colIdx = 1 ;
end
    
modules_lh = zeros(nVertices_lh,1);
modules_rh = zeros(nVertices_rh,1);
cmap_mod = cmap_spec(round(linspace(1,256,numModules)),:);    
        
for parcelID=1:numParcels

    disp(parcelID)

    if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'left'))
        tmpidx = parcel_vertexid_lh{parcelID};
        modules_lh(tmpidx) = communityAssign(parcelID,colIdx)+1;
    end
    
    if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'right'))
        tmpidx = parcel_vertexid_rh{parcelID};
        modules_rh(tmpidx) = communityAssign((parcelID)+numParcels,colIdx)+1;
    end
    
end
   
disp(plotWhat)

if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'left'))
    
    disp('rendering left')
    
    %render lh surfaces
    modules_lh = modules_lh + medialwall_lh;
    figure('Units', 'pixel', 'Position', [100 100 800 800]); 
    axis off
    SurfStatView(modules_lh, fs_lh, ' ', 'white', 'true'); 
    colormap([0 0 0; 0.5 0.5 0.5; cmap_mod]); 
    SurfStatColLim([-0.5 numModules+1.5]);
    set(gcf, 'PaperPositionMode', 'auto');

end

%figout = [ana_dir '/figures/' group{i} '.modules.lh.png'];
%print('-dpng', '-r300', figout); close

if or(strcmp(plotWhat,'both'),strcmp(plotWhat,'right'))
    
    disp('rendering right')
    
    %render rh surfaces
    modules_rh = modules_rh + medialwall_rh;
    figure('Units', 'pixel', 'Position', [100 100 800 800]); axis off
    SurfStatView(modules_rh, fs_rh, ' ', 'white', 'true'); 
    colormap([0 0 0; 0.5 0.5 0.5; cmap_mod]); 
    SurfStatColLim([-0.5 numModules+1.5]);
    set(gcf, 'PaperPositionMode', 'auto');
end

%figout = [ana_dir '/figures/' group{i} '.modules.rh.png'];
%print('-dpng', '-r300', figout); close

