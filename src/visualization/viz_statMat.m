function [] = viz_statMat(mat,colorbar_bool,colorbar_axis_vals,colorbar_title,plot_title,midValue)

% take care of possible empty vars
if ~exist('colorbar_bool','var') || isempty(colorbar_bool)
    colorbar_bool = 1;
end
    
if ~exist('colorbar_axis_vals','var') || isempty(colorbar_axis_vals)
    colorbar_axis_vals = 'auto' ;  
end

if ~exist('colorbar_title','var') || isempty(colorbar_title)
    colorbar_title = 'GnBu' ;  
end

if ~exist('plot_title','var') || isempty(plot_title)
    plot_title = '' ;  
end

% add paths and stuff
addpath('figures/')
addpath('figures/DrosteEffect-BrewerMap-04533de/')

% mat = mat;          
bttmInd = tril(ones(size(mat)),-1) ;

% set lower half to NaN so as to not viz
mat(~~bttmInd) = NaN ;

% initial imagesc command
im = imagesc(mat);            %# Create a colored plot of the matrix values

% do not render NaNs
set(im,'alphadata',~isnan(mat))

% colormap('default')
map = brewermap(50,colorbar_title); 
% map = flipud(map) ;
map(1,:) = [0.9 0.9 0.9];                        
colormap(map)                        
                         
textStrings = num2str(mat(:),'%.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
wantInd = mat>0; % only write text for non-zero vals
textStrings(wantInd == 0) = {''} ; % set the 0's to null
[x,y] = meshgrid(1:(size(mat,1)));   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');

if ~exist('midValue','var') || isempty(midValue)
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
end          

%midValue = 0.4 ; % set to make it easier to read
textColors = repmat(mat(:)> midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

% sets the color scale
caxis(colorbar_axis_vals)

if colorbar_bool 
    % get the original size
    colorbar
end
    
if ~isempty(plot_title)
    title(plot_title)
end

set(gca,'XTick',1:(size(mat,1)),...                         %# Change the axes tick mark
        'YTick',1:(size(mat,1)))