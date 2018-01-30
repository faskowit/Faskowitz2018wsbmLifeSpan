function [ ] = viz_statsInMat_wComm(mat,...
    commVec,commColorMapTitle,commColorMapNum,...
    valsColorBarBool,valsColorBarLims,valsColorBarTitle,...
    plot_title,midValue)

% take care of possible empty vars
if ~exist('commColorMapTitle','var') || isempty(commColorMapTitle)
   commColorMapTitle = 'paired' ;
end

if ~exist('valsColorBarBool','var') || isempty(valsColorBarBool)
    valsColorBarBool = 1 ;
end
    
if ~exist('valsColorBarLims','var') || isempty(valsColorBarLims)
    valsColorBarLims = 'auto' ; 
end

if ~exist('valsColorBarTitle','var') || isempty(valsColorBarTitle)
    valsColorBarTitle = 'GnBu'  ;
end

if ~exist('plot_title','var') || isempty(plot_title)
    plot_title = '' ;
end

% mat = mat;          
bttmInd = tril(ones(size(mat)),-1) ;

% set lower half to NaN so as to not viz
mat(~~bttmInd) = NaN ;

% initial imagesc command
im = imagesc(mat);            %# Create a colored plot of the matrix values

% do not render NaNs
set(im,'alphadata',~isnan(mat))

% colormap('default')
cmap = brewermap(50,valsColorBarTitle); 
% map = flipud(map) ;
% cmap(1,:) = [0.9 0.9 0.9];                        
colormap(cmap)                        
                         
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

textColors = repmat(mat(:)> midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

% sets the color scale
caxis(valsColorBarLims)

if valsColorBarBool 
    % get the original size
    colorbar
end
    
if ~isempty(plot_title)
    title(plot_title)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uniqueLabs = unique(commVec);
numComm = length(uniqueLabs);
set(gca,'XTick',1:numComm,...                         %# Change the axes tick mark
        'YTick',1:numComm,...
        'xticklabel',uniqueLabs,...
        'yticklabel',uniqueLabs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%uniqueLabs = unique(commVec);
numComm = length(uniqueLabs);
commCmap = brewermap(commColorMapNum,commColorMapTitle) ;

hold on

for idx = 1:numComm
    ind = find(uniqueLabs == uniqueLabs(idx));
    if ~isempty(ind)
        mn = min(ind) - 0.5;
        mx = max(ind) + 0.5;
        x = [0.5  0.5 ];
        y = [mn mx ];
        plot(x,y,'Color',commCmap(uniqueLabs(idx),:),'linewidth',6); 
        x = [ numComm+0.5 numComm+0.5 ] ;
        plot(y,x,'Color',commCmap(uniqueLabs(idx),:),'linewidth',6);
    end
end
    
hold off  
    
    