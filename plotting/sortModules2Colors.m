function [modules2colors, allColors, sortedModFreq] = sortModules2Colors(modules, colorNo, varargin)
%% Function to assign colors to modules
% 
% USAGE [modules2colors, allColors, sortedModFreq] = sortModules2Colors(modules, 
%                                                                       colorNo, 
%                                                                       n=10, 
%                                                                       cMethod='randColor')
%
% If there are enough colors to cover all modules uniquely, assignment is
% simply based on frequencies, so that colors reflect the same frequency 
% relations independent of precise module indices.
%
% If there are more modules than colors, the function tries to assign
% stable colors to the most freuqent "n" modules. For single-layer module 
% data, the rest of the colors are assigned NaN. For multi-layer module 
% data, by default, the function uses the rest of the colors repeatedly
% for the rest of the modules, preferably avoiding temporal overlap of 
% same color assignment for different modules. This behavior can be
% overridden by setting the optional input arg "cMethod" to 'nan', after
% which the function sets all modules above "n" to NaN. 
%
% Works for single- or multi-epoch module assignments as well.
%
% Mandatory inputs:
% modules   - Numeric vector or matrix, containing module indices for
%           ROIs/nodes in either a single epoch (vector) or across multiple
%           epochs (matrix). Size is [node no. X epoch no.]. Values must be
%           positive integers or zero - zeros is treated as a valid module
%           identifier.
% colorNo   - Positive integer value. Number of colors available for
%           assignment to modules.
%
% Optional inputs:
% n         - Positive integer value. Number of most frequent modules being
%           assigned a specific color. Defaults to 10.
% cMethod   - Char array, one of {'randColor', 'nan'}. If 'randColor', and
%           the module data ("modules") is a matrix (multilayer), the 
%           function tries to assign the colors above "n" to rest of the 
%           modules in a way that avoids confusing overlaps, that is, 
%           different modules colored the same way in the same or 
%           consecutive layers. If 'nan', or 'randColor', but module data 
%           for only one layer, every module above "n" is colored as NaN. 
%
% Outputs:
% modules2colors    - Numeric matrix sized [no. of unique modules X 2].
%           First column contains module indices, the second the
%           corresponding color number.
% allColors         - Numeric vector or matrix, same size as input arg
%           "modules". Contains color number for all module indices in
%           "modules".
% sortedModFreq     - Numeric matrix sized [no. of unique modules X 2].
%           First column contains the frequency (sum) of node memberships 
%           for each module, while the second column holds the module 
%           indices. The matrix is sorted by frequencies, that is, the 
%           first row contains the highest frequency value and the 
%           corresponding module index.    
%


%% Input checks

% no. of inputs
if ~ismember(nargin, 2:4)
    error(['Function sortModules2Colors requires input args "modules" and ',...
        '"colorNo" while args "n" and "colorMethod" are optional!']);
end
% mandatory inputs
if ~isnumeric(modules) || length(size(modules)) ~= 2 || any(mod(modules(:), 1)~=0) || any(modules(:)<0)
    error(['Input arg "modules" should be a numeric vector or matrix, ',...
        'containing only positive integers and zero!']);
end
if ~isnumeric(colorNo) || numel(colorNo)~=1 || mod(colorNo, 1)~=0 || colorNo < 0
    error('Input arg "colorNo" should be a positive integer value!');
end
% optional inputs
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && mod(varargin{v}, 1)==0 && varargin{v}>0 && varargin{v}<=colorNo && ~exist('n', 'var')
            n = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'randColor', 'nan'}) && ~exist('cMethod', 'var')
            cMethod = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to "n" or "cMethod"!');
        end  % if
    end  % for v
end  % if ~isempty
% default values
if ~exist('n', 'var')
    n = 10;
end
if ~exist('cMethod', 'var')
    cMethod = 'randColor';
end


%% Get basic info about modules

% single or multi-layer module data?
if ismatrix(modules)
    multiLayer = true;
else
    multiLayer = false;
end

% get array of unique module indices
modIndices = unique(modules(:));
% no. of modules
moduleNo = length(modIndices);
% get frequencies of modules
modIdxFreq = zeros(moduleNo, 1);
for i = 1:moduleNo
    modIdxFreq(i) = sum(modules(:)==modIndices(i));
end
% sort module indices based on frequencies
sortedModFreq = sortrows([modIdxFreq, modIndices], 1, 'descend');

% preallocate output vars
modules2colors = zeros(moduleNo, 2);
allColors = zeros(size(modules));


%% With enough colors, each module is assigned a color based on its freq

if moduleNo <= colorNo
    modules2colors = [sortedModFreq(:, 2), [1:moduleNo]'];
    % for output var "allColors", we just set colors based on
    % modules2colors for each value in "modules"
    for i = 1:moduleNo
        allColors(modules==modules2colors(i, 1)) = modules2colors(i, 2);
    end
    return;
end


%% Otherwise fixed colors are assigned only to first "n" modules
    
% check if we have a sensible value for "n"
if n > colorNo
    error(['Number of top frequent modules to be colored consistently ',...
        '("n") is larger than the number of available colors!']); 
end
if n >= moduleNo
    error(['Number of top frequent modules to be colored consistently ',...
        '("n") is equal to or larger than the number of unique colors!']); 
end    

% for output variable "modules2colors", in all cases (single or multilayer,
% 'randColor' or 'nan' for "cMethod") we assign first "n" colors to most
% frequent "n" modules
modules2colors(:, 1) = sortedModFreq(:, 2);
modules2colors(:, 2) = [[1:n]'; nan(moduleNo-n, 1)];

% basic color assignment is the same for single- and multilayer module
% data, irrespective of "cMethod" - go through rows of modules2colors and
% perform the assignment
for i = 1:moduleNo
    allColors(modules==modules2colors(i, 1)) = modules2colors(i, 2);
end

% if module data is multilayer and the requested coloring method
% ("cMethod") is 'randColor', we o through each layer and assign colors
% remaining after top "n" modules repeatedly to the smaller modules
if multiLayer && strcmp(cMethod, 'randColor')
    % "n" top frequent modules
    topFreqMods = modules2colors(1:n, 1);
    % color indices above "n"
    colorsLeft = n+1:colorNo;
    
    % loop through layers
    for layerIdx = 1:size(modules, 2)
        
        % for first layer we just 
        if layerIdx == 1
            % define modules present in layer "layerIdx" that are not among top frequent ones
            oldModulesLeft = setdiff(modules(:, layerIdx), topFreqMods);
            % assign colors to modules above "n" just based on their order
            oldLayerAssignment = [oldModulesLeft, colorsLeft(1:length(oldModulesLeft))'];
            % get color indices left after layer assignment
            colorsLeftForNext = setdiff(colorsLeft, oldLayerAssignment(:,2));
            % set allColors in layer "layerIdx" appropriately
            for i = 1:length(oldModulesLeft)
                allColors(modules(:, layerIdx)==oldModulesLeft(i), layerIdx) = oldLayerAssignment(i, 2);
            end
            
        % for the rest of the layers
        else

            % define modules present in layer "layerIdx" that are not among top frequent ones
            currentModulesLeft = setdiff(modules(:, layerIdx), topFreqMods);
            
            % go through the modules if there are any besides the top frequent "n"
            if ~isempty(currentModulesLeft)
                % check if we have enough colors left for coloring current
                % layer
                if numel(setdiff(currentModulesLeft, oldModulesLeft)) > numel(colorsLeftForNext)
                    error([char(10), 'We ran out of colors in layer, ', num2str(layerIdx),... 
                        ' when assigning colors to smaller modules. There are only ',... 
                        num2str(numel(colorsLeftForNext)), ' colors but ',... 
                        num2str(numel(setdiff(currentModulesLeft, oldModulesLeft))),... 
                        ' modules to color.']);
                end               
                % preallocate the new layer module-color assignment var
                newLayerAssignment = [currentModulesLeft, zeros(length(currentModulesLeft), 1)];
                % counter for using up the colors left for current layer
                colorCounter = 0;
                % go through each module in current layer not in top "n" 
                for c = 1:length(currentModulesLeft)
                    % if the module was also present in last layer, use the
                    % same color as there
                    if ismember(currentModulesLeft(c), oldModulesLeft)
                        newLayerAssignment(c, 2) = oldLayerAssignment(oldModulesLeft==currentModulesLeft(c), 2);
                    % otherwise assign a new color, not used either by top
                    % "n" modules or in last layer (oldColorsLeft)
                    else
                        colorCounter = colorCounter+1;
                        newLayerAssignment(c, 2) = colorsLeftForNext(colorCounter);
                    end  % if ismember
                end  % for c
                
                % set allColors in layer "layerIdx" appropriately
                for i = 1:length(currentModulesLeft)
                    allColors(modules(:, layerIdx)==currentModulesLeft(i), layerIdx) = newLayerAssignment(i, 2);
                end
                
            % if there are only the top "n" modules present in layer, we 
            % still define "newLayerAssignment"       
            else
                newLayerAssignment = modules2colors(1:n, 1:2);
                
            end  % if ~isempty
            
            % define "old" vars for the assignments in next layer
            oldModulesLeft = currentModulesLeft;
            oldLayerAssignment = newLayerAssignment;
            colorsLeftForNext = setdiff(colorsLeft, oldLayerAssignment(:,2));
            
        end  % if layerIdx == 1
        
    end  % for layerIdx
     
end  % if multiLayer &&
    


return
















