function [modules2colors, allColors] = sortModules2Colors(modules, colorNo, varargin)
%% Function to assign colors to modules
% 
% USAGE sortModules2Colors(modules, colorNo, n=10, cMethod='randColor')
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

% if module data is from one layer, or requested "cMethod" is 'nan', we 
% simply set allColors based on modules2colors, meaning NaN for any module 
% above "n" 
if ~multiLayer || strcmp(cMethod, 'nan')   
    for i = 1:moduleNo
        allColors(modules==modules2colors(i, 1)) = modules2colors(i, 2);
    end

% otherwise we go through layers and use the rest of the colors repeatedly
% to mark smaller modules
elseif multiLayer && strcmp(cMethod, 'randColor')
    
    
    
    
end
    

















