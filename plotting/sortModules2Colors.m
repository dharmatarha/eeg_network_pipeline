function sortModules2Colors(modules, colorNo, method)
%% Function to assign colors to modules
% 
% USAGE sortModules2Colors(modules, colorNo, method)
%
% If there are enough colors to cover all modules uniquely, assignment is
% simply based on frequencies, so that colors reflect the same frequency 
% relations independent of precise module indices.
%
% If there are more modules than colors, the function tries to assign
% stable colors to the most freuqent ten modules, and use recurring colors
% for the rest of the modules, preferably avoiding temporal overlap of same
% color assignment for different modules.
%
% Works for single- or multi-epoch module assignments as well.
%
% Inputs:
% modules   - Numeric vector or matrix, containing module indices for
%           ROIs/nodes in either a single epoch (vector) or across multiple
%           epochs (matrix). Size is [node no. X epoch no.]. Values must be
%           positive integers or zero - zeros is treated as a valid module
%           identifier.
% colorNo   - Positive integer value. Number of colors available for
%           assignment to modules.
%
%


%% Input checks

% no. of inputs
if nargin ~= 2
    error('Function sortModules2Colors requires input args "modules" and "colorNo"!');
end
% mandatory inputs
if ~isnumeric(modules) || length(size(modules)) ~= 2 || any(mod(modules(:), 1)~=0) || any(modules(:)<0)
    error(['Input arg "modules" should be a numeric vector or matrix, ',...
        'containing only positive integers and zero!']);
end
if ~isnumeric(colorNo) || numel(colorNo)~=1 || 


% get frequencies of modules
modLin = modules(:);
modIndices = unique(modLin);
modIdxFreq = zeros(size(modIndices));
for i = 1:length(modIndices)
    modIdxFreq(i) = sum(modLin==modIndices(i));
end