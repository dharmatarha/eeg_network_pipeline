function [normalizedMatrix] = normalizeMatrix(inputMatrix, stat, verbose)
%% Matrix normalization, with sum, mean, or median of all values
% 
% USAGE: normalizedMatrix = normalizeMatrix(inputMatrix, stat='mean', verbose = true)
%
% By normalization we simply mean rescaling by sum / mean / median of all elements.
%
% Mandatory input:
% inputMatrix       - The matrix to be normalized
%
% Optional input:
% stat              - Char array, one of {'sum', 'mean', 'median'}. 
%                     Statistic to normalize (rescale) the values with. 
%                     Defaults to 'mean'.
%
% Output:
% normalizedMatrix  - Normalized matrix
%


%% Input checks

if ~ismembertol(nargin, 1:3)
    error('Function "normalizeMatrix" requires input arg "inputMatrix" while input args "stat" and "verbose" are optional!');
end
if nargin == 1
    stat = 'mean';
    verbose = true;
elseif nargin == 2
    stat = 'mean';
    if ~islogical(verbose) || numel(verbose)~=1
        error('Optional input arg "verbose" must be a logical value!');
    end
else    
    if ~ismember(stat, {'sum', 'mean', 'median'})
        error('Optional input arg "stat" must be one of {''sum'', ''mean'', ''median''}!');
    end
    if ~islogical(verbose) || numel(verbose)~=1
        error('Optional input arg "verbose" must be a logical value!');
    end
end
if ~isnumeric(inputMatrix) || ~ismatrix(inputMatrix)
    error('Input arg "inputMatrix" shoud be a numeric matrix');
end


%% Get rescaling value

% linearize
inputLin = inputMatrix(:);
% check for NaN
if any(isnan(inputLin)) && verbose
    warning('There were NaN values in the input. We use ''omitnan'' flags, so the result will be correct.');
end
% rescaling constant
switch stat
    case 'sum'
        tmp = sum(inputLin, 'omitnan');
    case 'mean'
        tmp = mean(inputLin, 'omitnan');
    case 'median'
        tmp = median(inputLin, 'omitnan');
end
% normalize / rescale
normalizedMatrix = inputMatrix/tmp;


return

