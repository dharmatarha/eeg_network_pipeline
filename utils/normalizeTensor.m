function normalizedTensor = normalizeTensor(inputTensor, stat, verbose)
%% Simple multidimensional array normalization, with sum, mean, or median of all values
% 
% USAGE: normalizedTensor = normalizeTensor(inputTensor, stat='mean', verbose = true)
%
% By normalization we simply mean rescaling by sum / mean / median of all elements.
% Done by linearization of the multidimensional array.
%
% Mandatory input:
% inputTensor        - Numeric tensor to be normalized. 
%
% Optional input:
% stat               - Char array, one of {'sum', 'mean', 'median'}. 
%                      Statistic to normalize (rescale) the values with. 
%                       Defaults to 'mean'. 
%
% Output:
% normalizedTensor   - Numeric tensor, normalized.
%


%% Input checks

if ~ismembertol(nargin, 1:3)
    error('Function "normalizeTensor" requires input arg "inputTensor" while input args "stat" and "verbose" are optional!');
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
if ~isnumeric(inputTensor) || numel(inputTensor)<2
    error('Input arg "inputTensor" shoud be a numeric tensor/matrix/vector');
end


%% Get rescaling value

% linearize
inputLin = inputTensor(:);
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
normalizedTensor = inputTensor/tmp;


return

