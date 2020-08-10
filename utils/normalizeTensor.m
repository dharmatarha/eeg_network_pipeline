function normalizedTensor = normalizeTensor(inputTensor, stat, verbose)
%% Simple multidimensional array normalization, either with sum or mean of all values
% 
% USAGE: normalizedTensor = normalizeTensor(inputTensor, stat='mean', verbose = true)
%
% By normalization we simply mean rescaling by sum / mean of all elements.
% Done iteratively for all dimensions.
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


%% Check for NaN values

tmp = squeeze(any(isnan(inputTensor)));
% iterate while there are dimensions to squeeze
while numel(tmp)~=1
    tmp = squeeze(any(tmp));
end
if tmp && verbose
    warning('There were NaN values in the input. We use ''omitnan'' flags, so the result will be correct.');
end
    

%% Normalization 

% get first estimate for value to rescale with
switch stat
    case 'sum'
        tmp = squeeze(sum(inputTensor, 'omitnan'));
    case 'mean'
        tmp = squeeze(mean(inputTensor, 'omitnan'));
    case 'median'
        tmp = squeeze(median(inputTensor, 'omitnan'));
end
% iterate while there are dimensions to squeeze
while numel(tmp)~=1
    switch stat
        case 'sum'
            tmp = squeeze(sum(tmp, 'omitnan'));
        case 'mean'
            tmp = squeeze(mean(tmp, 'omitnan'));
        case 'median'
            tmp = squeeze(median(tmp, 'omitnan'));
    end
end
% rescaling
normalizedTensor = inputTensor/tmp;


return

