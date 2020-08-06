function normalizedTensor = normalizeTensor(inputTensor, stat)
%% Simple multidimensional array normalization, either with sum or mean of all values
% 
% USAGE: normalizedTensor = normalizeTensor(inputTensor, stat='mean')
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

if ~ismembertol(nargin, 1:2)
    error('Function "normalizeTensor" requires input arg "inputTensor" while input arg "stat" is optional!');
end
if nargin == 1
    stat = 'mean';
else
    if ~ismember(stat, {'sum', 'mean', 'median'})
        error('Optional input arg "stat" must be one of {''sum'', ''mean'', ''median''}!');
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
if tmp
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

