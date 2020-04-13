function prunedMask = getPrunedMask(prunedData, outputMode, thrs)
%% Calculates a mask of edges for a set of pruned connectivity matrices
%
% USAGE: prunedMask = getPrunedMask(prunedData, mode, thrs)
%
% A set of pruned connectivity matrices - even if from the same subject and
% task - is usually composed of partially overlapping elements. This 
% function derives a common mask based on the frequency of occurance of
% edges across the whole set.
% As our use cases are EEG datasets with size [ROIs, ROIs, epochs,
% conditions], inputs are either 4D or 5D arrays. In the latter case, the
% 5th dimension is for subjects. Masks can be either 2D (size [ROIs, ROIs])
% or 4D (size [ROIs, ROIs, epochs, conditions], only in the case of a 
% 5D input array). Arg "outputMode" controls which type of output is returned.
% Edge selection in the mask is controlled by simple thresholds (arg
% "thrs") - edges present in more than "thrs" ratio of matrices are kept in
% the mask (see details below).
%
% Important: Missing edges are expected to be marked with NaN.
%
% Expects undirected connectivity values and works only with the upper
% triangle of connectivity matrices.
%
% Inputs:
% prunedData    - 4D or 5D numeric array with the first two dims defining a
%           connectivity matrix for a given epoch, condition (and subject
%           in case of 5D). Pruned connectivity, meaning that the function
%           expects NaN values where edges were not significant.
% outputMode    - String, one of {'2D', '4D'}. Provides the dimensionality
%           of the desired output, i.e. edge mask. If the input is 4D, it
%           can only be '2D'.
% thrs          - Numeric vector with one or two values. 
%           If "prunedData" is 4D, "thrs" can only be one value and the 2D 
%           output "prunedMask" will have ones where there are edges in 
%           more than "thrs" ratio of all epochs in "prunedData". 
%           If "prunedData" is 5D and "outputMode" is '4D', "thrs" can
%           again be only one value and the returned mask is one where
%           there were edges in more than "thrs" ratio of all subjects (5th
%           dim of input data).
%           If "prunedData" is 5D and "outputMode" is '2D', "thrs" can be
%           two values, the first for selecting edges across subjects and
%           the second for selecting edges across all epochs (calculated in
%           two successive steps). If "thrs" is one value, it is treated 
%           as an overal threshold and edges present in more than "thrs" 
%           ratio of all epochs across all subjects are selected. One of
%           0:00.1:0.99.
% 
% Output:
% prunedMask    - Numeric array, either 2D or 4D, depending on
%           "outputMode". Ones for edges passing the threshold and zero
%           otherwise.
%


%% Input checks
if nargin ~= 3
    error('Function getPrunedMask requires three args: "prunedData", "outputMode" and "thrs"!');
end
if ~ismember(length(size(prunedData)), [4 5])
    error('Input arg "prunedData" is expected to be a 4D or 5D array!');
end
if ~isequal(size(prunedData, 1), size(prunedData, 2))
    error(['First two dimensions of input arg "prunedData" are expected to ',...
        'be equal (square matrices with connectivity values)!']);
end
if ~ismember(outputMode, {'2D', '4D'})
    error('Input arg "outputMode" must be one of "2D" or "4D"!');
end
for i = 1:numel(thrs)
    if ~ismember(thrs(i),0:0.01:0.99)
        error('Elements of input arg "thrs" must be in range 0:0.01:0.99!');
    end
end
if ~ismember(numel(thrs), [1 2])
    error('Input arg "thrs" is expected to contain 1 or 2 values!');
end
if length(size(prunedData))==4 && strcmp(outputMode, '4D')
    warning(['Input arg "prunedData" is a 4D array, so output "prunedMask"',...
        'can only be 2D! Setting "outputMode" to "2D".']);
    outputMode = '2D';
end
if length(size(prunedData))==4 && numel(thrs)==2
    error('Input arg "prunedData" is a 4D array, arg "thrs" can only contain one value!');
end
if length(size(prunedData))==5 && strcmp(outputMode, '4D') && numel(thrs)==2
    error('Input arg "prunedData" is a 5D array and "outputMode" is set to "4D", so arg "thrs" can only contain one value!');
end

% user message
disp([char(10), 'Called getPrunedMask with inputs: ',...
    char(10), 'Input data (pruned connectivity array) with size: ', num2str(size(prunedData)),...
    char(10), 'Output mask dimensionality: ', outputMode, ...
    char(10), 'Threshold(s) to use: ', num2str(thrs)])


%% Go through the cases, derive the mask

edgeFlag = 0;

% name input data dimensions
roiNo = size(prunedData, 1);
epochNo = size(prunedData, 3);
condNo = size(prunedData, 4);
if length(size(prunedData)) == 5
    subNo = size(prunedData, 5);
end


% 4D input, 2D output, one "thrs" value
if length(size(prunedData)) == 4
    
    % derive mask of edges present in more than thrs ratio of all epochs
    % (across all conditions)
    maskTmp = ~isnan(prunedData);  % mask of edges that survived pruning
    maskTmpSum = sum(sum(maskTmp, 4), 3);  % sum across all epochs
    prunedMask = maskTmpSum > round(thrs*epochNo*condNo);  % mask of edges that are present in more than thrs ratio of epochs
    edgeNo = sum(sum(prunedMask));  % number of remaining edges
%     % set everything outside the mask to NaN, and all NaN values to zero
%     maskTmp = repmat(connMaskFinal, [1, 1, 40, 4]);  % 4-D repetitions of mask
%     connData(~maskTmp) = NaN;
%     connData(isnan(connData)) = 0;


% 5D input, 2D output, one "thrs" value
elseif length(size(prunedData)) == 5 && strcmp(outputMode, '2D') && numel(thrs)==1
    
    % derive mask of edges present in more than thrs ratio of all epochs
    % (across all conditions and subjects)
    maskTmp = ~isnan(prunedData);  % mask of edges that survived pruning
    maskTmpSum = sum(sum(sum(maskTmp, 5), 4), 3);  % sum across all epochs
    prunedMask = maskTmpSum > round(thrs*epochNo*condNo*subNo);  % mask of edges that are present in more than thrs ratio of all epochs
    edgeNo = sum(sum(prunedMask));  % number of remaining edges    

% 5D input, 2D output, two "thrs" values
elseif length(size(prunedData)) == 5 && strcmp(outputMode, '2D') && numel(thrs)==2

    % First step:
    % derive mask of edges present in more than thrs(1) ratio of all subjects
    maskTmp = ~isnan(prunedData);  % mask of edges that survived pruning
    maskTmpSum = sum(maskTmp, 5);  % sum across all epochs
    prunedMask = maskTmpSum > round(thrs(1)*subNo);  % mask of edges that are present in more than thrs ratio of subjects
    % generate masked data 
    maskTmp = repmat(prunedMask, [1, 1, 1, 1, subNo]);  % 5D repetitions of mask
    prunedData(~maskTmp) = NaN;
    % Second step:
    % derive mask of edges present in more than thrs(2) ratio of all epochs
    % (across all conditions and subjects)
    maskTmp = ~isnan(prunedData);  % mask of edges that survived pruning
    maskTmpSum = sum(sum(sum(maskTmp, 5), 4), 3);  % sum across all epochs
    prunedMask = maskTmpSum > round(thrs(2)*epochNo*condNo*subNo);  % mask of edges that are present in more than thrs ratio of all epochs
    edgeNo = sum(sum(prunedMask));  % number of remaining edges     
    
    
% 5D input, 4D output, one "thrs" value
elseif length(size(prunedData)) == 5 && strcmp(outputMode, '4D') && numel(thrs)==1
    
    % derive mask of edges present in more than thrs ratio of all subjects
    maskTmp = ~isnan(prunedData);  % mask of edges that survived pruning
    maskTmpSum = sum(maskTmp, 5);  % sum across all epochs
    prunedMask = maskTmpSum > round(thrs*subNo);  % mask of edges that are present in more than thrs ratio of subjects
    edgeNo = squeeze(sum(sum(prunedMask)));  % number of remaining edges per epoch
    edgeNo = reshape(edgeNo, [epochNo*condNo, 1]);
    edgeFlag = 1;
    
end
    

%% Display information on number of edges, return

% message on edges
if edgeFlag
    m = ['There are ', num2str(mean(edgeNo)),... 
        ' (SD=', num2str(std(edgeNo)),') edges per epoch on average in the mask.'];
else
    m = ['There are ', num2str(edgeNo), ' edges in the mask.'];
end

disp([char(10), 'Done, returning mask. ',... 
    char(10), m]);


return










