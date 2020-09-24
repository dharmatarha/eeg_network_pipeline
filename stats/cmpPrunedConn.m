function [permRes, withinCondPermRes, connSim] = cmpPrunedConn(connData, varargin)
%% Compare pruned connectivity matrices across conditions
%
% USAGE: [permRes, withinCondPermRes, connSim] = cmpPrunedConn(connData, metric='corr', permNo=10000, permStat='mean', maskThr=0)
%
% Compares connectivity patterns across different groups of epoch-pairings 
% using permutation tests. Very similar to cmpFullConn, but works on 
% pruned connectivity data, that is, it expects NaN values in the upper
% triangles of connectivity matrices!
%
% Input data is expected to contain undirected connectivity values. From
% each connectivity matrix we only take the upper triangle into account.
%
% NaN values render the connectivity matrices different from
% epoch-to-epoch, posing a problem for similarity measures. We deal with
% this in two steps:
% (1) Only the edges present in more than maskThr ratio of epochs is taken
%       into account. All other edges are set to zero (not NaN anymore!).
%       By default this step is omitted (maskThr=0), as its results are
%       hard to interpret.
% (2) When calculating similarity, we delete the elements where both 
%       vectors (containing connectivity values) are zero. So, if, after 
%       step (1), one pruned connectivity matrix has an edge between two 
%       ROIs but the other matrix does not, we take those values into
%       account. However, when both values are zero we ignore (delete) 
%       them, to avoid inflating similarity due to missing edges.
%
% For comparisons, first, for each condition, the function tests the
% within-condition epoch-pairings against the across-condition epoch
% pairings. Second, for each condition-pairing, it tests if epoch-pairings
% in one condition have different values than in the other condition 
% (as in the case of pairwise comparisons). 
%
% For the tests we rely on permTest.m, input args "permNo" and "permStat"
% (if specified) are passed to permTest.m.
%
% Mandatory input:
% connData      - 4D numeric array, sets of connectivity matrices across
%               epochs and conditions (stimuli). First two dimensions have
%               equal size and determine a connectivity matrix, third
%               dimension is epochs, fourth is conditions (stimuli).
%
% Optional inputs:
% metric        - String specifying distance metric for connectivity
%               matrix comparisons. Matrices are first vectorized, 
%               then one of these distances is used: {'corr', 'eucl'}.  
% permNo        - Numeric value, the number of permutations to perform for
%               random permutation tests. One of 100:100:10^6.
% permStat      - String specifying the statistic we perform the random
%               permutaiton tests on. One of {'mean', 'median', 'std'} -
%               the ones supported by permTest.m
% maskThr       - Numeric value, threshold for frequency-based masking of
%               edges. Edges that are not present in more than maskThr ratio
%               of all epochs are set to 0. Defaults to 0, must be one 
%               of 0:0.01:1.
%
% Outputs:
% permRes       - Struct, random permutation test results. Each element of
%               permRes summerizes the results of the within-condition 
%               epoch pairings versus across-condition epoch-pairings test 
%               for a condition, i.e. its size is [1 "numberOfConditions"].
%               Has fields for the mean, median and SD values of both
%               epoch-pairing groups, and for the outcomes of permTest,
%               i.e. estimated p-value, real test stat difference and
%               permuted difference values.
% withinCondPermRes       - Struct, random permutation test results. 
%               Each element of withinCondPermRes summerizes the results of
%               a comparison of epoch-pairings across two conditions, i.e. 
%               its size is [1 nchoosek("numberOfConditions", 2)].
%               Has fields for the mean, median and SD values of both
%               epoch-pairing groups, and for the outcomes of permTest,
%               i.e. estimated p-value, real test stat difference and
%               permuted difference values.
% connSim       - Numeric matrix containing connectivity similarity values
%               for all epoch-pairings. Symmetric, as the upper triangle
%               values are calculated and mirrored to the lower half. Its
%               size is ["numbreOfEpochs"*"numberOfConditions",
%               "numbreOfEpochs"*"numberOfConditions"].
%


%% Input checks

% number of args
if ~ismember(nargin, 1:5)
    error('Wrong number of input args - "connData" is needed while "metric", "permNo", "permStat" and "maskThr" are optional!');
end
% loop through varargin to sort out input args
if nargin > 1
    for v = 1:length(varargin)
        if ischar(varargin{v})
            if ismember(varargin{v}, {'corr', 'eucl'})
                metric = varargin{v};
            elseif ismember(varargin{v}, {'mean', 'median', 'std'})
                permStat = varargin{v};
            else
                error('An input arg (string) could not be mapped to any optional arg!');
            end
        elseif isnumeric(varargin{v})
            if ismember(varargin{v}, 100:100:10^6)
                permNo = varargin{v};
            elseif ismember(varargin{v}, 0:0.01:1)
                maskThr = varargin{v};
            end
        else
            error('At least one input arg could not mapped to any optional arg!');
        end
    end
end    
% assign default values where necessary
if ~exist('metric', 'var')
    metric = 'corr';
end
if ~exist('permStat', 'var')
    permStat = 'mean';
end
if ~exist('permNo', 'var')
    permNo = 10^4;
end  
if ~exist('maskThr', 'var')
    maskThr = 0;
end  
% check size and dimensionality of mandatory arg "connData"
if ~isequal(size(connData, 1), size(connData, 2))
    error('First two dimensions of input arg "connData" need to have equal size!');
end
if length(size(connData)) ~= 4
    error('Input arg "connData" should have four dimensions!');
end

% user message
disp([char(10), 'Called cmpPrunedConn function with input args:',...
    char(10), 'Input data: array with size ', num2str(size(connData)),...
    char(10), 'Metric: ', metric,...
    char(10), 'No. of random permutations: ', num2str(permNo),...
    char(10), 'Test statistic for random permutations: ', permStat,...
    char(10), 'Edge masking threshold: ', num2str(maskThr)]);


%% Basics

% number of ROIS, conditions and epochs
roiNo = size(connData, 1);
condNo = size(connData, 4);
epochNo = size(connData, 3);
% preallocate for linearized data
dataLin = zeros((roiNo*roiNo-roiNo)/2, epochNo, condNo);


%% Create the mask of edges that are present in (maskThr*100)% of epochs

connMask = ~isnan(connData);  % mask of edges that survived pruning
connMaskSum = sum(sum(connMask, 4), 3);  % sum across all epochs
connMaskFinal = connMaskSum > round(maskThr*epochNo*condNo);  % mask of edges that are present in more than maskThr ratio of epochs
edgeNo = sum(sum(connMaskFinal));  % number of remaining edges

% set everything outside the mask to NaN, and all NaN values to zero
maskTmp = repmat(connMaskFinal, [1, 1, 40, 4]);  % 4-D repetitions of mask
connData(~maskTmp) = NaN;
connData(isnan(connData)) = 0;

% user message
disp([char(10), 'Only considering edges that are present in more than ',... 
    num2str(round(maskThr*epochNo*condNo)), ' epochs (maskThr=',... 
    num2str(maskThr), '), remaining no. of edges is ', num2str(edgeNo)]);


%% Linearize and reshape connectivity data

for cond = 1:condNo
    for epoch = 1:epochNo
        % get upper triangular part from epoch-level connectivity matrix in
        % a vector
        tmp = squeeze(connData(:, :, epoch, cond))';
        idx = tril(true(size(tmp)), -1);  % do not include data from main diagonal
        dataLin(:, epoch, cond) = tmp(idx)';      
    end
end

% conditions become epochs one after another, so e.g. epoch 2 at cond 3
% becomes epoch 22
dataLin = reshape(dataLin, [(roiNo*roiNo-roiNo)/2, epochNo*condNo]);

% delete rows where all values were set to zero (edges not present in
% enough epochs)
dataLin(~any(dataLin, 2), :) = [];

% user message
disp([char(10), 'Vectorized connectivity matrices (upper triangles)']);


%% Get similarity across all epoch-pairings

% preallocate
connSim = nan(epochNo*condNo);  % connectivity pattern similarity across epochs
pairCounter = 0;  % counter for pairings calculated
pastPairings = zeros((epochNo^2*condNo^2-epochNo*condNo)/2, 2);  % var storing epoch pairings calculated
% loops through all epochs
for epochOne = 1:epochNo*condNo
    for epochTwo = 1:epochNo*condNo
        % only consider epoch-pairing if the two epochs are not the same
        % and have not been encountered before
        if epochOne ~= epochTwo && ~ismember([epochOne, epochTwo], pastPairings, 'rows') && ~ismember([epochTwo, epochOne], pastPairings, 'rows')
            % store epoch-pairing
            pairCounter = pairCounter+1;
            pastPairings(pairCounter, :) = [epochOne, epochTwo];
            % clear cases where an edge is missing in both vectors, i.e.
            % both are zero - such cases could seriously inflate apparent proximity
            dataOne = dataLin(:, epochOne); 
            dataTwo = dataLin(:, epochTwo);
            zeroIdx = dataOne==0 & dataTwo==0;
            dataOne(zeroIdx) = [];
            dataTwo(zeroIdx) = [];           
            % calculate similarity according to the supplied metric
            switch metric
                case 'corr'
                    connSim(epochOne, epochTwo) = corr(dataOne, dataTwo);
                case 'eucl'
                    connSim(epochOne, epochTwo) = norm(dataOne-dataTwo);
            end
        end
    end
end

% user message
disp([char(10), 'Calculated similarity across all epoch-pairings']);


%% Compare similarity values between within- and across-condition pairings

% first mirror the upper triangular part of connSim matrix to the lower
% half
connSim = triu(connSim, 1) + triu(connSim, 1)'; 

% preallocate results struct
permRes = struct;
permRes.withinCondMean = nan;
permRes.withinCondSD = nan;
permRes.withinCondMedian = nan;
permRes.acrossCondMean = nan;
permRes.acrossCondSD = nan;
permRes.acrossCondMedian = nan;
permRes.realDiff = nan;
permRes.permDiff = nan(permNo, 1);
permRes.pEst = nan;

for cond = 1:condNo
    
    % within-cond epoch pairing similarities for given condition
    tmp = connSim((cond-1)*epochNo+1:cond*epochNo, (cond-1)*epochNo+1:cond*epochNo);
    tmp = triu(tmp, 1);  
    % linearize to vector
    tmpT = tmp';
    idx = tril(true(size(tmpT)), -1);
    withinCondSim = tmpT(idx)';
    
    % across-condition pairings for epochs in given condition
    tmpParts = cell(2, 1);
    if cond ~= 1
        tmpParts{1} = connSim((cond-1)*epochNo+1:cond*epochNo, 1:(cond-1)*epochNo);
    end
    if cond ~= condNo
        tmpParts{2} = connSim((cond-1)*epochNo+1:cond*epochNo, cond*epochNo+1:end);
    end
    % concatenate the parts
    if isempty(tmpParts{1})
        tmpAll = tmpParts{2};
    elseif isempty(tmpParts{2})
        tmpAll = tmpParts{1};
    else
        tmpAll = [tmpParts{1}, tmpParts{2}];
    end
    % linearize into a vector
    acrossCondSim = reshape(tmpAll', [1, epochNo^2*(condNo-1)]);
    
    % store descriptives in result struct
    permRes(cond).withinCondMean = mean(withinCondSim);
    permRes(cond).withinCondSD = std(withinCondSim);
    permRes(cond).withinCondMedian = median(withinCondSim);
    permRes(cond).acrossCondMean = mean(acrossCondSim);
    permRes(cond).acrossCondSD = std(acrossCondSim); 
    permRes(cond).acrossCondMedian = median(acrossCondSim);
    
    % permutation test across the two groups of connectivity similarity
    % values
    [permRes(cond).pEst, permRes(cond).realDiff,... 
        permRes(cond).permDiff] = permTest(withinCondSim, acrossCondSim, permNo, permStat);
    
end

% user message
disp([char(10), 'Compared similarity within- versus across-condition epoch-pairings']);


%% Compare similarity values between different within-condition pairings
% e.g. between epoch-pairings in within cond 1 and epoch-pairings within
% cond 2

% preallocate results struct
withinCondPermRes = struct;
withinCondPermRes.condOneMean = nan;
withinCondPermRes.condOneSD = nan;
withinCondPermRes.condOneMedian= nan;
withinCondPermRes.condTwoMean = nan;
withinCondPermRes.condTwoSD = nan;
withinCondPermRes.condTwoMedian = nan;
withinCondPermRes.realDiff = nan;
withinCondPermRes.permDiff = nan(permNo, 1);
withinCondPermRes.pEst = nan;

% counter for condition pairings
pairCounter = 0;
% var storing condition pairings already calculated
pastPairings = zeros(nchoosek(4, 2), 2);

% loops through conditions
for condOne = 1:condNo
    for condTwo = 1:condNo
        % only calculate if the two conditions are not the same and have
        % not been take into account yet
        if condOne~=condTwo && ~ismember([condOne, condTwo], pastPairings, 'rows') && ~ismember([condTwo, condOne], pastPairings, 'rows')
            pairCounter = pairCounter+1;
            pastPairings(pairCounter, :) = [condOne, condTwo];
            
            % within-cond epoch pairing similarities for condition one
            tmp = connSim((condOne-1)*epochNo+1:condOne*epochNo, (condOne-1)*epochNo+1:condOne*epochNo);
            tmp = triu(tmp, 1);  
            % linearize to vector
            tmpT = tmp';
            idx = tril(true(size(tmpT)), -1);
            withinCondSimOne = tmpT(idx)';

            % within-cond epoch pairing similarities for condition two
            tmp = connSim((condTwo-1)*epochNo+1:condTwo*epochNo, (condTwo-1)*epochNo+1:condTwo*epochNo);
            tmp = triu(tmp, 1);  
            % linearize to vector
            tmpT = tmp';
            idx = tril(true(size(tmpT)), -1);
            withinCondSimTwo = tmpT(idx)';
            
            % store descriptives in results struct
            withinCondPermRes(pairCounter).condOneMean = mean(withinCondSimOne);
            withinCondPermRes(pairCounter).condOneSD = std(withinCondSimOne);
            withinCondPermRes(pairCounter).condOneMedian = median(withinCondSimOne);
            withinCondPermRes(pairCounter).condTwoMean = mean(withinCondSimTwo);
            withinCondPermRes(pairCounter).condTwoSD = std(withinCondSimTwo);
            withinCondPermRes(pairCounter).condTwoMedian = median(withinCondSimTwo);
            
            % permutation test across the two groups of connectivity similarity
            % values
            [withinCondPermRes(pairCounter).pEst, withinCondPermRes(pairCounter).realDiff,... 
                withinCondPermRes(pairCounter).permDiff] = permTest(withinCondSimOne, withinCondSimTwo, permNo, permStat);
            
        end
        
    end
end

% user message
disp([char(10), 'Compared similarity across conditions for within-condition epoch-pairings']);
disp('Done with everything!');

return



