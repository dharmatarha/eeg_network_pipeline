function cmpFullConn(connData, metric, perm)
%% Compare full connectivity matrices across conditions
%
% USAGE: cmpFullConn(connData, metric='corr', surrNo=10000)
%
% Compares connectivity patterns across different conditions using
% permutation tests. 
%
% Input data is expected to contain undirected connectivity values. From
% each connectivity matrix we only take the upper triangle into account.
%
% Inputs:
% connData      - Numerical tensor, sets of connectivity matrices across
%               epochs and conditions (stimuli). First two dimensions have
%               equal size and determine a connectivity matrix, third
%               dimension is epochs, fourth is conditions (stimuli).
% metric        - String specifying distance metric for connectivity
%               matrix comparisons. Matrices are first vectorized, 
%               then one of these ditances is used: {'corr', 'eucl'}.  
%

%% Input checks

% number of args
if ~ismember(nargin, [1 2])
    error('Wrong number of input args - "connData" is needed while "metric" is optional!');
end
% assign default value to "metric" if no value was provided
if nargin == 1
    metric = 'corr';
else
    % check value provided for "metric"
    if ~ismember(metric, {'corr', 'eucl'})
        error('Input arg "metric" must be one of {"corr", "eucl"}!');
    end
end
% check dimensionality of arg "connData"
if ~isequal(size(connData, 1), size(connData, 2))
    error('First two dimensions of input arg "connData" need to have equal size!');
end
if length(size(connData)) ~= 4
    error('Input arg "connData" should have four dimensions!');
end

% user message
disp([char(10), 'Called cmpFullConn function with input args:',...
    char(10), 'Metric: ', metric,...
    char(10), 'Input data is array with size ', num2str(size(connData))]);


%% Basics

% number of ROIS, conditions and epochs
roiNo = size(connData, 1);
condNo = size(connData, 4);
epochNo = size(connData, 3);
% preallocate for linearized data
dataLin = zeros((roiNo*roiNo-roiNo)/2, epochNo, condNo);


%% Linearize and reshape connectivity data

for cond = 1:condNo
    for epoch = 1:epochNo
        % get upper triangular part from epoch-level connectivity matrix in
        % a vector
        tmp = squeeze(connData(:, :, epoch, cond))';
        idx = tril(true(size(tmp)), -1);  % do not include data from main diagonal
        dataLin(:, epoch, cond) = tmp(idx)';
        % warning if there are NaN values
        if any(isnan(dataLin(:, epoch, cond)))
            warning(['There is at least one NaN value among the ',...
                'connectivity values at condition ', num2str(cond),... 
                ', epoch ', num2str(epoch), '!']);
        end
    end
end

% conditions become epochs one after another, so e.g. epoch 2 at cond 3
% becomes epoch 22
dataLin = reshape(dataLin, [(roiNo*roiNo-roiNo)/2, epochNo*condNo]);


%% Get similarity across all epoch-pairings

% preallocate
connSim = nan(epochNo*epochNo);  % connectivity pattern similarity across epochs
pairCounter = 0;  % counter for pairings calculated
pastPairings = zeros((epochNo*epochNo-epochNo)/2, 2);  % var storing epoch pairings calculated
% loops through all epochs
for epochOne = 1:epochNo
    for epochTwo = 1:epochNo
        % only consider epoch-pairing if the two epochs are not the same
        % and have not been encountered before
        if epochOne ~= epochTwo && ~ismember([epochOne, epochTwo], pastPairings, 'rows') && ~ismember([epochTwo, epochOne], pastPairings, 'rows')
            % store epoch-pairing
            pairCounter = pairCounter+1;
            pastPairings(pairCounter, :) = [epochOne, epochTwo];
            % calculate similarity according to the supplied metric
            switch metric
                case 'corr'
                    connSim(epochOne, epochTwo) = corr(dataLin(:, epochOne), dataLin(:, epochTwo));
                case 'eucl'
                    connSim(epochOne, epochTwo) = norm(dataLin(:, epochOne)-dataLin(:, epochTwo));
            end
        end
    end
end


%% Compare similarity values between within- and across-condition pairings

% first mirror the upper triangular part of connSim matrix to the lower
% half
connSim = triu(connSim, 1) + triu(connSim, 1)'; 

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
    acrossConSim = reshape(tmpAll', [1, epochNo*(condNo-1)]);
    
    % permutation test across the two groups of connectivity similarity
    % values
    [pEst, realDiff, permDiff] = permTest(withinCondSim, acrossConSim, perm);
    
    % collect results into aggregate matrices
    
end





