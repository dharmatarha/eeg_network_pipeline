function [permRes, connSim] = cmpFullConn_edge(connData, varargin)
%% Compare full connectivity matrices across conditions
%
% USAGE: [permRes, connSim] = cmpFullConn_edge(connData, metric='corr', permNo=10000, permStat='mean')
%
% Evaluates the contribution of each edge to the difference between 
% within- vs. across-condition similarities calculated on the full 
% connectivity matrix.
% Compares connectivity (edge) values across different groups of 
% epoch-pairings using permutation tests. For each condition, it tests the
% within-condition epoch-pairings against the across-condition epoch
% pairings. 
%
% Input data is expected to contain undirected connectivity values. From
% each connectivity matrix we only take the upper triangle into account.
%
% For the tests we rely on permTest.m, input args "permNo" and "permStat"
% (if specified) are passed to permTest.m.
%
% Mandatory input:
% connData      - 4D Numerical tensor, sets of connectivity matrices across
%               epochs and conditions (stimuli). First two dimensions have
%               equal size and determine a connectivity matrix, third
%               dimension is epochs, fourth is conditions (stimuli).
%
% Optional inputs:
% metric        - String specifying distance metric for connectivity
%               matrix comparisons. Matrices are first vectorized, 
%               then one of these distances is used: {'corr', 'eucl'}. 
%               Defaults to 'corr'. 
% permNo        - Numeric value, the number of permutations to perform for
%               random permutation tests. One of 100:100:10^6. Defaults to
%               10^4.
% permStat      - String specifying the statistic we perform the random
%               permutation tests on. One of {'mean', 'median', 'std'} -
%               the ones supported by permTest.m. Defaults to 'mean'.
%
% Outputs:
% permRes       - Struct, random permutation test results. Its fields 
%               summerize the results of the within-condition 
%               epoch pairings versus across-condition epoch-pairings
%               comparison.
%               Has fields for the mean, median and SD values of both
%               epoch-pairing groups, and for the outcomes of permTest,
%               i.e. estimated p-value, real test stat difference, mean and
%               SD of
%               permuted difference values.
% connSim       - 3D Numeric array containing connectivity similarity values
%               for all edges and epoch-pairings. Each layer (as defined 
%               by the 2nd and 3rd dimensions is symmetric, as the upper 
%               triangle values are calculated and mirrored to the lower 
%               half. Its size is ["numberOfEdges", 
%               "numbreOfEpochs"*"numberOfConditions",
%               "numbreOfEpochs"*"numberOfConditions"].
%
%


%% Input checks

% number of args
if ~ismember(nargin, 1:4)
    error('Wrong number of input args - "connData" is needed while "metric", "permNo" and "permStat" are optional!');
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
        elseif isnumeric(varargin{v}) && ismember(varargin{v}, 100:100:10^6)
            permNo = varargin{v};
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
% check size and dimensionality of mandatory arg "connData"
if ~isequal(size(connData, 1), size(connData, 2))
    error('First two dimensions of input arg "connData" need to have equal size!');
end
if length(size(connData)) ~= 4
    error('Input arg "connData" should have four dimensions!');
end

% user message
disp([char(10), 'Called cmpFullConn_edge function with input args:',...
    char(10), 'Input data is array with size ', num2str(size(connData)),...
    char(10), 'Metric: ', metric,...
    char(10), 'No. of random permutations: ', num2str(permNo),...
    char(10), 'Test statistic for random permutations: ', permStat]);


%% Basics

% number of ROIS, conditions and epochs
roiNo = size(connData, 1);
condNo = size(connData, 4);
epochNo = size(connData, 3);
% number of edges = elements in the upper triangle excluding the diagonal
edgeNo = (roiNo*roiNo-roiNo)/2;
% preallocate for linearized data
dataLin = zeros((roiNo*roiNo-roiNo)/2, epochNo, condNo);


%% Linearize and reshape connectivity data

for condIdx = 1:condNo
    for epoch = 1:epochNo
        % get upper triangular part from epoch-level connectivity matrix in
        % a vector
        tmp = squeeze(connData(:, :, epoch, condIdx))';
        idx = tril(true(size(tmp)), -1);  % do not include data from main diagonal
        dataLin(:, epoch, condIdx) = tmp(idx)';
        % warning if there are NaN values
        if any(isnan(dataLin(:, epoch, condIdx)))
            warning(['There is at least one NaN value among the ',...
                'connectivity values at condition ', num2str(condIdx),... 
                ', epoch ', num2str(epoch), '!']);
        end
    end
end  % for condIdx

% conditions become epochs one after another, so e.g. epoch 2 at cond 3
% becomes epoch 22
dataLin = reshape(dataLin, [edgeNo, epochNo*condNo]);

% user message
disp([char(10), 'Vectorized connectivity matrices (upper triangles)']);


%% Normalize connectivity data if method is 'corr'

% for measuring the piecewise contributions of each edge to the total
% result, we look at covariances on normalized data
if strcmp(metric, 'corr')
    tmpMeans = mean(dataLin, 1);
    tmpSds = std(dataLin, 1, 1);
    dataLin = (dataLin-tmpMeans)./tmpSds; % it is still crazy we can do this
end


%% Get similarity across all epoch-pairings, separately for each edge

% preallocate
connSim = nan(edgeNo, epochNo*condNo, epochNo*condNo);  % connectivity pattern similarity across epochs, for each edge

% loop through all edges
for edgeIdx = 1:edgeNo
    % define edge data we work with
    edgeData = dataLin(edgeIdx, :);
    
    % loops through all epoch-pairings
    for epochOne = 1:epochNo*condNo
        for epochTwo = 1:epochNo*condNo
            % only consider epoch-pairings for the upper triangular of the
            % resulting connSim matrix
            if epochOne < epochTwo 
                % calculate similarity according to the supplied metric
                switch metric
                    case 'corr'
                        connSim(edgeIdx, epochOne, epochTwo) = edgeData(epochOne)*edgeData(epochTwo);  % "piecewise covariance" on normalized data
                    case 'eucl'
                        connSim(edgeIdx, epochOne, epochTwo) = abs(edgeData(epochOne)-edgeData(epochTwo));  % "piecewise norm-2"
                end
            end
        end  % for epochTwo
    end  % for epochOne
    
end  % for edgeIdx

% user message
disp([char(10), 'Calculated similarity across all epoch-pairings, all edges']);


%% Compare similarity values between within- and across-condition pairings 

% preallocate results struct

withinCondMean = nan(edgeNo, condNo);
withinCondSD = nan(edgeNo, condNo);
withinCondMedian = nan(edgeNo, 1);
acrossCondMean = nan(edgeNo, condNo);
acrossCondSD = nan(edgeNo, condNo);
acrossCondMedian = nan(edgeNo, condNo);
realDiff = nan(edgeNo, condNo);
permDiffMean = nan(edgeNo, condNo);
permDiffSD = nan(edgeNo, condNo);
pEst = nan(edgeNo, condNo);
cohend = nan(edgeNo, condNo);

% clock for measuring elapsed time
startClock = tic;

% loop through edges
for edgeIdx = 1:edgeNo

    % define data for edge
    edgeData = squeeze(connSim(edgeNo, :, :));
    % mirror the upper triangular part of connSim matrix to the lower
    % half
    edgeData = triu(edgeData, 1) + triu(edgeData, 1)'; 
    
    % loop through conditions / stimuli
    for condIdx = 1:condNo    
        
        % within-cond epoch pairing similarities for given condition
        tmp = edgeData((condIdx-1)*epochNo+1:condIdx*epochNo, (condIdx-1)*epochNo+1:condIdx*epochNo);
        tmp = triu(tmp, 1);  
        % linearize to vector
        tmpT = tmp';
        idx = tril(true(size(tmpT)), -1);
        withinCondSim = tmpT(idx)';

        % across-condition pairings for epochs in given condition
        tmpParts = cell(2, 1);
        if condIdx ~= 1
            tmpParts{1} = edgeData((condIdx-1)*epochNo+1:condIdx*epochNo, 1:(condIdx-1)*epochNo);
        end
        if condIdx ~= condNo
            tmpParts{2} = edgeData((condIdx-1)*epochNo+1:condIdx*epochNo, condIdx*epochNo+1:end);
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
        withinCondMean(edgeIdx, condIdx) = mean(withinCondSim);
        withinCondSD(edgeIdx, condIdx) = std(withinCondSim);
        withinCondMedian(edgeIdx, condIdx) = median(withinCondSim);
        acrossCondMean(edgeIdx, condIdx) = mean(acrossCondSim);
        acrossCondSD(edgeIdx, condIdx) = std(acrossCondSim); 
        acrossCondMedian(edgeIdx, condIdx) = median(acrossCondSim);

        % permutation test across the two groups of connectivity similarity
        % values
        [pEst(edgeIdx, condIdx), realDiff(edgeIdx, condIdx),... 
            permDiff, cohend(edgeIdx, condIdx)] = permTest(withinCondSim, acrossCondSim, permNo, permStat, 'silent');  % suppress permTest outputs to command prompt
        % save out the mean and SD of permited differences
        permDiffMean(edgeIdx, condIdx) = mean(permDiff);
        permDiffSD(edgeIdx, condIdx) = std(permDiff);

    end  % for condIdx
    
    % user message after each Nth edge
    if mod(edgeNo, 50) == 0
        elapsedTime = round(toc(startClock), 2);
        disp([char(10), 'Done with edge no. ', num2str(edgeNo),... 
            ', took ', num2str(elapsedTime), ' secs so far']);
    end
    
end  % for edgeIdx

% user message
disp([char(10), 'Compared similarity within- versus across-condition epoch-pairings']);


%% Save and return

permRes = struct;
permRes.withinCondMean = withinCondMean;
permRes.withinCondSD = withinCondSD;
permRes.withinCondMedian = withinCondMedian;
permRes.acrossCondMean = acrossCondMean;
permRes.acrossCondSD = acrossCondSD;
permRes.acrossCondMedian = acrossCondMedian;
permRes.pEst = pEst;
permRes.realDiff = realDiff;
permRes.permDiffMean = permDiffMean;
permRes.permDiffSD = permDiffSD;
permRes.cohend = cohend;



return














