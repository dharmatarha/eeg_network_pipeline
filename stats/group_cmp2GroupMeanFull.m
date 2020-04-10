function anovaRes = group_cmp2GroupMeanFull(groupData, metric)
%% Group-level stats for individual vs leave-one-out group connectivity comparisons
% Version for full connectivity matrices.
%
% USAGE anovaRes = group_cmp2GroupMeanFull(groupData, metric='corr')
%
% Same as cmp2GroupMeanFull.m but working on the group-level: it calculates
% the within-epoch and across-epoch connectivity similarities between
% individual and leave-one-out group datasets for all individuals and fits
% a linear model (anovan with subject as random effect) to the aggregate
% data.
%
% The anovan output is also displayed in a figure.
%
% Mandatory inputs:
% groupData     - Numeric array, sized [ROIs, ROIs, epochs, conditions, subjects]. 
%           Group-level full connectivity matrix. Undirected
%           connectivity, so we work only with values in the upper
%           triangle of each connectivity matrix (connectivity matrix for 
%           a given epoch, condition and subject is defined by 
%           dimensions 1 and 2).
%
% Optional inputs:
% metric        - String specifying distance metric for connectivity
%               matrix comparisons. Matrices are first vectorized, 
%               then one of these distances is used: {'corr', 'eucl'}.  
%
% Output:
% anovaRes      - Struct containing the results of anovan performed on the
%               whole dataset. Model is "full" with (1) fixed effect for
%               within- vs. across-epoch pairing similarities and (2) random
%               effect for subject numbers. Fields are "p", "table" and
%               "stats"
%


%% Input checks

% no. of args
if ~ismember(nargin, 1:2)
    error(['Wrong number of input args - "groupData" is needed ',...
        'while "metric" is optional!']);
end
% loop through varargin to sort out input args
if nargin == 1
    metric = 'corr';
else
    if ~ismember(metric ,{'corr', 'eucl'})
        error('Input arg "metric" must be one of {''corr'', ''eucl''}!');
    end
end
% check size and dimensionality of mandatory arg "connData"
if ~isequal(size(groupData, 1), size(groupData, 2))
    error('First two dimensions of input arg "groupData" need to have equal size!');
end
if length(size(groupData)) ~= 5
    error('Input arg "groupData" should have five dimensions!');
end

% user message
disp([char(10), 'Called cmp2GroupMeanFull function with input args:',...
    char(10), 'Group data is array with size ', num2str(size(groupData)),...
    char(10), 'Metric: ', metric,]);


%% Basics, settings

% number of ROIS, conditions and epochs
roiNo = size(groupData, 1);
epochNo = size(groupData, 3);
condNo = size(groupData, 4);
subNo = size(groupData, 5);

% preallocate results matrix
connSim = nan(epochNo*condNo, epochNo*condNo, subNo);


%% Loop through subjects, calculate subject vs leave-one-out group-mean similarity

for subIdx = 1:subNo
    
    % user message about progress
    disp([char(10), 'Calculating for subject ', num2str(subIdx)]);
    
    % Define subject- and leave-one-out mean data
    subData = squeeze(groupData(:, :, :, :, subIdx));
    meanData = groupData;
    meanData(:, :, :, :, subIdx) = [];
    meanData = mean(meanData, 5);
    
    % preallocate for linearized data
    subDataLin = zeros((roiNo*roiNo-roiNo)/2, epochNo, condNo);
    meanDataLin = subDataLin;
    
    % vectorize (linearize) data
    for cond = 1:condNo
        for epoch = 1:epochNo
            % get upper triangular part from epoch-level connectivity matrix in
            % a vector, for subject and group separately
            % subject
            tmp = squeeze(subData(:, :, epoch, cond))';
            idx = tril(true(size(tmp)), -1);  % do not include data from main diagonal
            subDataLin(:, epoch, cond) = tmp(idx)';
            % group
            tmp = squeeze(meanData(:, :, epoch, cond))';
            idx = tril(true(size(tmp)), -1);  % do not include data from main diagonal
            meanDataLin(:, epoch, cond) = tmp(idx)';        
            % warning if there are NaN values
            if any(isnan(subDataLin(:, epoch, cond))) || any(isnan(meanDataLin(:, epoch, cond)))
                warning(['There is at least one NaN value among the ',...
                    'connectivity values at condition ', num2str(cond),... 
                    ', epoch ', num2str(epoch), '!']);
            end
        end
    end    
    
    % conditions become epochs one after another, so e.g. epoch 2 at cond 3
    % becomes epoch 22
    subDataLin = reshape(subDataLin, [(roiNo*roiNo-roiNo)/2, epochNo*condNo]);
    meanDataLin = reshape(meanDataLin, [(roiNo*roiNo-roiNo)/2, epochNo*condNo]);

    % Get similarity score from all epoch-pairings
    % calculation depends on metric type
    switch metric

        case 'corr'  % simple correlation
            % get all column-pairwise correlations
            tmp = corr(subDataLin, meanDataLin);
            % keep the upper triangle, set the rest to NaN
            tmp(tril(true(epochNo*condNo), -1)) = nan;
            connSim(:, :, subIdx) = tmp;

        case 'eucl'  % 2-norm of differences

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THIS PART COULD BE SIMPLIFIED WITH BUILT-IN FUNCTION "VECNORM"
            % CALCULATION BELOW IS CURRENTLY PREFERRED FOR
            % BACKWARDS-COMPATIBILITY WITH OLDER MATLAB VERSIONS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % preallocate results matrix
            connSim = nan(epochNo*condNo);
            % loops through subject and group epochs
            for subEpoch = 1:epochNo*condNo
                for groupEpoch = 1:epochNo*condNo
                    if subEpoch <= groupEpoch  % only for upper triangle of pairings
                        % norm of diff vector
                        connSim(subEpoch, groupEpoch, subIdx) = norm(subDataLin(:, subEpoch)-meanDataLin(:, groupEpoch));              
                    end       
                end  % groupEpoch for
            end  % subEpoch for

    end  % switch metric
    
end


%% Stats on the difference between same-epoch vs different-epoch pairings

% user message
disp([char(10), 'Calling anovan on aggregated results for group-level stats']);

% preallocate
sameEpochSim = nan(epochNo*condNo, subNo);
diffEpochSim = nan((epochNo^2*condNo^2-epochNo*condNo)/2, subNo);
anovaRes = struct;

% get same-epoch pairing and different-epoch pairing similarities
for subIdx = 1:subNo
    sameEpochSim(:, subIdx) = diag(connSim(:, :, subIdx));
    tmp = connSim(:, :, subIdx)';
    idx = tril(true(size(tmp)), -1);
    diffEpochSim(:, subIdx) = tmp(idx)';
end
    
% reshape data for random-effects model (anovan)
tmpSame = reshape(sameEpochSim, [prod(size(sameEpochSim)), 1]);
tmpDiff = reshape(diffEpochSim, [prod(size(diffEpochSim)), 1]);
yData = [tmpSame; tmpDiff];

% grouping variables
sameVSdiff = [zeros(prod(size(sameEpochSim)), 1); ones(prod(size(diffEpochSim)), 1)];  % same-epoch vs different-epoch pairings
tmpSize1 = size(sameEpochSim, 1);
tmpSize2 = size(diffEpochSim, 1);
for subIdx = 1:subNo
    sIdxPart1(tmpSize1*(subIdx-1)+1:tmpSize1*subIdx, 1) = repmat(subIdx, [tmpSize1, 1]);
    sIdxPart2(tmpSize2*(subIdx-1)+1:tmpSize2*subIdx, 1) = repmat(subIdx, [tmpSize2, 1]);
end
subjectIdx = [sIdxPart1; sIdxPart2];  % subject numbers
vGrouping = {sameVSdiff, subjectIdx};

% linear model with random effect for subject no. 
[anovaRes.p, anovaRes.table, anovaRes.stats] = anovan(yData, vGrouping,... 
    'random', 2, 'model', 'full', 'varnames', {'EpochPairingType', 'SubjectNo'});

% user message
disp([char(10), 'Done, bye!']);

return









