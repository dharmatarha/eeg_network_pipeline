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
% data. See "help cmp2GroupMeanFull" for further details.
%
% Only for NON-DIRECTED connectivity, as only the upper triangle of
% connectivity matrices is used.
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
%               matrix comparisons. One of {'corr', 'eucl', 'deltaCon'} 
%               which stand for (1) pearson correlation, (2) eucledian 
%               (frobenius for matrix) norm and (3) DeltaCon, 
%               a network similarity measure (see "help deltaCon" for 
%               details). If 'corr', connectivity matrices are first 
%               vectorized, then simple Pearson correlation is used.
%               Defaults to 'corr'. 
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
    if ~ismember(metric ,{'corr', 'eucl', 'deltaCon'})
        error('Input arg "metric" must be one of {''corr'', ''eucl'', ''deltaCon''}!');
    end
end
% check size and dimensionality of mandatory arg "groupData"
if ~isequal(size(groupData, 1), size(groupData, 2))
    error('First two dimensions of input arg "groupData" need to have equal size!');
end
if length(size(groupData)) ~= 5
    error('Input arg "groupData" should have five dimensions!');
end

% user message
disp([char(10), 'Called group_cmp2GroupMeanFull function with input args:',...
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
    
    % reshape data so that epochs across conditions/stimuli come after each
    % other
    subData = reshape(subData, [roiNo, roiNo, epochNo*condNo]);
    meanData = reshape(meanData, [roiNo, roiNo, epochNo*condNo]);

    % Get similarity score from all epoch-pairings
    % calculation depends on metric type

    % simple correlation
    if strcmp(metric, 'corr')
        % extract upper triangles above the main diagonal
        subDataLin = linearizeTrius(subData, 1);  
        meanDataLin = linearizeTrius(meanData, 1);
        % get all column-pairwise correlations
        connSim(:, :, subIdx) = corr(subDataLin, meanDataLin);
        % keep the upper triangle, set the rest to NaN
        connSim(tril(true(epochNo*condNo), -1)) = nan;
    
    % frobenius norm of difference matrix ('eucl') or DeltaCon
    elseif ismember(metric, {'eucl', 'deltaCon'})
        % loops through subject and group epochs
        for subEpoch = 1:epochNo*condNo
            for groupEpoch = 1:epochNo*condNo
                if subEpoch <= groupEpoch  % only for upper triangle of pairings

                    % create symmetric adjacency matrices with zeros at diagonal
                    subEpochData = triu(subData(:, :, subEpoch), 1) + triu(subData(:, :, subEpoch), 1)';
                    meanEpochData = triu(meanData(:, :, groupEpoch), 1) + triu(meanData(:, :, groupEpoch), 1)';

                    % calculation depending on metric
                    switch metric
                        case 'eucl'
                            connSim(subEpoch, groupEpoch, subIdx) = norm(subEpochData-meanEpochData, 'fro');
                        case 'deltaCon'
                            connSim(subEpoch, groupEpoch, subIdx) = deltaCon(subEpochData, meanEpochData, false);  % verbosity of deltaCon is set to false
                    end

                end  % if
            end  % for groupEpoch
        end  % for subEpoch    

    end

    % user message about progress
    disp([char(10), 'Calculation finished for subject ', num2str(subIdx), char(10)]);

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
    tmp = connSim(:, :, subIdx);
    idx = triu(true(size(tmp)), 1);
    diffEpochSim(:, subIdx) = tmp(idx);
end
    
% reshape data for random-effects model (anovan)
tmpSame = reshape(sameEpochSim, [numel(sameEpochSim), 1]);
tmpDiff = reshape(diffEpochSim, [numel(diffEpochSim), 1]);
yData = [tmpSame; tmpDiff];

% grouping variables
sameVSdiff = [zeros(numel(sameEpochSim), 1); ones(numel(diffEpochSim), 1)];  % same-epoch vs different-epoch pairings
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









