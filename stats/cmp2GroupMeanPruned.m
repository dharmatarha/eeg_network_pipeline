function [permRes, withinSim, acrossSim] = cmp2GroupMeanPruned(groupData, subIdx, prunedMask, varargin)
%% Compare individual connectivity to leave-one-out group mean connectivity
% Version for pruned (sparse) connectivity matrices.
%
% USAGE [permRes, withinSim, acrossSim] = cmp2GroupMeanPruned(groupData, subIdx, prunedMask, metric='corr', permNo=10000, permTest='mean')
%
% Function to compare the pruned connectivity matrix of a subject to that
% of the whole group except that subject (leave-one-out mean). 
% We evaluate the null hypothesis that the similarity of individual and
% group connectivity in matching epoch-pairings is the same as in 
% non-matching epoch-pairings.
% That is, similarity calculated from epochs e.g. 1-1, 2-2, ..., 
% 10-10, 11-11,... for the individual
% and the group, respectively, is the same as from epochs e.g. 1-30, 2-5,
% etc.
%
% We calculate both groups of similarity values (within-epoch vs
% across-epochs) and compare the two sets of values with a permutation test
%
% Mandatory inputs:
% groupData     - Numeric array, sized [ROIs, ROIs, epochs, conditions, subjects]. 
%           Group-level pruned connectivity matrix. Undirected
%           connectivity, so we work only with values in the upper
%           triangle of each connectivity matrix (connectivity matrix for 
%           a given epoch, condition and subject is defined by 
%           dimensions 1 and 2).
% subIdx        - Numeric value. Index of the subject whose data we compare 
%           with the rest of the group. Obviously needs to be one of
%           1:numberOfSubjects.
% prunedMask    - Numeric binary matrix, normally the output of
%           getPrunedMask (called with groupData). Contains ones for edges
%           to consider and zeros for all other connections/edges.
%
% Optional inputs:
% metric        - String specifying distance metric for connectivity
%               matrix comparisons. Matrices are first vectorized, 
%               then one of these distances is used: {'corr', 'eucl'}.  
% permNo        - Numeric value, the number of permutations to perform for
%               random permutation tests. One of 100:100:10^6.
% permStat      - String specifying the statistic we perform the random
%               permutation tests on. One of {'mean', 'median', 'std'} -
%               the ones supported by permTest.m
%
% Outputs:
% permRes       - Struct containing the output of permTest.m. Its fields
%               are pEst (estimated prob.), realDiff (the real difference
%               between within-epoch and across-epoch pairing similarities
%               in terms of the test stat), permDiff (the permuted
%               difference values) and CohenD (effect size, Cohen's d).
% withinSim     - Numeric vector of within-epoch pairing connectivity
%               similarities. Its size is
%               [numberOfEpochs*numberOfConditions, 1]
% acrossSim     - Numeric vector of across-epoch pairing connectivity
%               similarities. Its size is
%               [(numberOfEpochs^2*numberOfConditions^2-numberOfEpochs*numberOfConditions)/2, 1]


%% Input checks

% no. of args
if ~ismember(nargin, 3:6)
    error(['Wrong number of input args - "groupData", "subIdx" and "prunedMask" are needed ',...
        'while "metric", "permNo" and "permStat" are optional!']);
end
% loop through varargin to sort out input args
if nargin > 3
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
% check size and dimensionality of mandatory arg "groupData"
tmp = size(groupData);
if ~isequal(tmp(1), tmp(2))
    error('First two dimensions of input arg "groupData" need to have equal size!');
end
if length(tmp) ~= 5
    error('Input arg "groupData" should have five dimensions!');
end
% check subject idx
subNo = tmp(5);
if ~ismember(subIdx, 1:subNo)
    error('Input arg "subIdx" is out of the range of subjects in input arg "groupData"!');
end
% check pruned data mask (size, binary)
if ~isequal(size(prunedMask), [tmp(1), tmp(2)])
    error('Size of input arg "prunedMask" does not match the size of the first two dims of "groupData"!');
end
if ~isequal(unique(prunedMask), [0 1]')
    error('Input arg "prunedMask" should be binary!');
end

% user message
disp([char(10), 'Called cmp2GroupMeanFull function with input args:',...
    char(10), 'Group data is array with size ', num2str(size(groupData)),...
    char(10), 'Subject index: ', num2str(subIdx),...
    char(10), 'Data mask is array with size ', num2str(size(prunedMask)),...
    char(10), 'Number of edges in mask: ', numstr(sum(sum(prunedMask))),...
    char(10), 'Metric: ', metric,...
    char(10), 'No. of random permutations: ', num2str(permNo),...
    char(10), 'Test statistic for random permutations: ', permStat]);


%% Define subject- and leave-one-out mean data

% define subject data and leave-out-one group mean
subData = squeeze(groupData(:, :, :, :, subIdx));
subData(isnan(subData)) = 0;  % turn NaN values to zero
groupData(:, :, :, :, subIdx) = [];
groupData(isnan(groupData)) = 0;  % turn Nan values to zero
groupData = mean(groupData, 5);

% number of ROIS, conditions and epochs
roiNo = size(groupData, 1);
condNo = size(groupData, 4);
epochNo = size(groupData, 3);

% preallocate for linearized data
subDataLin = zeros(sum(sum(prunedMask)), epochNo, condNo);
groupDataLin = subDataLin;

% user message
disp([char(10), 'Defined subject and leave-one-out group data']);


%% Linearize and reshape connectivity data

% for each epoch in each condition, use the mask to both select the proper
% values and to linearize subject and mean group data
for cond = 1:condNo
    for epoch = 1:epochNo
        tmp = subData(:, :, epoch, cond);  % we need to define a single matrix for properly indexing with mask
        subDataLin(:, epoch, cond) = tmp(prunedMask);
        tmp = groupData(:, :, epoch, cond);
        groupDataLin(:, epoch, cond) = tmp(prunedMask);
    end
end

% conditions become epochs one after another, so e.g. epoch 2 at cond 3
% becomes epoch 22
subDataLin = reshape(subDataLin, [sum(sum(prunedMask)), epochNo*condNo]);
groupDataLin = reshape(groupDataLin, [sum(sum(prunedMask)), epochNo*condNo]);

% user message
disp([char(10), 'Vectorized connectivity matrices (upper triangles), masked all connectivity data']);


%% Get similarity score from all epoch-pairings

% calculation depends on metric type
switch metric
    
    case 'corr'  % simple correlation
        % get all column-pairwise correlations
        connSim = corr(subDataLin, groupDataLin);
        % keep the upper triangle, set the rest to NaN
        connSim(tril(true(epochNo*condNo), -1)) = nan;
        
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
                    connSim(subEpoch, groupEpoch) = norm(subDataLin(:, subEpoch)-groupDataLin(:, groupEpoch));              
                end       
            end  % groupEpoch for
        end  % subEpoch for
        
end  % switch metric

% user message
disp([char(10), 'Calculated epoch-pairing similarities']);


%% Extract within-epoch and across-epoch values, compare them

% within-epoch and across-epoch values into vectors
withinSim = diag(connSim);
acrossSim = connSim(triu(true(epochNo*condNo), 1));

% results into a struct
permRes = struct;

% compare within- and across- with a random permutation test
[permRes.pEst, permRes.realDiff,... 
    permRes.permDiff, permRes.CohenD] = permTest(withinSim, acrossSim, permNo, permStat); 


return




