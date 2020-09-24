function [permRes, withinSim, acrossSim] = cmp2GroupMeanFull(groupData, subIdx, varargin)
%% Compare individual connectivity to leave-one-out group mean connectivity
% Version for full connectivity matrices.
%
% USAGE [permRes, withinSim, acrossSim] = cmp2GroupMeanFull(groupData, subIdx, metric='corr', permNo=10000, permTest='mean')
%
% Function to compare the full connectivity matrix of a subject to that
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
% Only for NON-DIRECTED connectivity, as only the upper triangle of
% connectivity matrices is used.
%
% Mandatory inputs:
% groupData     - Numeric array, sized [ROIs, ROIs, epochs, conditions, subjects]. 
%           Group-level full connectivity matrix. Undirected
%           connectivity, so we work only with values in the upper
%           triangle of each connectivity matrix (connectivity matrix for 
%           a given epoch, condition and subject is defined by 
%           dimensions 1 and 2).
% subIdx        - Numeric value. Index of the subject whose data we compare 
%           with the rest of the group. Obviously needs to be one of
%           1:numberOfSubjects.
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
% permNo        - Numeric value, the number of permutations to perform for
%               random permutation tests. One of 100:100:10^6. Defaults to
%               10^4.
% permStat      - String specifying the statistic we perform the random
%               permutation tests on. One of {'mean', 'median', 'std'} -
%               the ones supported by permTest.m. Defaults to 'mean'.
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
%


%% Input checks

% no. of args
if ~ismember(nargin, 2:5)
    error(['Wrong number of input args - "groupData" and "subIdx" are needed ',...
        'while "metric", "permNo" and "permStat" are optional!']);
end
% loop through varargin to sort out input args
if nargin > 2
    for v = 1:length(varargin)
        if ischar(varargin{v})
            if ismember(varargin{v}, {'corr', 'eucl', 'deltaCon'}) && ~exist('metric', 'var')
                metric = varargin{v};
            elseif ismember(varargin{v}, {'mean', 'median', 'std'}) && ~exist('permStat', 'var')
                permStat = varargin{v};
            else
                error('An input arg (string) could not be mapped to any optional arg!');
            end
        elseif isnumeric(varargin{v}) && ismember(varargin{v}, 100:100:10^6) && ~exist('permNo', 'var')
            permNo = varargin{v};
        else
            error(['At least one input arg could not mapped to nicely to ',...
                'any optional arg ("metric", "permStat" or "permNo")!']);
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
if ~isequal(size(groupData, 1), size(groupData, 2))
    error('First two dimensions of input arg "groupData" need to have equal size!');
end
if length(size(groupData)) ~= 5
    error('Input arg "groupData" should have five dimensions!');
end
% check subject idx
subNo = size(groupData, 5);
if ~ismember(subIdx, 1:subNo)
    error('Input arg "subIdx" is out of the range of subjects in input arg "groupData"!');
end

% user message
disp([char(10), 'Called cmp2GroupMeanFull function with input args:',...
    char(10), 'Group data is array with size ', num2str(size(groupData)),...
    char(10), 'Subject index: ', num2str(subIdx),...
    char(10), 'Metric: ', metric,...
    char(10), 'No. of random permutations: ', num2str(permNo),...
    char(10), 'Test statistic for random permutations: ', permStat]);


%% Define subject- and leave-one-out mean data

% define subject data and leave-out-one group mean
subData = squeeze(groupData(:, :, :, :, subIdx));
groupData(:, :, :, :, subIdx) = [];
groupData = mean(groupData, 5);

% number of ROIS, conditions and epochs
roiNo = size(groupData, 1);
condNo = size(groupData, 4);
epochNo = size(groupData, 3);

% reshape data so that epochs across conditions/stimuli come after each
% other
subData = reshape(subData, [roiNo, roiNo, epochNo*condNo]);
groupData = reshape(groupData, [roiNo, roiNo, epochNo*condNo]);

% user message
disp([char(10), 'Defined subject and leave-one-out group data']);


%% Get similarity score from all epoch-pairings

% calculation depends on metric type

% simple correlation
if strcmp(metric, 'corr')
    % extract upper triangles above the main diagonal
    subDataLin = linearizeTrius(subData, 1);  
    groupDataLin = linearizeTrius(groupData, 1);
    % get all column-pairwise correlations
    connSim = corr(subDataLin, groupDataLin);
    % keep the upper triangle, set the rest to NaN
    connSim(tril(true(epochNo*condNo), -1)) = nan;
    
% frobenius norm of difference matrix ('eucl') or DeltaCon
elseif ismember(metric, {'eucl', 'deltaCon'})
    % preallocate results matrix
    connSim = nan(epochNo*condNo);
    % loops through subject and group epochs
    for subEpoch = 1:epochNo*condNo
        for groupEpoch = 1:epochNo*condNo
            if subEpoch <= groupEpoch  % only for upper triangle of pairings
                
                % create symmetric adjacency matrices with zeros at diagonal
                subEpochData = triu(subData(:, :, subEpoch), 1) + triu(subData(:, :, subEpoch), 1)';
                groupEpochData = triu(groupData(:, :, subEpoch), 1) + triu(groupData(:, :, subEpoch), 1)';
                
                % calculation depending on metric
                switch metric
                    case 'eucl'
                        connSim(subEpoch, groupEpoch) = norm(subEpochData-groupEpochData, 'fro');
                    case 'deltaCon'
                        connSim(subEpoch, groupEpoch) = deltaCon(subEpochData, groupEpochData, false);  % verbosity of deltaCon is set to false
                end
                
            end  % if
        end  % for groupEpoch
    end  % for subEpoch    
    
end
       
% user message
disp([char(10), 'Calculated epoch-pairing distance / similarities']);


%% Extract within-epoch and across-epoch values, compare them

% within-epoch and across-epoch values into vectors
withinSim = diag(connSim);
acrossSim = connSim(triu(true(epochNo*condNo), 1));

% results into a struct
permRes = struct;

% compare within- and across- with a random permutation test
[permRes.pEst, permRes.realDiff,... 
    permRes.permDiff, permRes.CohenD] = permTest(withinSim, acrossSim, permNo, permStat, 'silent'); 

% user message
disp([char(10), 'Finished permutation test:']);
disp(permRes);


return








