function [permRes, withinCondPermRes, connSim] = cmpFullConn(connData, varargin)
%% Compare full connectivity matrices across conditions
%
% USAGE: [permRes, withinCondPermRes, connSim] = cmpFullConn(connData, metric='corr', permNo=10000, permStat='mean')
%
% Compares connectivity patterns across different groups of epoch-pairings 
% using permutation tests. First, for each condition, it tests the
% within-condition epoch-pairings against the across-condition epoch
% pairings. Second, for each condition-pairing, it tests if epoch-pairings
% in one condition have different values than in the other condition 
% (as in the case of pairwise comparisons). 
%
% Input data is expected to contain undirected connectivity values. From
% each connectivity matrix we only take the upper triangle into account.
%
% For the tests we rely on permTest.m, input args "permNo" and "permStat"
% (if specified) are passed to permTest.m. Tries to call permTest.m with
% 'studentized' flag (only supported for 'mean').
%
% Mandatory input:
% connData      - 4D numeric array, sets of connectivity matrices across
%               epochs and conditions (stimuli). First two dimensions have
%               equal size and determine a connectivity matrix, third
%               dimension is epochs, fourth is conditions (stimuli).
%
% Optional inputs:
% metric        - Char array specifying distance metric for connectivity
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
% permStat      - Char array specifying the statistic we perform the random
%               permutation tests on. One of {'mean', 'median', 'std'} -
%               the ones supported by permTest.m. Defaults to 'mean'.
%
% Outputs:
% permRes       - Struct, random permutation test results. The first 
%               "condition no." elements of permRes summerize the results 
%               of the within-condition epoch-pairings versus 
%               across-condition epoch-pairings test for each condition. 
%               The subsequent, last element contains the results for 
%               comparing all within-condition pairings to all 
%               across-condition pairings (across all stimuli). Thus, 
%               its length is "numberOfConditions"+1.
%               Has fields for the mean, median and SD values of both
%               epoch-pairing groups, and for the outcomes of permTest,
%               i.e. estimated p-value, real test stat difference, 
%               permuted difference values and effect size (Cohen's d).
% withinCondPermRes       - Struct, random permutation test results. 
%               Each element of withinCondPermRes summerizes the results of
%               a comparison of epoch-pairings across two conditions, i.e. 
%               its length is nchoosek("condition no.", 2).
%               Has fields for the mean, median and SD values of both
%               epoch-pairing groups, and for the outcomes of permTest,
%               i.e. estimated p-value, real test stat difference,
%               permuted difference values and effect size (Cohen's d).
% connSim       - Numeric matrix containing connectivity similarity values
%               for all epoch-pairings. Symmetric, as the upper triangle
%               values are calculated and mirrored to the lower half. Its
%               size is ["epoch no" X "condition no.",
%               "epoch no." X "condition no."].
%
%


%% Input checks

% number of args
if ~ismember(nargin, 1:4)
    error('Function cmpFullConn requires input arg "connData" while args "metric", "permNo" and "permStat" are optional!');
end
% check size and dimensionality of mandatory arg "connData"
if ~isnumeric(connData) || length(size(connData)) ~= 4
    error('Input arg "connData" should be 4D numeric array!');
end
if ~isequal(size(connData, 1), size(connData, 2))
    error('First two dimensions of input arg "connData" need to have equal sizes!');
end
% loop through varargin to sort out input args
if ~isempty(varargin)
    for v = 1:length(varargin)
        if ischar(varargin{v})
            if ismember(varargin{v}, {'corr', 'eucl', 'deltaCon'}) && ~exist('metric', 'var')
                metric = varargin{v};
            elseif ismember(varargin{v}, {'mean', 'median', 'std'}) && ~exist('permStat', 'var')
                permStat = varargin{v};
            else
                error('An input arg (char array) could not be mapped to either "metric" or "permStat"!');
            end
        elseif isnumeric(varargin{v}) && ismember(varargin{v}, 100:100:10^6) && ~exist('permNo', 'var')
            permNo = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to any optional arg!');
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

% user message
disp([char(10), 'Called cmpFullConn function with input args:',...
    char(10), 'Input data is array with size ', num2str(size(connData)),...
    char(10), 'Metric: ', metric,...
    char(10), 'No. of random permutations: ', num2str(permNo),...
    char(10), 'Test statistic for random permutations: ', permStat]);


%% Basics

% number of ROIS, conditions and epochs
[roiNo, ~, epochNo, condNo] = size(connData);

% reshape data so that epochs across conditions/stimuli come after each
% other
connData = reshape(connData, [roiNo, roiNo, epochNo*condNo]);


%% Get similarity score from all epoch-pairings

% calculation depends on metric type

% simple correlation
if strcmp(metric, 'corr')
    % extract upper triangles above the main diagonal
    dataLin = linearizeTrius(connData, 1);  
    % get all column-pairwise correlations
    connSim = corr(dataLin, dataLin);
    % keep the upper triangle, inlcuding main diagonal, set the rest to NaN
    connSim(tril(true(epochNo*condNo), -1)) = nan;
    
% frobenius norm of difference matrix ('eucl') or DeltaCon
elseif ismember(metric, {'eucl', 'deltaCon'})
    % preallocate results matrix
    connSim = nan(epochNo*condNo);
    % loops through epochs
    for epochOne = 1:epochNo*condNo
        % create symmetric adjacency matrices with zeros at diagonal
        epochOneData = triu(connData(:, :, epochOne), 1) + triu(connData(:, :, epochOne), 1)';
        for epochTwo = 1:epochNo*condNo
            if epochOne <= epochTwo  % only for upper triangle of pairings
                
                % create symmetric adjacency matrices with zeros at diagonal
                epochTwoData = triu(connData(:, :, epochTwo), 1) + triu(connData(:, :, epochTwo), 1)';
                
                % calculation depending on metric
                switch metric
                    case 'eucl'
                        connSim(epochOne, epochTwo) = norm(epochOneData-epochTwoData, 'fro');
                    case 'deltaCon'
                        connSim(epochOne, epochTwo) = deltaCon(epochOneData, epochTwoData, false);  % verbosity of deltaCon is set to false
                end
                
            end  % if
        end  % for groupEpoch
    end  % for subEpoch      
    
end

% user message
disp([char(10), 'Calculated epoch-pairing distance / similarities']);


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
permRes.studentDiff = nan;
permRes.permDiff = nan(permNo, 1);
permRes.pEst = nan;

% preallocate vars holding within- and across condition/stimulus similarities
withinCondSim = nan((epochNo^2-epochNo)/2, condNo);
acrossCondSim = nan(epochNo^2*(condNo-1), condNo);

for condIdx = 1:condNo
    
    % within-cond epoch pairing similarities for given condition
    tmp = connSim((condIdx-1)*epochNo+1:condIdx*epochNo, (condIdx-1)*epochNo+1:condIdx*epochNo); 
    withinCondSim(:, condIdx) = tmp(triu(true(epochNo), 1));
    
    % across-condition pairings for epochs in given condition
    tmpParts = cell(2, 1);
    if condIdx ~= 1
        tmpParts{1} = connSim((condIdx-1)*epochNo+1:condIdx*epochNo, 1:(condIdx-1)*epochNo);
    end
    if condIdx ~= condNo
        tmpParts{2} = connSim((condIdx-1)*epochNo+1:condIdx*epochNo, condIdx*epochNo+1:end);
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
    acrossCondSim(:, condIdx) = reshape(tmpAll', [epochNo^2*(condNo-1), 1]);
    
    % store descriptives in result struct
    permRes(condIdx).withinCondMean = mean(withinCondSim(:, condIdx));
    permRes(condIdx).withinCondSD = std(withinCondSim(:, condIdx));
    permRes(condIdx).withinCondMedian = median(withinCondSim(:, condIdx));
    permRes(condIdx).acrossCondMean = mean(acrossCondSim(:, condIdx));
    permRes(condIdx).acrossCondSD = std(acrossCondSim(:, condIdx)); 
    permRes(condIdx).acrossCondMedian = median(acrossCondSim(:, condIdx));
    
    % permutation test across the two groups of connectivity similarity
    % values
    [permRes(condIdx).pEst,... 
     permRes(condIdx).realDiff,... 
     permRes(condIdx).permDiff,... 
     permRes(condIdx).cohenD,...
     permRes(condIdx).studentDiff] = permTest(withinCondSim(:, condIdx),... 
                                         acrossCondSim(:, condIdx),... 
                                         permNo,... 
                                         permStat,...
                                         'studentized',...
                                         'silent');
    
end

% user message
disp([char(10), 'Compared similarity within- versus across-condition epoch-pairings']);


%% Comparison between within- and across-condition pairings across all conditions / stimuli

% store descriptives in result struct
permRes(condNo+1).withinCondMean = mean(withinCondSim(:));
permRes(condNo+1).withinCondSD = std(withinCondSim(:));
permRes(condNo+1).withinCondMedian = median(withinCondSim(:));
permRes(condNo+1).acrossCondMean = mean(acrossCondSim(:));
permRes(condNo+1).acrossCondSD = std(acrossCondSim(:)); 
permRes(condNo+1).acrossCondMedian = median(acrossCondSim(:));

% permutation test on vectorized withinCondSim and acrossCondSim matrices
[permRes(condNo+1).pEst,... 
 permRes(condNo+1).realDiff,... 
 permRes(condNo+1).permDiff,... 
 permRes(condNo+1).cohenD,...
 permRes(condNo+1).studentDiff] = permTest(withinCondSim(:),... 
                                      acrossCondSim(:),... 
                                      permNo,... 
                                      permStat,...
                                      'studentized',...
                                      'silent');

% user message
disp([char(10), 'Compared within- vs across-cond similarity across all stimuli']);    
    

%% Compare similarity values between different within-condition pairings
% e.g. between epoch-pairings in within cond 1 and epoch-pairings within
% cond 2

% preallocate results struct
withinCondPermRes = struct;

% loops through conditions
counter = 0;
for condOne = 1:condNo
    for condTwo = 1:condNo
        % only calculate if the two conditions are not the same
        if condTwo > condOne
            
            % adjust counter
            counter = counter + 1;
            
            % permutation test across the two groups of connectivity similarity
            % values
            [withinCondPermRes(counter).pEst,... 
                withinCondPermRes(counter).realDiff,... 
                withinCondPermRes(counter).permDiff,...
                withinCondPermRes(counter).cohenD,...
                withinCondPermRes(counter).studentDiff] = permTest(withinCondSim(:, condOne),... 
                                                              withinCondSim(:, condTwo),... 
                                                              permNo,... 
                                                              permStat,...
                                                              'studentized',...
                                                              'silent');
            
        end  % if 
    end  % for condTwo
end  % for condOne

% user message
disp([char(10), 'Compared similarity across conditions for within-condition epoch-pairings']);
disp('Done with everything!');


return





