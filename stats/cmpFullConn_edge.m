function [permRes, withinCondAnova, connSim] = cmpFullConn_edge(connData, varargin)
%% Compare edge contributions to connectivity matrice differences across conditions
%
% USAGE: [permRes, withinCondAnova, connSim] = cmpFullConn_edge(connData, 
%                                                               metric='corr', 
%                                                               permNo=10000, 
%                                                               permStat='mean')
%
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
% (if specified) are passed to permTest.m. Tries to call permTest.m with
% 'studentized' flag (only supported for 'mean').
%
% USES PARFOR FOR THE LOOP ACROSS EDGES!
%
% Mandatory input:
% connData      - 4D numeric array, sets of connectivity matrices across
%               epochs and conditions (stimuli). First two dimensions have
%               equal size and determine a connectivity matrix, third
%               dimension is epochs, fourth is conditions (stimuli).
%
% Optional inputs:
% metric        - Char array specifying distance metric for connectivity
%               matrix comparisons. One of {'corr', 'eucl'} 
%               which stand for (1) pearson correlation, (2) eucledian 
%               (frobenius for matrix) norm. If 'corr', connectivity 
%               matrices are first vectorized, then simple Pearson 
%               correlation is used. Defaults to 'corr'. 
% permNo        - Numeric value, the number of permutations to perform for
%               random permutation tests. One of 100:100:10^6. Defaults to
%               10^4.
% permStat      - String specifying the statistic we perform the random
%               permutation tests on. One of {'mean', 'median', 'std'} -
%               the ones supported by permTest.m. Defaults to 'mean'.
%
% Outputs:
% permRes       - Struct, random permutation test results for each edge. 
%               It has fields for the results of permutation tests 
%               comparing within-condition epoch-pairings versus 
%               across-condition epoch-pairings for each condition:
%               withinCondMean, withinCondSD, withinCondMedian,
%               acrossCondMean, acrossCondSD, acrossCondMedian, realDiff, 
%               permDiffMean, permDiffSD, pEst, cohend. 
%               Each of these fields is a numeric matrix sized 
%               ["no. of edges" X "no of conditions"].
%               Also contains fields for the results of within vs across
%               comparisons across all stimuli/conditions:
%               withinAllMean, withinAllSD, withinAllMedian,
%               acrossAllMean, acrossAllSD, acrossAllMedian, realDiffAll, 
%               permDiffAllMean,  permDiffAllSD, pEstAll, cohendAll.
%               These are all numeric vectors sized ["no of edges" X 1].
% withinCondAnova       - Struct, anova results for comparing within-cond 
%               epoch-pairings across stimuli/conditions, separately for 
%               each edge. 
%               Its length is nchoosek("condition no.", 2).
%               Has fields for the mean, median and SD values of both
%               epoch-pairing groups, and for the outcomes of permTest,
%               i.e. estimated p-value, real test stat difference,
%               permuted difference values and effect size (Cohen's d).
% connSim       - 3D numeric array containing connectivity similarity values
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
    error('Function cmpFullConn_edge requires input arg "connData" while args "metric", "permNo" and "permStat" are optional!');
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
            if ismember(varargin{v}, {'corr', 'eucl'}) && ~exist('metric', 'var')
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
disp([char(10), 'Called cmpFullConn_edge function with input args:',...
    char(10), 'Input data is numeric array with size ', num2str(size(connData)),...
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

% clock for measuring time for similarity measurements
simClock = tic;

% preallocate
% connectivity pattern similarity across epochs, for each edge
edgeNo = roiNo*(roiNo-1)/2;
connSim = nan(edgeNo, epochNo*condNo, epochNo*condNo);  

% extract upper triangles above the main diagonal
dataLin = linearizeTrius(connData, 1); 

% extra preprocessing steps are needed for correlations - normalizations
% before looking at covariance
if strcmp(metric, 'corr')
    % for measuring the piecewise contributions of each edge to the total
    % result, we look at covariances on normalized data
    tmpMeans = mean(dataLin, 1);
    tmpSds = std(dataLin, 1, 1);
    dataLin = (dataLin-tmpMeans)./tmpSds; % it is still crazy we can do this    
end

% calculation depends on metric type
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
simTime = round(toc(simClock), 2);
disp([char(10), 'Calculated similarity for all epoch-pairings, ',...
    'for all edges, took ', num2str(simTime), ' secs!']);
    

%% Compare similarity values between within- and across-condition pairings, for all edges

% preallocate results vars
withinCondMean = nan(edgeNo, condNo);
withinCondSD = nan(edgeNo, condNo);
withinCondMedian = nan(edgeNo, condNo);
acrossCondMean = nan(edgeNo, condNo);
acrossCondSD = nan(edgeNo, condNo);
acrossCondMedian = nan(edgeNo, condNo);
realDiff = nan(edgeNo, condNo);
studentDiff = nan(edgeNo, condNo);
permDiffMean = nan(edgeNo, condNo);
permDiffSD = nan(edgeNo, condNo);
pEst = nan(edgeNo, condNo);
cohend = nan(edgeNo, condNo);
withinAllMean = nan(edgeNo, 1);
withinAllSD = nan(edgeNo, 1);
withinAllMedian = nan(edgeNo, 1);
acrossAllMean = nan(edgeNo, 1);
acrossAllSD = nan(edgeNo, 1);
acrossAllMedian = nan(edgeNo, 1);
realDiffAll = nan(edgeNo, 1);
studentDiffAll = nan(edgeNo, 1);
permDiffAllMean = nan(edgeNo, 1);
permDiffAllSD = nan(edgeNo, 1);
pEstAll = nan(edgeNo, 1);
cohendAll = nan(edgeNo, 1);
anovaP = nan(edgeNo, 1);
anovaTab = cell(edgeNo, 1);
anovaStats = cell(edgeNo, 1);
anovaComp = nan(edgeNo, condNo*(condNo-1)/2, condNo*(condNo-1)/2);

% clock for measuring elapsed time
startClock = tic;

%%%%%%%%%%%%%%%%%%%%%
%%%%% PARFOR!!! %%%%%
%%%%%%%%%%%%%%%%%%%%%
% loop through edges
parfor edgeIdx = 1:edgeNo

    % define data for edge
    edgeData = squeeze(connSim(edgeIdx, :, :));
    % mirror the upper triangular part of connSim matrix to the lower
    % half
    edgeData = triu(edgeData, 1) + triu(edgeData, 1)'; 

    % preallocate vars holding within- and across condition/stimulus similarities
    withinCondSim = nan((epochNo^2-epochNo)/2, condNo);
    acrossCondSim = nan(epochNo^2*(condNo-1), condNo);

    
    %% Comparisons separately for each condition / stimulus
    
    for condIdx = 1:condNo

        % within-cond epoch pairing similarities for given condition
        tmp = edgeData((condIdx-1)*epochNo+1:condIdx*epochNo, (condIdx-1)*epochNo+1:condIdx*epochNo); 
        withinCondSim(:, condIdx) = tmp(triu(true(epochNo), 1));

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
        acrossCondSim(:, condIdx) = reshape(tmpAll', [epochNo^2*(condNo-1), 1]);

        % store descriptives in result struct
        withinCondMean(edgeIdx, condIdx) = mean(withinCondSim(:, condIdx));
        withinCondSD(edgeIdx, condIdx) = std(withinCondSim(:, condIdx));
        withinCondMedian(edgeIdx, condIdx) = median(withinCondSim(:, condIdx));
        acrossCondMean(edgeIdx, condIdx) = mean(acrossCondSim(:, condIdx));
        acrossCondSD(edgeIdx, condIdx) = std(acrossCondSim(:, condIdx)); 
        acrossCondMedian(edgeIdx, condIdx) = median(acrossCondSim(:, condIdx));

        % permutation test across the two groups of connectivity similarity
        % values
        [pEst(edgeIdx, condIdx),... 
         realDiff(edgeIdx, condIdx),... 
         tmpPermDiff,... 
         cohend(edgeIdx, condIdx),...
         studentDiff(edgeIdx, condIdx)] = permTest(withinCondSim(:, condIdx),... 
                                             acrossCondSim(:, condIdx),... 
                                             permNo,... 
                                             permStat,...
                                             'studentized',...
                                             'silent');
         % get only mean and SD from the distribution of permuted stats                                
         permDiffMean(edgeIdx, condIdx) = mean(tmpPermDiff);
         permDiffSD(edgeIdx, condIdx) = std(tmpPermDiff);

    end
    
    
    %% Comparison between within- and across-condition pairings across all conditions / stimuli
    
    % store descriptives in result struct
    withinAllMean(edgeIdx) = mean(withinCondSim(:));
    withinAllSD(edgeIdx) = std(withinCondSim(:));
    withinAllMedian(edgeIdx) = median(withinCondSim(:));
    acrossAllMean(edgeIdx) = mean(acrossCondSim(:));
    acrossAllSD(edgeIdx) = std(acrossCondSim(:)); 
    acrossAllMedian(edgeIdx) = median(acrossCondSim(:));

    % permutation test on vectorized withinCondSim and acrossCondSim matrices
    [pEstAll(edgeIdx),... 
     realDiffAll(edgeIdx),... 
     tmpPermDiff,... 
     cohendAll(edgeIdx),...
     studentDiffAll(edgeIdx)] = permTest(withinCondSim(:),... 
                                          acrossCondSim(:),... 
                                          permNo,... 
                                          permStat,...
                                          'studentized',...
                                          'silent');
    % get only mean and SD from the distribution of permuted stats                                
    permDiffAllMean(edgeIdx) = mean(tmpPermDiff);
    permDiffAllSD(edgeIdx) = std(tmpPermDiff);          
    
    
    %% ANOVA for comparing within-cond similarities across conditions/stimuli
    
    [anovaP(edgeIdx), anovaTab{edgeIdx}, stats] = anova1(withinCondSim, [], 'off');
    anovaComp(edgeIdx, :, :) = multcompare(stats, 'ctype', 'lsd', 'display', 'off');
    anovaStats{edgeIdx} = stats;
    
    %% user message about progress
    elapsedTime = round(toc(startClock), 2);
    disp(['Done with edge no. ', num2str(edgeIdx), ', took ', num2str(elapsedTime), ' secs so far']);     
     
    
end  % parfor edgeIdx

% user message
disp([char(10), 'Finished comparing similarity within- versus across-condition epoch-pairings!']);


%% Collect resutls 

% all permutation test results into one struct
permRes = struct;
permRes.withinCondMean = withinCondMean;
permRes.withinCondSD = withinCondSD;
permRes.withinCondMedian = withinCondMedian;
permRes.acrossCondMean = acrossCondMean;
permRes.acrossCondSD = acrossCondSD;
permRes.acrossCondMedian = acrossCondMedian;
permRes.pEst = pEst;
permRes.realDiff = realDiff;
permRes.studentDiff = studentDiff;
permRes.permDiffMean = permDiffMean;
permRes.permDiffSD = permDiffSD;
permRes.cohend = cohend;
permRes.withinAllMean = withinAllMean;
permRes.withinAllSD = withinAllSD;
permRes.withinAllMedian = withinAllMedian;
permRes.acrossAllMean = acrossAllMean;
permRes.acrossAllSD = acrossAllSD;
permRes.acrossAllMedian = acrossAllMedian;
permRes.pEstAll = pEstAll;
permRes.realDiffAll = realDiffAll;
permRes.studentDiffAll = studentDiffAll;
permRes.permDiffAllMean = permDiffAllMean;
permRes.permDiffAllSD = permDiffAllSD;
permRes.cohendAll = cohendAll;

% all anova results into other struct
withinCondAnova = struct;
withinCondAnova.anovaP = anovaP;
withinCondAnova.anovaTab = anovaTab;
withinCondAnova.anovaStats = anovaStats;
withinCondAnova.anovaComp = anovaComp;

% user message
disp([char(10), 'Done with everything, returning results']);


return















