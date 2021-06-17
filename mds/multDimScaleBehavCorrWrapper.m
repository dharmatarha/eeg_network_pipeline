%% Specify parameters
%
% isThresholdingApplied                       - 1 for thresholded data, 0 for unthresholded data
% connMeasure                                 - connectivity measure, one of {'plv', 'iplv'}
% graphDistanceMetric                         - graph distance metrci, one of {'adjacencySpectral', 'LaplacianSpectral'}
% numberOfIterations                          - number of iterations for surrogate data generation, between 1 and 10000
% numberOfDimensionsToConsider                - number of dimensions from multidimensional scaling to consider, between 1 and 5
% surrogateDataType                           - surrogate data type, one of {'preserving', 'existingEdgesRandom', 'edgesRandom'}
% pathName_behavData                          - path name for behavior data
% pathName_connData                           - path name for connectivity data
% fileName_behavData                          - file name for behavior data
% fileNamePrefix_connData_withThresholding    - file name prefix for thresholded connectivity data
% fileNamePrefix_connData_withoutThresholding - file name prefix for unthresholded connectivity data
%

isThresholdingApplied = 0;
connMeasure = 'plv';
graphDistanceMetric = 'adjacencySpectral';
numberOfIterations = 1000;
numberOfDimensionsToConsider = 5;
surrogateDataType = 'edgesRandom';
% pathName_behavData = 'D:\mult_dim_scale\';
pathName_behavData = '/home/adamb/Downloads/';
% pathName_connData = 'D:\mult_dim_scale\';
pathName_connData = '/media/adamb/bonczData/EEG_resting_state/alpha/';
fileName_behavData = 'Big_data_all_behav.xlsx';
fileNamePrefix_connData_withThresholding = 'surrConn_alpha_';
fileNamePrefix_connData_withoutThresholding = 'group_alpha_';

%
if ~ismember(isThresholdingApplied, 0:1)
    isThresholdingApplied = 1;
end
if ~ismember(connMeasure, {'plv', 'iplv'})
    connMeasure = 'plv';
end
if ~ismember(graphDistanceMetric, {'adjacencySpectral', 'LaplacianSpectral'})
    graphDistanceMetric = 'adjacencySpectral';
end
if ~ismember(numberOfIterations, 1:10000)
    numberOfIterations = 10000;
end
if ~ismember(numberOfDimensionsToConsider, 1:5)
    numberOfDimensionsToConsider = 5;
end
if ~ismember(surrogateDataType, {'preserving', 'existingEdgesRandom', 'edgesRandom'})
    surrogateDataType = 'preserving';
end


%% Read behavior data
%
% (1) Read behavior data from a table
% (2) Delete subject 139
% (3) Store behavior variables to vectors
%
behaviorTable = readtable([pathName_behavData fileName_behavData]);
if isThresholdingApplied
    behaviorTable(139, :) = [];
    numberOfSubjects = 199;
else
    numberOfSubjects = 200;
end
numberOfBehaviorVariables = 7;
behaviorVectors = nan(numberOfSubjects, numberOfBehaviorVariables);
behaviorVectors(:, 1) = behaviorTable.HIT;
behaviorVectors(:, 2) = behaviorTable.FA;
behaviorVectors(:, 3) = behaviorTable.DIST;
behaviorVectors(:, 4) = behaviorTable.dp;
behaviorVectors(:, 5) = behaviorTable.RT_avg;
behaviorVectors(:, 6) = behaviorTable.RT_med;
behaviorVectors(:, 7) = behaviorTable.Memory_acc;

disp('Loaded behavioral data');


%% Read connectivity data
%
% (1) Read connectivity data (thresholded or unthresholded) from file
% (2) Average epochs for unthresholded data (for thresholded data, it is already done)
% (3) Delete subject 139
%
if isThresholdingApplied
    fileName_connData = ([fileNamePrefix_connData_withThresholding connMeasure '.mat']);
    dataStructure = open([pathName_connData fileName_connData]);
    connectivityTensor = dataStructure.acrossEpochs.maskedConn;
else
    fileName_connData = ([fileNamePrefix_connData_withoutThresholding connMeasure '.mat']);
    dataStructure = open([pathName_connData fileName_connData]);
    connectivityTensor = dataStructure.connData;
    connectivityTensor = squeeze(mean(connectivityTensor, 2));
end
if isThresholdingApplied
    connectivityTensor(139, :, :, :) = [];
end

disp('Loaded connectivity data');


%% Surrogate data generation and behavior correlation computation

% Define tensors for the surrogate distance matrices and surrogate connectivity matrices
[subNo, roiNo, ~] = size(connectivityTensor);
surrogateDistanceTensor = nan(numberOfIterations, subNo, subNo);
surrogageCorrelationTensor = nan(numberOfIterations, numberOfBehaviorVariables, numberOfDimensionsToConsider);

randLoopClock = tic;
% Perform surrogate iterations the required times
parfor iterationIndex = 1 : numberOfIterations
    % Create surrogate connectiviry tensor (edge rewiring for each subject separately)
    connectivityTensor_randomized = nan(subNo, roiNo, roiNo);
    for subIdx = 1:subNo
        connectivityMatrixUnderTest = squeeze(connectivityTensor(subIdx, :, :));
        connectivityMatrixUnderTest = triu(connectivityMatrixUnderTest, 1) + triu(connectivityMatrixUnderTest, 1)';
        switch surrogateDataType
            case 'preserving'
                connectivityTensor_randomized(subIdx, :, :) = null_model_und_sign_mod(connectivityMatrixUnderTest);
            case 'existingEdgesRandom'
                connectivityTensor_randomized(subIdx, :, :) = edgeRandomization(connectivityMatrixUnderTest, true);
            case 'edgesRandom'
                connectivityTensor_randomized(subIdx, :, :) = edgeRandomization(connectivityMatrixUnderTest, false);
        end
    end
    % Create distance matrix (distance between all subjects) for surrogate data
    [distRes_randomized] = connDistanceTest_betweenSubject_epochAveraged(connectivityTensor_randomized, graphDistanceMetric, 'silent');
    distRes_randomized(isnan(distRes_randomized)) = 0;
    distRes_randomized = distRes_randomized + distRes_randomized';
    surrogateDistanceTensor(iterationIndex, :, :) = distRes_randomized;
    % Calculate subject coordinates by classical multidimensional scaling for surrogate data
    [Y_randomized, eigvals_randomized] = cmdscale(distRes_randomized);
    coordinateVectors_randomized = Y_randomized(:, 1:numberOfDimensionsToConsider);
    % Calculate coordinates between surrogate subject coordinates and behavior data
    correlationValues_randomized = corr(behaviorVectors, coordinateVectors_randomized, 'Type', 'Spearman');
    surrogageCorrelationTensor(iterationIndex, :, :) = correlationValues_randomized;
end
elapsedTime = toc(randLoopClock);
disp(['Elapsed time: ', num2str(round(elapsedTime, 4)), ' secs']);

% Create distance matrix (distance between all subjects) for real data
[distRes_real] = connDistanceTest_betweenSubject_epochAveraged(connectivityTensor, graphDistanceMetric, 'silent');
distRes_real(isnan(distRes_real)) = 0;
distRes_real = distRes_real + distRes_real';
% Calculate subject coordinates by classical multidimensional scaling 
[Y_real, eigvals_real] = cmdscale(distRes_real);
coordinateVectors_real = Y_real(:, 1:numberOfDimensionsToConsider);
% Calculate coordinates between subject coordinates and behavior data
correlationValues_real = corr(behaviorVectors, coordinateVectors_real, 'Type', 'Spearman');


