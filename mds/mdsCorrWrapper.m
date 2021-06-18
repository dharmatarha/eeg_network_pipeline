function mdsCorrWrapper(dirName, connMeasure, thresholding, graphDistMetric, surrType, surrNo)

%% Wrapper function to correlate behavioral measures with surrogate (randomized) connectivity matrix MDS dimensions in the RS data set
%
% USAGE: mdsCorrWrapper(dirName=pwd, 
%                       connMeasure='plv', 
%                       thresholding='unthr', 
%                       grapDistMetric='adjacencySpectral', 
%                       surrType='edgesRandom', 
%                       surrNo=10)
%
% Workflow: 
% (1) Loading behavioral data and appropriate version of connectivity data
% (controlled by args "dirname", "connMeasure", "thresholding")
% (2) In a loop, generate surrogate data by randomizing connectivity 
% matrices ("surrType", "surrNo")
% (3) Estimate distances across surrogate connectivity matrices (surrogate
% distance matrix across subjects) ("graphDistMetric")
% (3) Classical multidimensional scaling (MDS) is calculated on surrogate
% connectivity data
% (4) Behavioral measures are correlated with the first few dimensions from
% MDS for each surrogate data set
%
% The results of the correlations and other potentially useful vars are 
% saved out in "dirName" to a file named:
% RsBehavCorr_CONNMEASURE_THRESHOLDING_GRAPHDISTMETRIC_SURRTYPE.mat
%
% NOTE THAT BY DEFAULT THE MAIN LOOP IS WITH PARFOR!
% CHANGE IT TO SIMPLE FOR LOOP IF NECESSARY!
%
% Optional inputs:
% dirName           - Char array, path to folder where the behavioral and
%               connectivity data files are located. Also the location for
%               results file. Defaults to current working directory.
% connMeasure       - Char array, one of 
%               {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'}. Connectivity
%               measure. The connectivity data file for the selected
%               connectivity measure is loaded for surrogate data
%               generation. Defaults to 'plv'.
% thresholding      - Char array, one of {'thr', 'unthr'}. Flag for loading
%               thresholded or unthresholded data. Defaults to 'unthr'
%               (unthresholded data).
% graphDistMetric   - Char array, one of 
%               {'adjacencySpectral', 'LaplacianSpectral'}. Distance
%               measure to be used on surrogate connectivity data for
%               generating the surrogate distance matric for MDS. Defaults
%               to 'adjacencySpectral'. 
% surrType          - Char array, one of 
%               {'edgesRandom', 'existingEdgesRandom', 'preserving'}.
%               Connectivity matrix randomization method, corresponding to
%               (1) weight randomization across all potential edges; (2)
%               weight randomization across existing edges; and (3)
%               randomization preserving weight, degree and strength
%               distributions, respectively. Defaults to 'edgesRandom'. 
% surrNo            - Numeric value in range 1:10^5. Number of surrogates
%               to generate. Defaults to 10 (testing purposes).
%
%
% Input args are simply positional. Any arg left empty is assinged its
% default value. 
%


%% Simplified input arg handling

% surrNo
if nargin < 6 || isempty(surrNo)
    surrNo = 10;
elseif ~isnumeric(surrNo) || ~ismember(surrNo, 1:10^5)
    error('Input arg "surrNo" should be numeric value in range 1:10^5!');
end
% surrType
if nargin < 5 || isempty(surrType)
    surrType = 'edgesRandom';
elseif ~ischar(surrType) || ~ismember(surrType, {'edgesRandom', 'existingEdgesRandom', 'preserving'})
    error('Input arg "surrType" should be one of {"edgesRandom", "existingEdgesRandom", "preserving"}!');
end
% graphDistMetric
if nargin < 4 || isempty(graphDistMetric)
    graphDistMetric = 'adjacencySpectral';
elseif ~ischar(graphDistMetric) || ~ismember(graphDistMetric, {'adjacencySpectral', 'LaplacianSpectral'})
    error('Input arg "graphDistMetric" should be one of {"adjacencySpectral", "LaplacianSpectral"}!');
end
% thresholding
if nargin < 3 || isempty(thresholding)
    thresholding = 'thr';
elseif ~ischar(thresholding) || ~ismember(thresholding, {'thr', 'unthr'})
    error('Input arg "thresholding" should be one of {"thr", "unthr"}!');
end
% connMeasure
if nargin < 2 || isempty(connMeasure)
    connMeasure = 'plv';
elseif ~ischar(connMeasure) || ~ismember(connMeasure, {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'})
    error('Input arg "connMeasure" should be one of {"plv", "iplv", "ampCorr", "orthAmpCorr"}!');
end
% dirName
if nargin < 1 || isempty(dirName)
    dirName = pwd;
elseif ~ischar(dirName) || ~exist(dirName, 'dir')
    error('Input arg "dirName" should be a valid path to a folder!');
end


%% Additional settings

%%%%%%%%%%%%%%% HARDCODED PARAMS !!! %%%%%%%%%%%%%%%%%%%%%%
mdsDimNo = 5;  % number of MDS dimensions to consider in correlations
fileName_behavData = 'Big_data_all_behav.xlsx';  % file name of behavioral data
freq = 'alpha';
fileNameStartThr = ['surrConn_', freq, '_'];  % file name start for thresholded connectivity data
fileNameStartUnthr = ['group_', freq, '_'];  % file name start for unthresholded connectivity data
behavVarNo = 7;  % number of behavioral variables to work with
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get boolean flag from thresholding, keep char array form as well for
% results file
thresholdingText = thresholding;
if strcmp(thresholding, 'unthr')
    thresholding = false;
else
    thresholding = true;
end


% user message
disp([char(10), 'Called mdsCorrWrapper with input args: ',...
    char(10), 'Connectivity measure: ', connMeasure,...
    char(10), 'Thresholding: ', num2str(thresholding),...
    char(10), 'Distance metric: ', graphDistMetric,...
    char(10), 'Randomization type: ', surrType,...
    char(10), 'Number of surrogate data sets: ', num2str(surrNo)]);



%% Read behavioral data

% load into table
behaviorTable = readtable([dirName fileName_behavData]);
behavSubNo = size(behaviorTable, 1);  % number of subjects based on behavioral data

% collect behavioral data into numeric matrix
behaviorVectors = nan(behavSubNo, behavVarNo);
behaviorVectors(:, 1) = behaviorTable.HIT;
behaviorVectors(:, 2) = behaviorTable.FA;
behaviorVectors(:, 3) = behaviorTable.DIST;
behaviorVectors(:, 4) = behaviorTable.dp;
behaviorVectors(:, 5) = behaviorTable.RT_avg;
behaviorVectors(:, 6) = behaviorTable.RT_med;
behaviorVectors(:, 7) = behaviorTable.Memory_acc;

disp('Loaded behavioral data');



%% Read connectivity data

% (1) Read connectivity data (thresholded or unthresholded) from file
% (2) Average epochs for unthresholded data (for thresholded data, it is already done)
% (3) Check for connectivity matrices filled only with zeros

if thresholding
    fileName_connData = fullfile(dirName, [fileNameStartThr, connMeasure, '.mat']);
    dataStructure = open(fileName_connData);
    connTensor = dataStructure.acrossEpochs.maskedConn;
else
    fileName_connData = fullfile(dirName, [fileNameStartUnthr, connMeasure, '.mat']);
    dataStructure = open(fileName_connData);
    connTensor = dataStructure.connData;
    connTensor = squeeze(mean(connTensor, 2));
end

% number of subjects in connectivity data
[connSubNo, roiNo, ~] = size(connTensor);
% Sanity check
if behavSubNo ~= connSubNo
    error('Subject numbers for behavioral data and connectivity data do not match!');
end

% change NaN values to zeros
connTensor(isnan(connTensor)) = 0;

disp('Loaded connectivity data');



%% Handle cases where there is missing data / no edge survived thresholding

% NaN values in behavioral data are handled at correlation estimation with
% param "rows" set to "complete" (only data without missing NaNs is taken
% into account).
% However, for connectivity data we might have thresholded data without
% surviving edges:

% check for connectivity matrices filled with zeros
zeroConns = false(connSubNo, 1);
for subIdx = 1:connSubNo
    if nnz(connTensor(subIdx, :, :)) == 0
        zeroConns(subIdx, 1) = true;
    end
end

% if there is any, report it and remove both connectivity and behavioral
% data
if sum(zeroConns) ~= 0
    % remove corresponding data
    behaviorVectors(zeroConns, :) = [];
    connTensor(zeroConns, :, :) = [];
    % list affected subjects
    removedSubs = find(zeroConns);
    % new subject number after remove
    subNo = connSubNo - sum(zeroConns);
    % report to user
    disp([char(10), 'Found ', num2str(sum(zeroConns)), ' empty connectivity matrix(ces), for subjects: ']);
    disp(removedSubs);
    disp(['Corresponding behavioral and connectivity data are removed, new sample size (N) is ', num2str(subNo)]);
else
    subNo = connSubNo;
    removedSubs = [];
end



%% Surrogate data generation and behavior correlation computation

% Define tensors for the surrogate distance matrices and surrogate connectivity matrices
surrDistTensor = nan(surrNo, subNo, subNo);
surrCorrTensor = nan(surrNo, behavVarNo, mdsDimNo);

% Clock to measure the elapsed time for the main loop
randLoopClock = tic;

% Surrogate generation loop
for surrIdx = 1:surrNo
% parfor surrIdx = 1:surrNo
    
    % Create surrogate connectivity tensor (edge rewiring for each subject separately)
    connTensorRand = nan(subNo, roiNo, roiNo);
    for subIdx = 1:subNo
        subConnMatrix = squeeze(connTensor(subIdx, :, :));
        subConnMatrix = triu(subConnMatrix, 1) + triu(subConnMatrix, 1)';
        switch surrType
            case 'preserving'
                connTensorRand(subIdx, :, :) = null_model_und_sign(subConnMatrix);
            case 'existingEdgesRandom'
                connTensorRand(subIdx, :, :) = edgeRandomization(subConnMatrix, true);
            case 'edgesRandom'
                connTensorRand(subIdx, :, :) = edgeRandomization(subConnMatrix, false);
        end
    end
    
    % Create distance matrix (distance between all subjects) for surrogate data
    surrDistMatrix = connDistanceTest_betweenSubject_epochAveraged(connTensorRand, graphDistMetric, 'silent');
    surrDistMatrix(isnan(surrDistMatrix)) = 0;
    surrDistMatrix = surrDistMatrix + surrDistMatrix';
    surrDistTensor(surrIdx, :, :) = surrDistMatrix;
    
    % Calculate subject coordinates by classical multidimensional scaling for surrogate data
    [Y_surr, ~] = cmdscale(surrDistMatrix);
    surrCoordVectors = Y_surr(:, 1:mdsDimNo);
    
    % Calculate coordinates between surrogate subject coordinates and behavior data
    surrCorrValues = corr(behaviorVectors, surrCoordVectors, 'type', 'Spearman', 'rows', 'complete');  % NaN handling is set
    surrCorrTensor(surrIdx, :, :) = surrCorrValues;
    
end  % parfor/for surrIdx

% Report elapsed time for the main loop
elapsedTime = toc(randLoopClock);
disp(['Elapsed time for surrogate generator loop: ', num2str(round(elapsedTime, 4)), ' secs']);



%% Get correlations between real connectivity data MDS dimensions and behavioral data

% Create distance matrix (distance between all subjects) for real data
realDistMatrix = connDistanceTest_betweenSubject_epochAveraged(connTensor, graphDistMetric, 'silent');
realDistMatrix(isnan(realDistMatrix)) = 0;
realDistMatrix = realDistMatrix + realDistMatrix';

% Calculate subject coordinates by classical multidimensional scaling 
[Y_real, eigvals_real] = cmdscale(realDistMatrix);
realCoordVectors = Y_real(:, 1:mdsDimNo);

% Calculate coordinates between subject coordinates and behavior data
realCorrValues = corr(behaviorVectors, realCoordVectors, 'type', 'Spearman', 'rows', 'complete');  % NaN handling is set



%% Compare surrogate data distance matrices, test their similarity

% We randomly select surrogate distance matrices and correlate them,
% "surrNo" times.

% var holding correlation coeffs between randomly selected surrogate
% distance matrices
distMatrixCorrs = nan(surrNo, 1); 
% loop for comparisons (correlations) across two randomly selected surrogate distance
% matrices
for i=1:surrNo
    tmp=randi(surrNo, 2, 1);  % get two random indices for the surrogate distance matrices 
    d1=squeeze(surrDistTensor(tmp(1), :, :)); 
    d2=squeeze(surrDistTensor(tmp(2), :, :));
    distMatrixCorrs(i) = corr(d1(:), d2(:)); 
end



%% Compare surrogate data distance matrices to real distance matrix

% var holding correlation coeffs between surrogate data distance matrices
% and the real distance matrix
surrRealDistCorrs = nan(surrNo, 1); 
% loop for comparisons (correlations) across two randomly selected surrogate distance
% matrices
for i=1:surrNo
    tmp = squeeze(surrDistTensor(i, :, :));
    surrRealDistCorrs(i) = corr(realDistMatrix(:), tmp(:)); 
end



%% Get p-values comparing surrogate and real correlations

% Simply call "estimatedP" for a two-tailed test,
% in a loop for each MDS dimension
pValues = nan(behavVarNo, mdsDimNo);
for dimNo = 1:mdsDimNo
    tmpRealValues = realCorrValues(:, dimNo);
    tmpSurrValues = surrCorrTensor(:, :, dimNo);
    pEst = estimatedP(tmpRealValues, tmpSurrValues, 2);
    pValues(:, dimNo) = pEst(:, 1);
end



%% Save out main results

saveFile = [dirName, 'RsBehavCorr_', connMeasure, '_',... 
    thresholdingText, '_', graphDistMetric, '_', surrType, '.mat'];
save(saveFile, 'behaviorVectors', 'realDistMatrix', 'Y_real', 'eigvals_real',...
    'realCoordVectors', 'realCorrValues', 'distMatrixCorrs',...
    'surrRealDistCorrs', 'surrCorrTensor', 'pValues', 'graphDistMetric',...
    'connMeasure', 'surrType', 'surrNo', 'thresholding', 'removedSubs');

return



