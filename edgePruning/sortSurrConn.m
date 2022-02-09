function sortSurrConn(selectedConnFile, surrDataDir, freq, method, truncated)
%% Function to assign surrogate-data derived statistics to connectivity data
%
% USAGE: sortSurrConn(selectedConnFile, surrDataDir, freq, method, truncated=[])
%
% In the resting-state EEG dataset, we calculate connectivity values for
% the data of each subject with "connectivityWrapperReal" then randomly 
% select non-overlapping epochs using "selectEpochs". Separately, we
% also generate surrogate (phase-scrambled) data and estimate probability
% distribution parameters for the connectivity values in the surrogate 
% data (using "surrEdgeEstimationReal"). 
%
% This function assigns the surrogate connectivity distribution
% parameters to the selected epochs of the real connectivity data. The
% input "selectedConnFile" is the path to the file generated by
% "selectConnEpochs" (e.g. "group_alpha_ampCorr.mat")  and contains the 
% full connectivity dataset in var "connData", the correspondonding 
% subject IDs in "subjects" and the indices of the selected epochs in 
% "epochIndices". Input arg surrDataDir points to the directory containing 
% the "surrEdgeEstimationReal" output files 
% (e.g. "l3_s01_alpha_surrEdgeEstReal_ampCorr.mat"). 
%
% The outputs are saved out into a file
% 'surrConn_FREQUENCYBAND_METHOD.mat'.
%
% Mandatory inputs:
% selectedConnFile      - Char array, path to group-level connectivity data
%                       file (e.g. 'group_alpha_orthAmpCorr.mat'). Contains
%                       the following variables:
%                       subjects:       cell array with subject ids
%                       epochIndices:   cell array with indeices of
%                                       selected epochs for each subject
%                       connData:       numeric array containing
%                                       connectivity data, sized [subjects
%                                       X epochs X rois X rois]
% surrDataDir           - Char array, folder containing the surrogate
%                       connectivity files for each subject (e.g.
%                       'l3_s01_alpha_surrEdgeEstReal_orthAmpCorr.mat'). 
% freq                  - Char array, one of {'delta', 'theta', 'alpha', 
%                       'beta', 'gamma'}. Frequency band which is reflected
%                       in surrogate connectivity file names.
% method                - Char array, one of {'plv', 'iplv', 'pli', 
%                       'ampCorr', 'orthAmpCorr'}. Connectivity method,
%                       reflected in surrogate connectivity file names.
%
% Optional inputs:
% truncated             - Char array, one of {'truncated', 'nontruncated'}.
%                       Defines if the normal distributions derived from
%                       the surrogate data should be truncated to [0 1]. If
%                       left empty, the function tries to get this
%                       information from the surrogate data files (from var
%                       "truncated"). Defaults to empty.
%
% Outputs:
%
% The output is saved into a file named "surrConn_FREQUENCYBAND_METHOD.mat".
% It contains the following variables:
%
% surrNormalMu          - 4D numeric array, contains the "mu"s of the
%                       normal distributions fitted to the surrogate 
%                       connectivity values. Sized [subjects X epochs X
%                       rois X rois].
% surrNormalSigma       - 4D numeric array, contains the "sigma"s of the
%                       normal distributions fitted to the surrogate 
%                       connectivity values. Sized [subjects X epochs X
%                       rois X rois].
% surrNormalP           - 4D numeric array, contains the probabilities 
%                       from the Kolmogorov-Smirnov tests on the normal 
%                       fits of the surrogate connectivity values. Sized 
%                       [subjects X epochs X rois X rois].
% realConnP             - 4D numeric array, contains the probabilities 
%                       of the real connectivity values based on the normal
%                       distributions of the surrogate values. Sized 
%                       [subjects X epochs X rois X rois].
% maskedConn            - 4D numeric array, connectivity values after
%                       thresholding based on surrogate connectivity
%                       values.
% criticalP             - 2D numeric array, contains the critical p values
%                       for each epoch from FDR (q = .05), used for 
%                       thresholding. Sized [subjects X epochs].
% survivalRate          - 2D numeric array, contains the ratio of edges
%                       surviving the FDR-based thresholding in each epoch.
%                       Sized [subjects X epochs].
%


%% Input checks

% check no. of inputs
if ~ismember(nargin, 4:5)
    error('Function sortSurrConn requires input args "selectedConnFile", "surrDataDir", "freq" and "method", while input arg "truncated" is optional!');
end
% check each mandatory input
if ~exist(selectedConnFile, 'file')
    error(['Input arg "selectedConnFile" should point to a file containing ',...
        'the connectivity data from the randomly selected epochs (output of "selectConnEpochs")!']);
end
if ~exist(surrDataDir, 'dir')
    error(['Input arg "surrDataDir" should point to a directory containing ',...
        'the surrogate connectivity data (outputs of "surrEdgeEstimationReal")!']);
end
if ~ismember(freq, {'alpha', 'beta', 'gamma', 'delta', 'theta'})
    error('Input arg "freq" should be one of {"alpha", "beta", "gamma", "delta", "theta"}!');
end
if ~ismember(method, {'plv', 'iplv', 'pli', 'ampCorr', 'orthAmpCorr'})
    error('Input arg "freq" should be one of {"plv", "iplv", "pli", "ampCorr", "orthAmpCorr"}!');
end
% check optional input
if nargin == 4 || isempty(truncated)
    truncated = 'filebased';
elseif nargin == 5
    if ~ismember(truncated, {'truncated', 'nontruncated'})
        error('Input arg "truncated" should be one of {"truncated", "nontruncated"}!');
    end
end


% user message
disp([char(10), 'Called function sortSurrConn with inputs: ',...
    char(10), 'Connectivity file: ', selectedConnFile,...
    char(10), 'Surrogate data folder: ', surrDataDir,...
    char(10), 'Frequency band: ', freq,...
    char(10), 'Connectivity method: ', method,...
    char(10), 'Truncate normals: ', truncated]);


%% Settings, params

% strip ending '/' from folder path
if strcmp(surrDataDir(end), '/')
    surrDataDir = surrDataDir(1:end-1);
end


%% Load connectivity data

% connectivity data file should contain cell arrays for subject IDs and
% epoch indices, plus a numeric array holding all connectivity data
tmp = load(selectedConnFile);
subjects = tmp.subjects;  % subjects X 1 cell array, each cell contains a char array as subject ID
epochIndices = tmp.epochIndices;  % subjects X 1 cell array, each cell contains a vector of epoch indices
connData = tmp.connData;  % 4D array, subects X epochs X rois X rois
[subNo, epochNo, roiNo, ~] = size(connData);


%% Loop through subjects

% preallocate result vars
surrNormalMu = nan(subNo, epochNo, roiNo, roiNo);
surrNormalSigma = surrNormalMu;
surrNormalP = surrNormalMu;
realConnP = surrNormalMu;
maskedConn = surrNormalMu;
criticalP = nan(subNo, epochNo);
survivalRate = nan(subNo, epochNo);
% preallocate a struct for holding the group-level thresholding results
acrossEpochs = struct;
acrossEpochs.surrNormalMu = nan(subNo, roiNo, roiNo);
acrossEpochs.surrNormalSigma = acrossEpochs.surrNormalMu;
acrossEpochs.realConnP = acrossEpochs.surrNormalMu;
acrossEpochs.maskedConn = acrossEpochs.surrNormalMu;
acrossEpochs.criticalP = nan(subNo, 1);
acrossEpochs.survivalRate = nan(subNo, 1);


for subIdx = 1:subNo
    
    
    %% Load surrogate data, preparations
    
    subClock = tic;
    
    % get current subject's surrogate connectivity data file
    subID = subjects{subIdx};
    subSurrFile = dir([surrDataDir, '/', subID, '_', freq, '_surrEdgeEstReal_*', method, '*.mat']);
    if numel(subSurrFile) ~= 1
        error(['There are none or too many surrogate data file(s) for subject ', subID, '!']);
    end
    subSurrPath = fullfile(subSurrFile.folder, subSurrFile.name);
    tmp = load(subSurrPath);
    
    % check if it is a file containing multiple surrogates, get method
    % index, if necessary
    methodFlag = false;
    if iscell(tmp.method) && numel(tmp.method) > 1
        methodFlag = true;
        methodIdx = find(strcmp(method, tmp.method));
        if isempty(methodIdx)
            error(['Cannot find requested method in surrogate connectivity file for subject ', subID, '!']);
        end
    end
    
    % get epoch indices for current subject
    subEpochs = epochIndices{subIdx};
    if isempty(subEpochs)
        subEpochs = 1:epochNo;
    end

    % select mu, sigma and fit p-value for selected epochs
    % change the order of their dimensions first as they are rois X rois X
    % epochs
    if methodFlag
        tmpMu = permute(squeeze(tmp.surrNormalMu(methodIdx, :, :, :)), [3 1 2]);
        tmpSigma = permute(squeeze(tmp.surrNormalSigma(methodIdx, :, :, :)), [3 1 2]);
        tmpP = permute(squeeze(tmp.surrNormalP(methodIdx, :, :, :)), [3 1 2]);        
    else
        tmpMu = permute(tmp.surrNormalMu, [3 1 2]);
        tmpSigma = permute(tmp.surrNormalSigma, [3 1 2]);
        tmpP = permute(tmp.surrNormalP, [3 1 2]);
    end
    
    % collect subject-level (truncated) normal params to group-level var
    surrNormalMu(subIdx, :, :, :) = tmpMu(subEpochs, :, :);
    surrNormalSigma(subIdx, :, :, :) = tmpSigma(subEpochs, :, :);
    surrNormalP(subIdx, :, :, :) = tmpP(subEpochs, :, :);
    
    % check for information about truncation: 
    % truncation info should be a logical var "truncated"
    if strcmp(truncated, 'filebased')
        if methodFlag
            if tmp.truncated(methodIdx)
                subTruncated = 'truncated';
            else
                subTruncated = 'nontruncated';
            end
        else
            if tmp.truncated
                subTruncated = 'truncated';
            else
                subTruncated = 'nontruncated';
            end
        end
    else
        subTruncated = truncated;
    end
    
    
    %% Estimate p-values for real connectivity in each epoch separately
    % Logic:
    % - go through each edge, get the p-value of the corresponding
    % connectivity value from the normal distribution of surrogate data
    % - FDR across all edges in epoch, q = 0.05
    % - Thresholding
    
    for epochIdx = 1:epochNo
        for roi1 = 1:roiNo
            for roi2 = 1:roiNo
                % work only in upper triangle of connectivity matrix
                % (assume symmetric values)
                if roi1 < roi2
                    
                    % Case of using simple normal distribution
                    if strcmp(subTruncated, 'nontruncated')
                    
                        % get p-value based on normal distribution of surrogate
                        % values
                        normalP = normcdf(connData(subIdx, epochIdx, roi1, roi2),... 
                                            surrNormalMu(subIdx, epochIdx, roi1, roi2),... 
                                            surrNormalSigma(subIdx, epochIdx, roi1, roi2));
                        
                    % Case of using truncated normal distribution
                    elseif strcmp(subTruncated, 'truncated')
                        
                        % get a truncated normal distribution first
                        pd = makedist('normal', 'mu', surrNormalMu(subIdx, epochIdx, roi1, roi2), 'sigma', surrNormalSigma(subIdx, epochIdx, roi1, roi2));  % returns a probability distribution object
                        
                        % the following is a necessity for 
                        pd = truncate(pd, 0, 1);
                        
                        % get p-value from the cumulative version of the
                        % distribution function
                        normalP = pd.cdf(connData(subIdx, epochIdx, roi1, roi2));
                        
                    end
                    
                    % convert for a 2-two-sided test, apply correction as
                    % well
                    if normalP > 0.5
                        normalP = 1-normalP;
                    end                                         
                    realConnP(subIdx, epochIdx, roi1, roi2) = normalP*2;
                    
                end
            end
        end
        
        % FDR correction on epoch-level
        realPs = realConnP(subIdx, epochIdx, :, :); 
        realPs = realPs(:); realPs(isnan(realPs)) = [];
        [~, criticalP(subIdx, epochIdx)] = fdr(realPs, 0.05, 'bh'); 
        
        % create masked connectivity tensor
        epochConnData = squeeze(connData(subIdx, epochIdx, :, :));
        epochP = squeeze(realConnP(subIdx, epochIdx, :, :));
        epochConnData(epochP > criticalP(subIdx, epochIdx)) = 0;
        maskedConn(subIdx, epochIdx, :, :) = epochConnData;
        % get ratio of edges below the threshold (edges surviving the
        % pruning)
        edgesBelowCrit = epochP <= criticalP(subIdx, epochIdx);
        survivalRate(subIdx, epochIdx) = sum(edgesBelowCrit(:), 'omitnan')/(roiNo*(roiNo-1)/2);
        
    end  % for epochIdx
    
    
    %% Estimate p-values for average connectivity across epochs
    % Logic:
    % - get mean connectivity across epochs
    % - for each edge, compare the real connectivity value to the normal
    % distribution derived from the per-epoch surrogate normals
    % - FDR across edges, q = 0.05
    % - Thresholding
    
    % mean connectivity across epochs for current subject
    subConnData = squeeze(mean(connData(subIdx, :, :, :), 2));
    
    % loops across ROIs (defining edges)
    for roi1 = 1:roiNo
        for roi2 = 1:roiNo
            % work only in upper triangle of connectivity matrix
            % (assume symmetric values)
            if roi1 < roi2            
                
                % Derive normal distribution of surrogate data for the
                % "group of epochs"
                % NOTE: Linear combinations of normal distributions are
                % just linear combinations of parameters
                mu = mean(squeeze(surrNormalMu(subIdx, :, roi1, roi2)), 'omitnan');  % mu is simple the mean of all mu values
                sigma = mean(squeeze(surrNormalSigma(subIdx, :, roi1, roi2)), 'omitnan')/sqrt(epochNo);  % sigma is the mean sigma divided by sqrt(n)
                
                % Case of using simple normal distribution
                if strcmp(subTruncated, 'nontruncated')
                    normalP = normcdf(subConnData(roi1, roi2), mu, sigma);
                    
                % Case of using truncated normal distribution
                elseif strcmp(subTruncated, 'truncated')
                    % get a truncated normal distribution first
                    pd = makedist('normal', 'mu', mu, 'sigma', sigma);  % returns a probability distribution object
                    pd = truncate(pd, 0, 1);
                    % get p-value from the cumulative version of the
                    % distribution function
                    normalP = pd.cdf(subConnData(roi1, roi2));  
                end
                
                % convert for a 2-two-sided test, apply correction as well
                if normalP > 0.5
                    normalP = 1-normalP;
                end                                         
                acrossEpochs.realConnP(subIdx, roi1, roi2) = normalP*2;
                
                % store normal params
                acrossEpochs.surrNormalMu = mu;
                acrossEpochs.surrNormalSigma = sigma;
                
            end  % if
        end  % for roi1
    end  % for roi2
    
    % FDR correction on average connectivity matrix
    realPs = acrossEpochs.realConnP(subIdx, :, :); 
    realPs = realPs(:); realPs(isnan(realPs)) = [];
    [~, acrossEpochs.criticalP(subIdx)] = fdr(realPs, 0.05, 'bh');     
    
    % create masked connectivity tensor
    pMask = squeeze(acrossEpochs.realConnP(subIdx, :, :)) <= acrossEpochs.criticalP(subIdx);
    tmpData = subConnData;
    tmpData(~pMask) = 0;  % values not surviving the threshold are set to zero
    acrossEpochs.maskedConn(subIdx, :, :) = tmpData;
    
    % get rate of edges surviving pruning
    acrossEpochs.survivalRate(subIdx) = sum(pMask(:), 'omitnan')/(roiNo*(roiNo-1)/2);
    
    % display elapsed time 
    elapsedT = round(toc(subClock), 2);
    disp([char(10), 'Finished with subject ', subID, '. Took ', num2str(elapsedT), ' secs.']);
    
    
end  % for subIdx


%% Save out

saveF = [surrDataDir, '/', freq, '_', method, '_surrConn','.mat'];
save(saveF, 'surrNormalMu', 'surrNormalSigma', 'surrNormalP',... 
    'realConnP', 'maskedConn', 'criticalP', 'survivalRate',... 
    'acrossEpochs', 'truncated');






