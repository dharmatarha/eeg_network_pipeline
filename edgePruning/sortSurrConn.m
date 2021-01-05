function sortSurrConn(selectedConnFile, surrDataDir, freq, method)
%% Function to assign surrogate-data derived statistics to connectivity data
%
% USAGE: sortSurrConn(selectedConnFile, surrDataDir, freq, method)
%
% In the resting-state EEG dataset, we calculate connectivity values for
% the data of each subject with "connectivityWrapperReal" then randomly 
% select non-overlapping epochs using "selectConnEpochs". Separately, we
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
% The output of the function is three variables containing the normal
% distribution parameters (mu and sigma) and the empirical p-value
% estimates.
%
%
%


%% Input checks

if nargin ~= 4
    error('Function sortSurrConn requires input args "selectedConnFile", "surrDataDir", "freq" and "method"!');
end
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

if strcmp(surrDataDir(end), '/')
    surrDataDir = surrDataDir(1:end-1);
end


%% Load connectivity data

tmp = load(selectedConnFile);
subjects = tmp.subjects;
epochIndices = tmp.epochIndices;
connData = tmp.connData;
[subNo, epochNo, roiNo, ~] = size(connData);


%% Loop through subjects

% result vars
surrNormalMu = nan(subNo, epochNo, roiNo, roiNo);
surrNormalSigma = surrNormalMu;
surrNormalP = surrNormalMu;
realConnP = surrNormalMu;
maskedConn = surrNormalMu;
criticalP = nan(subNo, epochNo);
survivalRate = nan(subNo, epochNo);

for subIdx = 1:subNo
    
    subClock = tic;
    
    % get current subject's surrogate connectivity data file
    subID = subjects{subIdx};
    subSurrFile = [surrDataDir, '/', subID, '_', freq, '_surrEdgeEstReal_', method, '.mat'];
    if ~exist(subSurrFile, 'file')
        error(['Cannot find surrogate data file for subject ', subID, '!']);
    end
    tmp = load(subSurrFile);
    
    % get epoch indices for current subject
    subEpochs = epochIndices{subIdx};
    if isempty(subEpochs)
        subEpochs = 1:epochNo;
    end

    % select mu, sigma and fit p-value for selected epochs
    tmpMu = permute(tmp.surrNormalMu, [3 1 2]);
    tmpSigma = permute(tmp.surrNormalSigma, [3 1 2]);
    tmpP = permute(tmp.surrNormalP, [3 1 2]);
    
    % !!! IMPORTANT !!!
    % "selectConnEpochs" epochIndices index epochs after every second was
    % discarded!
    tmpMu = tmpMu(1:2:end, :, :);
    tmpSigma = tmpSigma(1:2:end, :, :);
    tmpP = tmpP(1:2:end, :, :);
    
    % collect subject-level (truncated) normal params to group-level var
    surrNormalMu(subIdx, :, :, :) = tmpMu(subEpochs, :, :);
    surrNormalSigma(subIdx, :, :, :) = tmpSigma(subEpochs, :, :);
    surrNormalP(subIdx, :, :, :) = tmpP(subEpochs, :, :);
    
    % estimate p-values for real connectivity
    for epochIdx = 1:epochNo
        for roi1 = 1:roiNo
            for roi2 = 1:roiNo
                if roi1 < roi2
                    
                    % get p-value based on normal distribution of surrogate
                    % values, multiply by 2 for 2-sided test
                    normalP = normcdf(connData(subIdx, epochIdx, roi1, roi2),... 
                                        surrNormalMu(subIdx, epochIdx, roi1, roi2),... 
                                        surrNormalSigma(subIdx, epochIdx, roi1, roi2));
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
        epochConnData(epochConnData < criticalP(subIdx)) = 0;
        maskedConn(subIdx, epochIdx, :, :) = epochConnData;
        % survived edges
        edgesAboveCrit = epochConnData(:, :) >= criticalP(subIdx);
        survivalRate(subIdx, epochIdx) = sum(edgesAboveCrit(:), 'omitnan')/(roiNo*(roiNo-1)/2);
        
    end
    
    elapsedT = round(toc(subClock), 2);
    disp([char(10), 'Finished with subject ', subID, '. Took ', num2str(elapsedT), ' secs.']);
    
end









