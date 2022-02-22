
load('/home/peternagy/NAS502/EEG_resting_state/theta/selectedEpochsSelectedSubjects.mat', 'subjectIDs');
frequencyBands = {'theta', 'alpha', 'beta', 'gamma', 'delta'};
numberOfFrequencyBands = numel(frequencyBands);
connMetrics = {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'};
numberOfConnMetrics = numel(connMetrics);
varname = 'connRes';
epochDim = 1;
epochNo = 75;
epochMask = zeros(1,300); epochMask(1:2:end) = 1;
subjectNo = 200;
roiNo = 62;

dirNameBase = '../NAS502/EEG_resting_state/';
for freqBandIndex = 1 : numberOfFrequencyBands
    for connMetricIndex = 1 : numberOfConnMetrics
        
        freqBandString = frequencyBands{freqBandIndex};
        connMetricString = connMetrics{connMetricIndex};
        dirName = [dirNameBase freqBandString];
        filePattern = [freqBandString, '_', connMetricString, '.mat'];
        frequencyBand = freqBandString;
        connMetric = connMetricString;
        
        [selectedFiles, selectedEpochs, distSimResult] = epochDistSim_subject(dirName, ...
            filePattern, varname, epochDim, epochNo, subjectIDs, epochMask);

        %% Sanity check

        % Check if dimensions are correct
        if numel(size(distSimResult)) ~= 4 || size(distSimResult, 1) ~= subjectNo ||...
                size(distSimResult, 2) ~= 2 || size(distSimResult, 3) ~= epochNo ||...
                size(distSimResult, 4) ~= epochNo
            error('Dimensions returned by epochDistSim_subject are incorrect!');
        end

        % User message
        disp([char(10), 'Sanity check after epoch distance and similarity calculation performed. No issues found.']);

        % Save variables to file
        distSimFileName = [freqBandString, '_', connMetricString, '_epochDistSim', '.mat'];
        save([dirName, '/', distSimFileName], 'selectedFiles', 'selectedEpochs', 'distSimResult');
        disp([char(10), 'Epoch distance-similarity file: ', distSimFileName,...
            char(10), 'saved with the variable: "distSimResult"']);
        
    end
end