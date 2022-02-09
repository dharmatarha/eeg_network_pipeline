
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
        
        if freqBandIndex == 1 && connMetricIndex == 1
            [subjectFiles, epochIndices, connData] = selectEpochs(dirName, ...
                filePattern, varname, epochDim, epochNo, epochMask);
            epochIndicesConstraint = epochIndices;
            subjectIDs = cell(subjectNo, 1);
            for subjectIndex = 1 : size(connData, 1)
                fileNameunderTest = subjectFiles{subjectIndex};
                subjectIDs{subjectIndex} = fileNameunderTest(numel([dirName, '/'])+1 : end-numel(['_' filePattern]));
            end
            save([dirName, '/', 'selectedEpochsSelectedSubjects.mat'], 'epochIndicesConstraint', 'subjectIDs');
        else
            [subjectFiles, epochIndices, connData] = selectEpochs(dirName, ...
                filePattern, varname, epochDim, subjectIDs, epochIndicesConstraint);
        end

        %% Sanity check

        % Check if dimensions are correct
        if size(subjectFiles, 1) ~= subjectNo || size(subjectFiles, 2) ~= 1
            error('Dimensions returned by epoch selection is incorrect!');
        end

        if size(epochIndices, 1) ~= subjectNo || size(epochIndices, 2) ~= 1
            error('Dimensions returned by epoch selection is incorrect!');
        end

        if size(connData, 1) ~= subjectNo || size(connData, 2) ~= epochNo || ...
                size(connData, 3) ~= roiNo || size(connData, 4) ~= roiNo
            error('Dimensions returned by epoch selection is incorrect!');
        end

        % Check if all subjects are unique
        uniqueElementsOfSubjectFiles = unique(subjectFiles);
        if size(uniqueElementsOfSubjectFiles, 1) ~= subjectNo
            error('Not all subjects returned by epoch selection are unique!')
        end

        % Check if all subjects are real
        selectedFileNames = cell(subjectNo, 1);
        for subjectIndex = 1 : size(connData, 1)
            fileNameunderTest = subjectFiles{subjectIndex};
            selectedFileNames{subjectIndex} = fileNameunderTest(numel([dirName, '/'])+1 : end);
        end

        if dirName(end) == '/'
            dirName = dirName(1:end-1);
        end
        filesInFolder = dir([dirName, '/*', filePattern, '*']); 
        fileNamesInFolder = cell(size(filesInFolder, 1), 1);
        for index = 1 : size(filesInFolder, 1)
            fileNamesInFolder{index} = filesInFolder(index).name;
        end

        for subjectIndex = 1 : subjectNo
            if ~any(contains(fileNamesInFolder, selectedFileNames{subjectIndex}))
                error('Not all subjects returned by epoch selection are real subjects!');
            end
        end

        subjects = cell(subjectNo, 1);
        for subjectIndex = 1 : subjectNo
            fileNameUnderTest = selectedFileNames{subjectIndex};
            subjects{subjectIndex} = fileNameUnderTest(1 : end-numel(['_', filePattern]));
        end

        % Check if all selected epochs are uniqe
        if freqBandIndex == 1 && connMetricIndex == 1
            for subjectIndex = 1 : subjectNo
                if numel(unique(epochIndices{subjectIndex})) ~= epochNo
                    error('Number of uniqe epochs selected is less than expected!');
                end
            end
        else
            for subjectIndex = 1 : subjectNo
                if epochIndices{subjectIndex} ~= epochIndicesConstraint{subjectIndex}
                    error('Selected epochs do not correspond to the selection of the first run!');
                end
            end
        end

        % Check if all epochs contain non-Nan values in their upper-triangular matrix
        % and NaN values in their lower-triangular matrix (including the diagonal)
        for subjectIndex = 1 : subjectNo
            for epochIndex = 1 : epochNo
                connDataUnderTest = squeeze(connData(subjectIndex, epochIndex, :, :));
                linTriuConnData = connDataUnderTest(triu(true(roiNo), 1));
                if any(isnan(linTriuConnData))
                    error('At least one NaN value is found where real data should stay!');
                end
                if any(linTriuConnData == 0)
                    error('At least one zero value is found where real data should stay!');
                end
                if ~isreal(linTriuConnData)
                    error('At least one zero value is found where real data should stay!');
                end
                linTrilConnData = connDataUnderTest(tril(true(roiNo), 0));
                if any(~isnan(linTrilConnData))
                    error('At least one non-NaN value is found where only NaN values should stay!');
                end
            end
        end

        % User message
        disp([char(10), 'Sanity check after epoch selection performed. No issues found.']);

        % Save variables to file
        connectivityFileName = [freqBandString, '_', connMetricString, '_group', '.mat'];
        save([dirName, '/', connectivityFileName], 'subjects', 'epochIndices', 'connData');
        disp([char(10), 'Connectivity file: ', connectivityFileName,...
            char(10), 'saved with variables: "subjects", "epochIndices", "connData"']);

        if ismember(connMetric, {'plv', 'iplv'})
            truncated = 'truncated';
        else
            truncated = 'nontruncated';
        end
        sortSurrConn([dirName, '/', connectivityFileName], dirName, frequencyBand, connMetric, truncated);
    end
end