
frequencyBands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
numberOfFrequencyBands = numel(frequencyBands);
connMetrics = {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'};
numberOfConnMetrics = numel(connMetrics);
varname = 'connRes';
epochNo = 75;
subjectNo = 200;
roiNo = 62;

dirNameBase = '../NAS502/EEG_resting_state/';
for freqBandIndex = 1 : numberOfFrequencyBands
    
    groupLevelConnVectors = nan(subjectNo*epochNo*roiNo*(roiNo-1)/2, numberOfConnMetrics);
    for connMetricIndex = 1 : numberOfConnMetrics
        
        freqBandString = frequencyBands{freqBandIndex};
        connMetricString = connMetrics{connMetricIndex};
        dirName = [dirNameBase freqBandString];
        fileName = [freqBandString, '_', connMetricString, '_group.mat'];
        frequencyBand = freqBandString;
        connMetric = connMetricString;
        
        load([dirName, '/', fileName], 'connData');
        connectivityTensor_unthr = connData;
        epochLevelConnVectors = nan(subjectNo, epochNo, roiNo*(roiNo-1)/2);
        for subjectIndex = 1 : subjectNo
            for epochIndex = 1 : epochNo
                epochLevelConnMatrix = squeeze(connectivityTensor_unthr(subjectIndex, epochIndex, :, :));
                epochLevelConnVector = epochLevelConnMatrix(triu(true(roiNo), 1));
                epochLevelConnVectors(subjectIndex, epochIndex, :) = epochLevelConnVector;
            end
        end
        groupLevelConnVector = epochLevelConnVectors(:);
        groupLevelConnVectors(:, connMetricIndex) = groupLevelConnVector;
    end % for connMetricIndex = 1 : numberOfConnMetrics
    
    matrixOfCorrelationValues = corr(groupLevelConnVectors);
    figure();
    imagesc(matrixOfCorrelationValues);
    xticks([1 2 3 4]);
    xticklabels({'PLV' 'iPLV' 'ampCorr' 'orthAmpCorr'});
    yticks([1 2 3 4]);
    yticklabels({'PLV' 'iPLV' 'ampCorr' 'orthAmpCorr'});
    title(['No avg. corr. ', freqBandString, ' band, unthr']);    
    colormap('jet');
    caxis([-0.8 1]);
    colorbar;
    set(gca, 'FontSize', 14);
    
end % for freqBandIndex = 1 : numberOfFrequencyBands
