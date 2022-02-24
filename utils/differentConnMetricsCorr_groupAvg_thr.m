
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
    
    groupLevelAvgConnVectors = nan(roiNo*(roiNo-1)/2, numberOfConnMetrics);
    for connMetricIndex = 1 : numberOfConnMetrics
        
        freqBandString = frequencyBands{freqBandIndex};
        connMetricString = connMetrics{connMetricIndex};
        dirName = [dirNameBase freqBandString];
        fileName = [freqBandString, '_', connMetricString, '_surrConn.mat'];
        frequencyBand = freqBandString;
        connMetric = connMetricString;
        
        load([dirName, '/', fileName], 'acrossEpochs');
        connectivityTensor_thr = acrossEpochs.maskedConn;
        groupLevelAvgConnMatrix = squeeze(mean(connectivityTensor_thr, 1));
        groupLevelAvgConnVector = groupLevelAvgConnMatrix(triu(true(roiNo), 1));
        groupLevelAvgConnVectors(:, connMetricIndex) = groupLevelAvgConnVector;
    end % for connMetricIndex = 1 : numberOfConnMetrics
    
    matrixOfCorrelationValues = corr(groupLevelAvgConnVectors);
    figure();
    imagesc(matrixOfCorrelationValues);
    xticks([1 2 3 4]);
    xticklabels({'PLV' 'iPLV' 'ampCorr' 'orthAmpCorr'});
    yticks([1 2 3 4]);
    yticklabels({'PLV' 'iPLV' 'ampCorr' 'orthAmpCorr'});
    title(['Group avg. corr. ', freqBandString, ' band, thr']);    
    colormap('jet');
    caxis([-0.8 1]);
    colorbar;
    set(gca, 'FontSize', 14);
    
end % for freqBandIndex = 1 : numberOfFrequencyBands
