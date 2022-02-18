frequencyBands = {'theta', 'alpha', 'beta', 'gamma', 'delta'};
numberOfFrequencyBands = numel(frequencyBands);
connMetrics = {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'};
numberOfConnMetrics = numel(connMetrics);
dirNameBase = '/home/peternagy/NAS502/EEG_resting_state/';

permutationNo = 1000;
simMetric = 'corr';

for freqBandIndex = 1 : numberOfFrequencyBands
    for connMetricIndex = 1 : numberOfConnMetrics
        
        freqBandString = frequencyBands{freqBandIndex};
        connMetricString = connMetrics{connMetricIndex};
        dirName = [dirNameBase freqBandString];
        connectivityFileName = [freqBandString, '_', connMetricString, '_surrConn', '.mat'];
        
        dataStructure = open([dirName, '/', connectivityFileName]);
        connectivityTensor = dataStructure.acrossEpochs.maskedConn;
        [~, edgeContr_correlation] = connSimTest_group_edgeContr_thr(connectivityTensor, 'corr');
        
        realConnMean = squeeze(mean(edgeContr_correlation, 3, 'omitnan'));
            [equalFlag, matchingSetsFlag, roiLabelsPlotting] = roiLabelMatching_mod();
        
        subjectAveragedConnArray = squeeze(mean(connectivityTensor, 1));
        
        % further checks
        if length(roiLabelsPlotting) ~= size(realConnMean, 1)
            error('Length of ROI / node labels cell array does not match the number of nodes in "realConn"!');
        end
        if ~iscolumn(roiLabelsPlotting)
            roiLabelsPlotting = roiLabelsPlotting';
        end
        
        % define trimmingThr and group2color
        trimmingThr = [0, 0.002];
        group2color = [1, 1];
    
        % make sure we have symmetric matrices with zeros at diagonals before reordering rows/columns
        realConnMean = triu(realConnMean, 1) + triu(realConnMean, 1)';
    
        edgeMembership = ones(size(realConnMean));
    
        % define final connectivity matrix with the right name
        connMatrix = realConnMean;
        % define final label cell array with the correct name
        labels = roiLabelsPlotting;
    
        % only upper triangles
        connMatrix(tril(true(size(connMatrix, 1)))) = nan;
        edgeMembership(tril(true(size(edgeMembership, 1)))) = nan;
    
        colorTriplets = [0.25 0.25 0.25];
    
        % call the main plotting function
        mainFig = circleGraphPlot_edges_mod_thr(connMatrix, subjectAveragedConnArray, edgeMembership, colorTriplets,...
        group2color, trimmingThr, labels, 'draw');
        
        switch connMetricString
            case 'plv'
                title(['PLV, ' freqBandString, ', thresholding']);
            case 'iplv'
                title(['iPLV, ', freqBandString, ', thresholding']);
            case 'ampCorr'
                title(['ampCorr, ', freqBandString, ', thresholding']);
            case 'orthAmpCorr'
                title(['orthAmpCorr, ', freqBandString, ', thresholding']);
            otherwise
        end
        
    end
end

