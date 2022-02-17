
frequencyBands = {'theta', 'alpha', 'beta', 'gamma', 'delta'};
numberOfFrequencyBands = numel(frequencyBands);
connMetrics = {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'};
numberOfConnMetrics = numel(connMetrics);
dirNameBase = '/home/peternagy/NAS502/EEG_resting_state/';

permutationNo = 1000;
simMetric = 'corr';

arrayOfGroupNumbers = [1 10 20 30 40 50 60 70 80 90 100];
numberOfGroupSizes = numel(arrayOfGroupNumbers);

for freqBandIndex = 1 : numberOfFrequencyBands
    
    matrixOfCorrValues = nan(permutationNo, numberOfGroupSizes);
    figure();
    for connMetricIndex = 1 : numberOfConnMetrics
        
        freqBandString = frequencyBands{freqBandIndex};
        connMetricString = connMetrics{connMetricIndex};
        dirName = [dirNameBase freqBandString];
        connectivityFileName = [freqBandString, '_', connMetricString, '_group', '.mat'];
        
        dataStructure = open([dirName, '/', connectivityFileName]);
        connectivityTensor = dataStructure.connData;
        
        for groupSizeIndex = 1 : numberOfGroupSizes
            simRes_correlation = connSimTest_group_var_size(connectivityTensor,...
                arrayOfGroupNumbers(groupSizeIndex), permutationNo, simMetric);
            matrixOfCorrValues(:, groupSizeIndex) = simRes_correlation';
        end
        
        subplot(2, 2, connMetricIndex);
        switch connMetricString
            case 'plv'
                [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10,...
                    'mc', [], 'medc', 'k', 'facecolor', [0, 0.4470, 0.7410], 'facealpha', 1);
                title(['PLV, ' freqBandString, ', no thresholding']);
            case 'iplv'
                [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10,...
                    'mc', [], 'medc', 'k', 'facecolor', [0.9290, 0.6940, 0.1250], 'facealpha', 1);
                title(['iPLV, ', freqBandString, ', no thresholding']);
            case 'ampCorr'
                [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10,...
                    'mc', [], 'medc', 'k', 'facecolor', [0.4660, 0.6740, 0.1880], 'facealpha', 1);
                title(['ampCorr, ', freqBandString, ', no thresholding']);
            case 'orthAmpCorr'
                [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10,...
                    'mc', [], 'medc', 'k', 'facecolor', [0.6350, 0.0780, 0.1840], 'facealpha', 1);
                title(['orthAmpCorr, ', freqBandString, ', no thresholding']);
            otherwise
        end
    
        xticks([1 10 20 30 40 50 60 70 80 90 100]/10);
        xticklabels({'1', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
        xlabel('Number of subjects in groups');
        yticks([-0.2 0 0.2 0.4 0.6 0.8 1]);
        yticklabels({'-0.2', '0', '0.2', '0.4','0.6', '0.8', '1'});
        ylabel('Between-group correlation');
        ylim([-0.2 1])
        hold on;
        connectingLineXvalues = [1 10 20 30 40 50 60 70 80 90 100]/10;
        connectingLineYvalues = MED;
        plot(connectingLineXvalues, connectingLineYvalues, 'k', 'LineWidth', 1);
        set(gca, 'FontSize', 10);
        set(gcf, 'Color', 'w');
        grid on;
    
    end
end

