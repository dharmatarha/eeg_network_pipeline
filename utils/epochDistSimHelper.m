
% x_values = squeeze(distSimResult(:, 2, :));
% x_values = x_values(:);
% y_values = squeeze(distSimResult(:, 1, :));
% y_values = y_values(:);
% figure();
% scatter(x_values, y_values);

% load('/home/peternagy/NAS502/EEG_resting_state/delta/lh_s04_delta_ampCorr.mat');
% epochMatrixA = squeeze(connRes(37, :, :));
% epochMatrixB = squeeze(connRes(39, :, :));
% load('/home/peternagy/NAS502/EEG_resting_state/delta/lvid_s08_delta_ampCorr.mat');
% epochMatrixA = squeeze(connRes(137, :, :));
% epochMatrixB = squeeze(connRes(139, :, :));
% linA = epochMatrixA(triu(true(roiNo), 1));
% linB = epochMatrixB(triu(true(roiNo), 1));
% neighboringEpochCorrelation = corr(linA, linB)
% meanEpochMatrixA = squeeze(mean(connRes(1:2:75, :, :), 1));
% meanEpochMatrixB = squeeze(mean(connRes(77:2:151, :, :), 1));
% linA = meanEpochMatrixA(triu(true(roiNo), 1));
% linB = meanEpochMatrixB(triu(true(roiNo), 1));
% correlationOfAveragedEpochs = corr(linA, linB)

frequencyBands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
numberOfFrequencyBands = numel(frequencyBands);
connMetrics = {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'};
numberOfConnMetrics = numel(connMetrics);
dirNameBase = '/home/peternagy/NAS502/EEG_resting_state/';

for freqBandIndex = 1 : numberOfFrequencyBands
    figure();
    hold on;
    freqBandString = frequencyBands{freqBandIndex};
    for connMetricIndex = 1 : numberOfConnMetrics
        
        connMetricString = connMetrics{connMetricIndex};
        dirName = [dirNameBase freqBandString];
        fileName = [freqBandString, '_', connMetricString, '_', 'epochDistSim.mat'];
        load([dirName '/', fileName]);
        
        groupLevelMeanValues = nan(74, 1);
        groupLevelSEValues = nan(74, 1);

        for epochDistance = 1 : 74
            corrValuesTmpArray = [];
            for subjectIndex = 1 : size(distSimResult, 2)
                for rowIndex = 1 : size(distSimResult, 3)
                    for columnIndex = 1 : size(distSimResult, 4)
                        if distSimResult(subjectIndex, 1, rowIndex, columnIndex) == epochDistance
                            corrValuesTmpArray = [corrValuesTmpArray distSimResult(subjectIndex, 2, rowIndex, columnIndex)];
                        end
                    end
                end
            end
            groupLevelMeanValues(epochDistance) = mean(corrValuesTmpArray);
            groupLevelSEValues(epochDistance) = std(corrValuesTmpArray)/sqrt(length(corrValuesTmpArray));
        end

        SECurve1 = groupLevelMeanValues + groupLevelSEValues;
        SECurve2 = groupLevelMeanValues - groupLevelSEValues;
        epochDistances = 1 : numel(groupLevelMeanValues);
        epochDistances2 = [epochDistances fliplr(epochDistances)];
        SEInBetween = [SECurve1' fliplr(SECurve2')];
        h = fill(epochDistances2, SEInBetween, [128 128 128]/255);
        set(h, 'edgealpha', 0);
        switch connMetricString
            case 'plv'
                plot(epochDistances, groupLevelMeanValues', 'Color', [0, 0.4470, 0.7410]);
            case 'iplv'
                plot(epochDistances, groupLevelMeanValues', 'Color', [0.9290, 0.6940, 0.1250]);
            case 'ampCorr'
                plot(epochDistances, groupLevelMeanValues', 'Color', [0.4660, 0.6740, 0.1880]);
            case 'orthAmpCorr'
                plot(epochDistances, groupLevelMeanValues', 'Color', [0.6350, 0.0780, 0.1840]);
        end
    end % connMetricIndex = 1 : numberOfConnMetrics
    switch freqBandString
        case 'delta'
            title('Epoch similarity, delta band');
        case 'theta'
            title('Epoch similarity, theta band');
        case 'alpha'
            title('Epoch similarity, alpha band');
        case 'beta'
            title('Epoch similarity, beta band');
        case 'gamma'
            title('Epoch similarity, gamma band');
    end
    xlabel('Epoch distance (correlation)');
    ylabel('Epoch similiarty');
    f=get(gca,'Children');
    f = flipud(f);
    legend([f(2), f(4), f(6), f(8)], 'PLV', 'iPLV', 'ampCorr', 'orthAmpCorr');
end % freqBandIndex = 1 : numberOfFrequencyBands

