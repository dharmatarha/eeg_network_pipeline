
frequencyBands = {'theta', 'alpha', 'beta', 'gamma', 'delta'};
connMetrics = {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'};
dirNameBase = '/home/peternagy/NAS502/EEG_resting_state/';
subjectNo = 200;
permutationNo = 1000;

freqBandString = frequencyBands{1};
connMetricString = connMetrics{2};
dirName = [dirNameBase freqBandString];
connectivityFileName = [freqBandString, '_', connMetricString, '_group', '.mat'];
if ismember(connMetricString, {'plv', 'iplv'})
	truncated = 'truncated';
else
	truncated = 'nontruncated';
end
simResCorr = connSimTest_subject_thresholding([dirName, '/', connectivityFileName], ...
dirName, freqBandString, connMetricString, truncated);
        
simResCorrFileName = [freqBandString, '_', connMetricString, '_simResCorr_thr', '.mat'];
save([dirName, '/', simResCorrFileName], 'simResCorr');

% simResCorr_vector = reshape(simResCorr, [], 1);
% tensorOfCorrValues(freqBandIndex, connMetricIndex, :) = simResCorr_vector;
% for freqBandIndex = 1 : numberOfFrequencyBands
%     matrixOfCorrValues = squeeze(tensorOfCorrValues(freqBandIndex, :, :))';
%     figure();
%     [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 2 3 4], 'mc', [], ...
%         'medc', 'k', 'facecolor', [0, 0.4470, 0.7410; 0.9290, 0.6940, ...
%         0.1250; 0.4660, 0.6740, 0.1880; 0.6350, 0.0780, 0.1840], 'facealpha', 1);
%     xticks([1 2 3 4]);
%     xticklabels({'PLV' 'iPLV' 'ampCorr' 'orthAmpCorr'});
%     xlim([0 5]);
%     ylabel('Within-subject correlation');
%     switch freqBandIndex
%         case 1
%             title('Theta band');
%         case 2
%             title('Alpha band');
%         case 3
%             title('Beta band');
%         case 4
%             title('Gamma band');
%         case 5
%             title('Delta band');
%     end
%     set(gca, 'FontSize', 20);
%     set(gcf, 'Color', 'w');
%     grid on;
% end

