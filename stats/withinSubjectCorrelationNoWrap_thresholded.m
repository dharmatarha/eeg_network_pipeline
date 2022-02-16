function withinSubjectCorrelationNoWrap_thresholded(freq, connMethod)
%% Wrapper for calling the function connSimTest_subject_thresholding
%
% For a specific frequency band and connectivity method. Saves out
% the output to a file named FREQ_CONNMETHOD_simResCorr_thr.mat. 
% 
% Contains hardcoded path variables.
%


%% Input checks

if nargin ~= 2
    error('Requires input args "freq" and "connMethod"!');
end
if ~ischar(freq) || ~ismember(freq, {'theta', 'alpha', 'beta', 'gamma', 'delta'})
    error('Input arg "freq" should be a char array, one of {"delta", "theta", "alpha", "beta", "gamma"}!');
end
if ~ischar(freq) || ~ismember(connMethod, {'plv', 'iplv', 'ampCorr', 'orthAmpCorr'})
    error('Input arg "connMethod" should be a char array, one of {"plv", "iplv", "ampCorr", "orthAmpCorr"}!');
end


%% Basic params, hardcoded values

dirNameBase = '/home/peternagy/NAS502/EEG_resting_state/';
dirName = [dirNameBase freq];

% group connectivity data file
connectivityFileName = [freq, '_', connMethod, '_group', '.mat'];
% results file
simResCorrFileName = [freq, '_', connMethod, '_simResCorr_thr', '.mat'];

permutationNo = 1000;
simMetric = 'corr';

if ismember(connMethod, {'plv', 'iplv'})
	truncated = 'truncated';
else
	truncated = 'nontruncated';
end


%% Call the within-subject similarity calculator function

simResCorr = connSimTest_subject_thresholding([dirName, '/', connectivityFileName], ...   % connectivity data file
                                            dirName,...                                   % surrogate data folder  
                                            freq,...                                      % frequency band
                                            connMethod,...                                % connectivity method
                                            truncated,...                                 % truncated or non-truncated normals for stats
                                            permutationNo,...                             % number of permutations for testing within-subject connectivity similarity
                                            simMetric);                                   % similarity metric


%% Save out resutls

save([dirName, '/', simResCorrFileName], 'simResCorr');


return




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

