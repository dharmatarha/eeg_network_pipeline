clear all;
close all;

matrixOfCorrValues = nan(200000, 4);

load('D:\psych\rs_fc\alpha\connSimRes_alpha_plv.mat');
simRes_correlation_vector = reshape(withinSub_corr, [], 1);
matrixOfCorrValues(:, 1) = simRes_correlation_vector;

load('D:\psych\rs_fc\alpha\connSimRes_alpha_iplv.mat');
simRes_correlation_vector = reshape(withinSub_corr, [], 1);
matrixOfCorrValues(:, 2) = simRes_correlation_vector;

load('D:\psych\rs_fc\alpha\connSimRes_alpha_ampCorr.mat');
simRes_correlation_vector = reshape(withinSub_corr, [], 1);
matrixOfCorrValues(:, 3) = simRes_correlation_vector;

load('D:\psych\rs_fc\alpha\connSimRes_alpha_orthAmpCorr.mat');
simRes_correlation_vector = reshape(withinSub_corr, [], 1);
matrixOfCorrValues(:, 4) = simRes_correlation_vector;

[h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 2 3 4], 'mc', [], 'medc', 'k', ...
    'facecolor', [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.4660, 0.6740, 0.1880; 0.6350, 0.0780, 0.1840]);
xticks([1 2 3 4]);
xticklabels({'PLV' 'iPLV' 'ampCorr' 'orthAmpCorr'});
xlim([0 5]);
ylabel('Within-subject correlation');
title('Alpha band');
set(gca, 'FontSize', 16);
grid on;