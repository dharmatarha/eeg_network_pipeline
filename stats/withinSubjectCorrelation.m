clear all;
close all;

pathName = ('D:\psych\rs_fc\alpha\');

matrixOfCorrValues = nan(200000, 4);
matrixOfCorrValues_patAvg = nan(200, 4);

for index = 1 : 4
    switch index
        case 1
            fileName = ('group_alpha_plv.mat');
            dataStructure = open([pathName fileName]);
            connectivityTensor = dataStructure.connData;
        case 2
            fileName = ('group_alpha_iplv.mat');
            dataStructure = open([pathName fileName]);
            connectivityTensor = dataStructure.connData;
        case 3
            fileName = ('group_alpha_ampCorr.mat');
            dataStructure = open([pathName fileName]);
            connectivityTensor = dataStructure.connData;
        case 4
            fileName = ('group_alpha_orthAmpCorr.mat');
            dataStructure = open([pathName fileName]);
            connectivityTensor = dataStructure.connData;
        otherwise
    end
    
    simRes_correlation = connSimTest_subject(connectivityTensor, 'corr');
    matrixOfCorrValues_patAvg(:, index) = mean(simRes_correlation, 2);
    simRes_correlation_vector = reshape(simRes_correlation, [], 1);
    matrixOfCorrValues(:, index) = simRes_correlation_vector;
end

[h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 2 3 4]);
xticks([1 2 3 4]);
xticklabels({'PLV' 'iPLV' 'ampCorr' 'orthAmpCorr'});
xlim([0 5]);

% figure();
% [h2,L2,MX2,MED2] = violin(matrixOfCorrValues_patAvg, 'x', [1 2]);
% xticks([1 2]);
% xticklabels({'iplv' 'plv'});
% xlim([0 3]);

% hold on;
% connectingLineXvalues = [1 2];
% connectingLineYvalues = MX1;
% plot(connectingLineXvalues, connectingLineYvalues, 'k');

