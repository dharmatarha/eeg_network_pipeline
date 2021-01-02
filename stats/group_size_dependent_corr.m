clear all;
close all;

pathName = ('D:\psych\rs_fc\alpha\');

numberOfPermutations = 1000;
arrayOfGroupNumbers = [1 10 20 30 40 50 60 70 80 90 100];
numberOfGroupSizes = numel(arrayOfGroupNumbers);
matrixOfCorrValues = nan(numberOfPermutations, numberOfGroupSizes);

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
    
    for iterationIndex = 1 : numberOfGroupSizes
        simRes_correlation = connSimTest_group_var_size(connectivityTensor, numberOfPermutations, arrayOfGroupNumbers(iterationIndex), 'corr');
        matrixOfCorrValues(:, iterationIndex) = simRes_correlation';
    end
    
    subplot(2, 2, index);
    switch index
        case 1
            [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10, 'mc', [], 'medc', 'k', 'facecolor', [0, 0.4470, 0.7410], 'facealpha', 1);
        case 2
            [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10, 'mc', [], 'medc', 'k', 'facecolor', [0.9290, 0.6940, 0.1250], 'facealpha', 1);
        case 3
            [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10, 'mc', [], 'medc', 'k', 'facecolor', [0.4660, 0.6740, 0.1880], 'facealpha', 1);
        case 4
            [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10, 'mc', [], 'medc', 'k', 'facecolor', [0.6350, 0.0780, 0.1840], 'facealpha', 1);
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
    connectingLineYvalues = MX;
    plot(connectingLineXvalues, connectingLineYvalues, 'k', 'LineWidth', 1);
    set(gca, 'FontSize', 10);
    grid on;
    
    switch index
        case 1
            title('PLV');
        case 2
            title('iPLV');
        case 3
            title('ampCorr');
        case 4
            title('orthAmpCorr');
        otherwise
    end

end


% for iterationIndex = 1 : numberOfGroupSizes
%     simRes_correlation = connSimTest_group_var_size(connectivityTensor, numberOfPermutations, arrayOfGroupNumbers(iterationIndex), 'corr');
%     matrixOfCorrValues(:, iterationIndex) = simRes_correlation';
% end

% [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10);
% xticks([1 10 20 30 40 50 60 70 80 90 100]/10);
% xticklabels({'1', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
% xlabel('Number of subjects in groups');
% ylabel('Between-group correlation');
% hold on;
% connectingLineXvalues = [1 10 20 30 40 50 60 70 80 90 100]/10;
% connectingLineYvalues = MX;
% plot(connectingLineXvalues, connectingLineYvalues, 'k');

