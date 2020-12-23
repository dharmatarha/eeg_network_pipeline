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
    [h,L,MX,MED] = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10);
    xticks([1 10 20 30 40 50 60 70 80 90 100]/10);
    xticklabels({'1', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
    xlabel('Number of subjects in groups');
    ylabel('Between-group correlation');
    ylim([-0.2 1])
    hold on;
    connectingLineXvalues = [1 10 20 30 40 50 60 70 80 90 100]/10;
    connectingLineYvalues = MX;
    plot(connectingLineXvalues, connectingLineYvalues, 'k');
    
    switch index
        case 1
            title('Group Alpha PLV');
        case 2
            title('Group Alpha iPLV');
        case 3
            title('Group Alpha ampCorr');
        case 4
            title('Group Alpha orthAmpCorr');
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

