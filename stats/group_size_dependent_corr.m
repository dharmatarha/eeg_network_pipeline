clear all;
close all;

pathName = ('D:\psych\rs_fc\alpha\');
fileName = ('group_alpha_plv.mat');

dataStructure = open([pathName fileName]);

connectivityTensor = dataStructure.connData;

numberOfPermutations = 1000;
arrayOfGroupNumbers = [1 10 20 30 40 50 60 70 80 90 100];
numberOfGroupSizes = numel(arrayOfGroupNumbers);
matrixOfCorrValues = nan(numberOfPermutations, numberOfGroupSizes);

for iterationIndex = 1 : numberOfGroupSizes
    simRes_correlation = connSimTest_group_var_size(connectivityTensor, numberOfPermutations, arrayOfGroupNumbers(iterationIndex), 'corr');
    matrixOfCorrValues(:, iterationIndex) = simRes_correlation';
end

h = violin(matrixOfCorrValues, 'x', [1 10 20 30 40 50 60 70 80 90 100]/10);
xticks([1 10 20 30 40 50 60 70 80 90 100]/10);
xticklabels({'1', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
xlabel('Number of subjects in groups');
ylabel('Between-group correlation');


