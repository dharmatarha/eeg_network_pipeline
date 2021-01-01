clear all;
close all;

pathName = ('D:\psych\rs_fc\alpha\');

meanConnTensor = nan(4, 62, 62);
matrixOfCorrValues = nan(4, 4);

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
    
    meanConnMatrix = meanConnMatrix_group(connectivityTensor);
    meanConnTensor(index, :, :) = meanConnMatrix;
    
end

[numberOfMeasures, roiNo, ~] = size(meanConnTensor);

for firstMeasureIndex = 1 : numberOfMeasures
    for secondMeasureIndex = 1 : numberOfMeasures
        meanConnMatrix_A = squeeze(meanConnTensor(firstMeasureIndex, :, :));
        meanConnMatrix_B = squeeze(meanConnTensor(secondMeasureIndex, :, :));
        linA = meanConnMatrix_A(triu(true(roiNo), 1));
        linB = meanConnMatrix_B(triu(true(roiNo), 1)); 
        matrixOfCorrValues(firstMeasureIndex, secondMeasureIndex) = corr(linA, linB);        
    end
end

imagesc(matrixOfCorrValues);
xticks([1 2 3 4]);
xticklabels({'PLV' 'iPLV' 'ampCorr' 'orthAmpCorr'});
yticks([1 2 3 4]);
yticklabels({'PLV' 'iPLV' 'ampCorr' 'orthAmpCorr'});
colormap('jet'),
colorbar;

