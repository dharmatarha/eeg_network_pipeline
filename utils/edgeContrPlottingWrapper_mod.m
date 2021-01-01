clear all;
close all;

pathName = ('D:\psych\rs_fc\alpha\');

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
    
    [~, edgeContr_correlation] = connSimTest_group_edgeContr(connectivityTensor, 'corr');
    
    switch index
        case 1
            edgeContr_correlation_PLV = edgeContr_correlation;
            edgeContr_correlation_mean_PLV = squeeze(mean(edgeContr_correlation, 3, 'omitnan'));
        case 2
            edgeContr_correlation_iPLV = edgeContr_correlation;
            edgeContr_correlation_mean_iPLV = squeeze(mean(edgeContr_correlation, 3, 'omitnan'));
        case 3
            edgeContr_correlation_ampCorr = edgeContr_correlation;
            edgeContr_correlation_mean_ampCorr = squeeze(mean(edgeContr_correlation, 3, 'omitnan'));
        case 4
            edgeContr_correlation_orthAmpCorr = edgeContr_correlation;
            edgeContr_correlation_mean_orthAmpCorr = squeeze(mean(edgeContr_correlation, 3, 'omitnan'));
        otherwise
    end
    
    
    realConnMean = squeeze(mean(edgeContr_correlation, 3, 'omitnan'));
    
    
    [equalFlag, matchingSetsFlag, roiLabelsPlotting] = roiLabelMatching_mod();
    % further checks
    if length(roiLabelsPlotting) ~= size(realConnMean, 1)
        error('Length of ROI / node labels cell array does not match the number of nodes in "realConn"!');
    end
    if ~iscolumn(roiLabelsPlotting)
        roiLabelsPlotting = roiLabelsPlotting';
    end
    
    % define trimmingThr and group2color
    trimmingThr = [0, 0.001];
    group2color = [1, 1];
    
    % make sure we have symmetric matrices with zeros at diagonals before reordering rows/columns
    realConnMean = triu(realConnMean, 1) + triu(realConnMean, 1)';
    
    edgeMembership = ones(size(realConnMean));
    
    % define final connectivity matrix with the right name
    connMatrix = realConnMean;
    % define final label cell array with the correct name
    labels = roiLabelsPlotting;
    
    % only upper triangles
    connMatrix(tril(true(size(connMatrix, 1)))) = nan;
    edgeMembership(tril(true(size(edgeMembership, 1)))) = nan;
    
    colorTriplets = [0.25 0.25 0.25];
    
    % call the main plotting function
    mainFig = circleGraphPlot_edges_mod(connMatrix, edgeMembership, colorTriplets,...
        group2color, trimmingThr, labels, 'draw');
    
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






