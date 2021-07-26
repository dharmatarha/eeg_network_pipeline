%% Script for Figure 2A
% The script generates the group-averaged connectivity matrix, 
% reordered for easy-to-read labeling.
%
% 

%% Base params

method = 'plv'; 
freq = 'alpha';

% type of thresholding (if any), one of {'unthr', 'thrSub', 'thrGroup'}
thr = 'thrSub';

baseDir = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';


%% Load data, determined by "method" and "freq"

fileP = fullfile(baseDir, freq, ['group_surrResults_', freq, '_', method, '.mat']);
res = load(fileP);


%% Get mean connectivity, depending on "thr"

switch thr
    case 'unthr'
        connData = res.meanConn;
    case 'thrSub'
        % get the data from the averaged, subject-level thresholded matrices, both
        % negative and positive diffs
        connData = res.meanMaskedConnPos + res.meanMaskedConnNeg;   
    case 'thrGroup'
        connData = res.groupMaskedConnPos + res.groupMaskedConnNeg;
end  % switch
        
% get averaged conn matrix
meanConn = squeeze(mean(mean(connData, 1), 2));


%% Load ROI label names

labels = load('/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat');
labels = labels.roisShort;  % short names


%% Get reordered heatmap

[connMatrixRo, labelsRo, h, hs] = heatmapReordered(meanConn, labels);


%% Save out, quit

saveas(gcf, fullfile(baseDir, ['meanConn_', thr, '_', freq, '_', method, '.svg']));
saveas(gcf, fullfile(baseDir, ['meanConn_', thr, '_', freq, '_', method, '.png']));


