%% Base vars, loading data

baseDir = '/media/adamb/bonczData/EEG_resting_state/alpha';
freq = 'alpha';
method = 'orthAmpCorr';
gammaParam = 1.05;

% load connectivity data
connF = [baseDir, '/group_', freq, '_', method, '.mat'];
tmp = load(connF);
connData = tmp.connData;

% load modularity data
connModF = [baseDir, '/fullConnMod_', freq, '_', method, '_', num2str(gammaParam), '.mat'];
tmp = load(connModF);
consPart = tmp.consPart;


%% Get consensus of all subjects' data

% get consensus partition based on nodal assoc matrix
[S2, Q2, ~, QPC] = consensus_iterative(consPart);
% select best one
tmpIdx = find(Q2==max(Q2), 1);        
groupP = S2(tmpIdx, :);        

% group connectivity average
groupConn = squeeze(mean(mean(connData, 1), 2));


%% Plot group-level connectivity with group-level modules

% get labels
roiF = '/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat';
load(roiF, 'roisShort');
roiLabels = roisShort;

% get colors
colorF = '/home/adamb/eeg_network_pipeline/utils/colorTriplets.mat';
load(colorF, 'colorTriplets24');
colorTriplets = colorTriplets24;
myColors = [0, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.3660, 0.6740, 0.1];

% trimming threshold
trimmingThr = [0.24];

% figure title
figTitle = 'Group-level (consensus partition) modularity';

% see if the supplied labels match any of the standard sets
[equalFlag, matchingSetsFlag, roiLabelsPlotting] = roiLabelMatching(roiLabels);

if ~equalFlag && matchingSetsFlag
    % reorder connectivity data to match new label order
    [connMatrix, old2new] = matrixReorder(groupConn, roiLabels', roiLabelsPlotting);
    % apply the same re-ordering to ROI/node module indices
    modIndicesVector = groupP(old2new);

    % plotting
    [mainFig, subFig] = circleGraphPlot(connMatrix,... 
                                        modIndicesVector,... 
                                        myColors,... 
                                        trimmingThr, ...
                                        roiLabelsPlotting,... 
                                        figTitle);

end