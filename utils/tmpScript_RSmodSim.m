
baseDir = '/media/adamb/bonczData/EEG_resting_state/';
freq = 'alpha';
method = 'iplv';
gammaParam = 1.5;
repNo = 50;

% load data
connF = [baseDir, '/surrConn_', freq, '_', method, '.mat'];
tmp = load(connF);
connData = tmp.acrossEpochs.maskedConn;

% get sizes
[subNo, roiNo, ~] = size(connData);

% preallocate
consPart = nan(subNo, roiNo);
consQ = nan(subNo, 1);
qpc = nan(subNo, 1);
% for similarity
zrandRes = nan(200,200);

% loop through subjects
for s = 1:subNo
    
    % average over epochs
    subData = squeeze(connData(s, :, :));
    
    % normalize
    subData = normalizeMatrix(subData, 'mean', false);
    
    % symm
    subData = triu(subData, 1) +  triu(subData, 1)';
    
    % get B (null model)
    [B, twom] = modularity(subData, gammaParam);
    
    % iterated louvain, repNo repetitions
    part = nan(repNo, roiNo);
    for n = 1:repNo
        part(n, :) = iterated_genlouvain(B, 5000, 0, [], 'moverandw');
    end
    
    % get consensus partition based on nodal assoc matrix
    [S2, Q2, ~, QPC] = consensus_iterative(part);
    % select best one
    tmpIdx = find(Q2==max(Q2), 1);
    if isempty(tmpIdx)
        tmpIdx = 1;
    end
    
    % collect results
    consPart(s, :) = S2(tmpIdx, :);
    consQ(s) = Q2(tmpIdx);
    qpc(s) = QPC;
    
end

% z rand
for s1 = 1:subNo
    for s2 = 1:subNo
        if s1<s2
            zrandRes(s1, s2) = zrand(consPart(s1,:), consPart(s2,:));
        end
    end
end
        

%% save out        
        
saveF = [baseDir, '/fullConnMod_', freq, '_', method, '_', num2str(gammaParam), '.mat'];
save(saveF, 'consPart', 'consQ', 'qpc', 'zrandRes');
            


%% Get consensus of all subjects' data

% get consensus partition based on nodal assoc matrix
[S2, Q2, ~, QPC] = consensus_iterative(consPart);
% select best one
tmpIdx = find(Q2==max(Q2), 1);        
groupP = S2(tmpIdx, :);        

% group connectivity average
groupConn = squeeze(mean(connData, 1));


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
trimmingThr = [0.07];

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
                                        colorTriplets,... 
                                        trimmingThr, ...
                                        roiLabelsPlotting,... 
                                        figTitle);

end




