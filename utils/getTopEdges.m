
%% params

gammaOmegaTargets = [1 1];
nodeNo = 62;
layerNo = 40;
stimNo = 4;
topN = 1:100;


%% Load data for edge contributions

% get edge contribution data
fileP = '/media/stim/bonczData/hyperscan/edgeContributions/alpha_iplv_groupPruning_edgeContrib_edgeRm_corr.mat';
% load
e = load(fileP);
% cohen d values after removing each edge
edgeRmD = e.permRes.cohendAll;

% get "original" within- vs. across-stimuli similarity results
fileP = '/media/stim/bonczData/hyperscan/alpha_iplv_GroupPruned_cmpFullConnRes_corr.mat';
% load
o = load(fileP);
% base cohen D value
baseD = o.permRes(5).cohenD;  % fifth test was performed on data across all stimuli pairings

% variable of interest: changes relative to the base cohen D value
edgeContr = baseD-edgeRmD;

% matrix version
edgeContrM = nan(nodeNo);
edgeContrM(triu(true(nodeNo), 1)) = edgeContr;


%% Load modularity data

% preallocate
fileP = cell(stimNo, 1);
% var holding all module partitions
allModules = nan(nodeNo, layerNo*stimNo);

for stim = 1:stimNo
    
    fileP{stim} = ['/media/stim/bonczData/hyperscan/modularityResults/tensorNormRes_groupPruned/alpha_stim', num2str(stim), '_iplv_modRes_iter_postpr.mat'];

    modRes = load(fileP{stim});

    % get modularity around target gamma, omega pairing
    [~, gIdx] = min(abs(modRes.gammaValues-gammaOmegaTargets(1)));
    [~, oIdx] = min(abs(modRes.omegaValues-gammaOmegaTargets(2)));
    % provide feedback about selected gamma & omega values
    disp([char(10), 'Requested [gamma, omega] values were: ', num2str(gammaOmegaTargets),...
        char(10), 'Closest match in the data was: ',... 
        num2str([modRes.gammaValues(gIdx), modRes.omegaValues(oIdx)])]);

    modules = modRes.res.consSim(gIdx, oIdx, :);
    modules = double(squeeze(modules));  % res.consSim is uint16 by default if coming from multiCommDetectWrapper
    % reshape modularity memberships into one module vector per layer/epoch
    allModules(:, (stim-1)*layerNo+1:stim*layerNo) = reshape(modules, [nodeNo, layerNo]);


end


%% Get consensus similarities

% [cons, cs, ps] = consensus_similarity(allModules');
% [s2, q2, X, qpc] = consensus_iterative(allModules');
% i = find(q2==max(q2));
% cons2 = s2(i(1), :);

[cons, cs, ps] = consensus_similarity(allModules(:, 121:160)');
[s2, q2, X, qpc] = consensus_iterative(allModules(:, 121:160)');
i = find(q2==max(q2));
cons2 = s2(i(1), :);


%% Select top N edges & test within-module connectivities 

edgeSort = sort(edgeContr, 'descend');

aW = nan(28,2,16);
bW = nan(16,2);
% for n = topN
counter = 1;
for n = 10:5:80
    % get critical value of contribution for top n edges
    critValue = edgeSort(n+1);
    % get binary matrix containing ones for top n edges 
    critM = double(edgeContrM>critValue);
    
    % get within- and between-module connectivities on critical matrix
    [withinConn, betweenConn, withinEdgeIndices, betweenEdgeIndices] = modConnections(cons2, critM);
    
    % perm test
    pBetweenC = nan(1000,1);
    pWithinC = nan(1000, length(unique(cons2)));
    for z = 1:10000
        % permute modularity partition
        pCons = cons2(randperm(length(cons2)));
        % get connectivities
        [pWithinC(z,:), pBetweenC(z), ~, ~] = modConnections(pCons, critM);
    end
    
    
    pEstimates = estimatedP(withinConn, pWithinC);
    pEstimates(isnan(withinConn)) = nan;
    pEstBetween = estimatedP(betweenConn, pBetweenC);
    
    disp(['n = ', num2str(n)]);
    disp(pEstimates);
    disp(pEstBetween);
    
    aW(:, :, counter) = pEstimates;
    bW(counter, 1:2) = pEstBetween;
    counter = counter + 1;
    
end
   

%% Load connectivity data for plotting

fileP = ['/media/stim/bonczData/hyperscan/edgePruningResults_iplv/alpha/avg_alpha_edgePruningInfo.mat'];
r = load(fileP);
realConn = r.meanConn.realConnAvg;
connMatrix = realConn;
clearvars r

myModule = cons2;


labels = load('/home/stim/eeg_network_pipeline/utils/roiNamesInOrder.mat');
labels = labels.roisShort;
c = load('/home/stim/eeg_network_pipeline/utils/colorTriplets.mat');
colorTriplets = c.colorTriplets24;

% new colors
myColors = [colorTriplets(5,:); repmat([0.5,0.5,0.5], [62, 1])];

% mod2color
mod2color = [10, 1; [1:9, 11:49]', [2:49]'];

% get fdr mask
gp = ['/media/stim/bonczData/hyperscan/surrEdgeEstimates/truncatedNormal/alpha_groupEdgePruningInfo.mat'];
gp = load(gp);
fdrMask = gp.fdrMask;
connMatrix(~fdrMask) = 0;

% average across layers
connMatrix = mean(mean(connMatrix, 4), 3);

% relabel
% see if the supplied labels match any of the standard sets
[equalFlag, matchingSetsFlag, roiLabelsPlotting] = roiLabelMatching(labels);

% if the labels were not an exact match but overlapped with standard
% sets, rearrange labels + data
if ~equalFlag && matchingSetsFlag
    % rearrange connectivity matrix based on new ROI label order
    [connMatrixNew, old2new] = matrixReorder(connMatrix, labels', roiLabelsPlotting);
    
    % apply the same re-ordering to ROI/node module indices
    myModuleNew = myModule(old2new);
    myModuleNew2 = myModuleNew;
    counter = 1;
    for z = 1:62
        if myModuleNew2(z)~=10
            myModuleNew2(z) = counter;
            counter = counter+1;
        end
    end

    % plotting
    [mainFig, subFig] = circleGraphPlot(connMatrixNew, myModuleNew2,... 
                            myColors, mod2color,... 
                            [0, 0.08], roiLabelsPlotting, 'Module 10, alpha');

end




%% Plot top n edges




% rearrange edge grouping matrix based on new ROI label order
[critMnew, old2new2] = matrixReorder(critM, labels', roiLabelsPlotting);

critMnew(tril(true(size(critMnew, 1)))) = nan;

group2color = [1, 1];

% call the main plotting function
mainFig = circleGraphPlot_edges(connMatrixNew, critMnew, myColors,...
                                group2color, [0, 0.08], roiLabelsPlotting, 'draw');











