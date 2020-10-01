
%% params

gammaOmegaTargets = [1 1];
nodeNo = 62;
layerNo = 40;
stimNo = 4;
topN = 1:100;


%% Load data for edge contributions

% get edge contribution data
fileP = '/media/adamb/bonczData/hyperscan/edgeContributions/alpha_iplv_groupPruning_edgeContrib_edgeRm_corr.mat';
% load
e = load(fileP);
% cohen d values after removing each edge
edgeRmD = e.permRes.cohendAll;

% get "original" within- vs. across-stimuli similarity results
fileP = '/media/adamb/bonczData/hyperscan/alpha_iplv_GroupPruned_cmpFullConnRes_corr.mat';
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
    
    fileP{stim} = ['/media/adamb/bonczData/hyperscan/modularityResults/tensorNormRes_groupPruned/alpha_stim', num2str(stim), '_iplv_modRes_iter_postpr.mat'];

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

% for n = topN
for n = [5:5:50, 60:10:100]
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
    
    
end
   







