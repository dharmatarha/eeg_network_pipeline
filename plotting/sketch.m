% scratch 

% basic data
meanConnFile = '/home/adamb/eeg_network_pipeline/dev/edgePruningResults/alpha/avg_alpha_edgePruningInfo.mat';
membershipFile = '/home/adamb/eeg_network_pipeline/dev/edgePruningResults/alpha/avg_alpha_edgePruningInfo_modules.mat';
colorFile = '/home/adamb/eeg_network_pipeline/utils/colorTriplets.mat';
load(meanConnFile, 'prunedConnAvg');
load(membershipFile, 'memberships');
load(colorFile, 'colorTriplets');

% select graph (and corresponding memberships)
epochIdx = 1;
condIdx = 1;
g = prunedConnAvg(:, :, epochIdx, condIdx);
m = memberships(:, epochIdx, condIdx);

% load labels / ROI names
rois = '/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat';
l = load(rois); labels = l.roisShort;


newLabels = {'lateralorbitofrontal L', 'medialorbitofrontal L', 'parsorbitalis L', 'parstriangularis L', 'parsopercularis L', 'rostralmiddlefrontal L', 'caudalmiddlefrontal L', 'superiorfrontal L', 'precentral L',...
    'lateralorbitofrontal R', 'medialorbitofrontal R', 'parsorbitalis R', 'parstriangularis R', 'parsopercularis R', 'rostralmiddlefrontal R', 'caudalmiddlefrontal R', 'superiorfrontal R', 'precentral R',...
    'rostralanteriorcingulate L', 'caudalanteriorcingulate L', 'posteriorcingulate L', 'isthmuscingulate L',... 
    'rostralanteriorcingulate R', 'caudalanteriorcingulate R', 'posteriorcingulate R', 'isthmuscingulate R',...
    'transversetemporal L', 'superiortemporal L', 'middletemporal L', 'inferiortemporal L', 'entorhinal L', 'parahippocampal L', 'fusiform L', 'insula L',...
    'transversetemporal R', 'superiortemporal R', 'middletemporal R', 'inferiortemporal R', 'entorhinal R', 'parahippocampal R', 'fusiform R', 'insula R',...
    'supramarginal L', 'inferiorparietal L', 'superiorparietal L', 'postcentral L', 'paracentral L', 'precuneus L',...
    'supramarginal R', 'inferiorparietal R', 'superiorparietal R', 'postcentral R', 'paracentral R', 'precuneus R',...
    'lateraloccipital L', 'pericalcarine L', 'lingual L', 'cuneus L',...
    'lateraloccipital R', 'pericalcarine R', 'lingual R', 'cuneus R'};

newLabels = {'lateralorbitofrontal L', 'medialorbitofrontal L', 'parsorbitalis L', 'parstriangularis L', 'parsopercularis L', 'rostralmiddlefrontal L', 'caudalmiddlefrontal L', 'superiorfrontal L', 'precentral L',...
    'rostralanteriorcingulate L', 'caudalanteriorcingulate L', 'posteriorcingulate L', 'isthmuscingulate L',...
    'transversetemporal L', 'superiortemporal L', 'middletemporal L', 'inferiortemporal L', 'entorhinal L', 'parahippocampal L', 'fusiform L', 'insula L',...
    'supramarginal L', 'inferiorparietal L', 'superiorparietal L', 'postcentral L', 'paracentral L', 'precuneus L',...
    'lateraloccipital L', 'pericalcarine L', 'lingual L', 'cuneus L',...
    'cuneus R', 'lingual R', 'pericalcarine R', 'lateraloccipital R',...
    'precuneus R', 'paracentral R', 'postcentral R', 'superiorparietal R', 'inferiorparietal R', 'supramarginal R',...
    'insula R', 'fusiform R', 'parahippocampal R', 'entorhinal R', 'inferiortemporal R', 'middletemporal R', 'superiortemporal R', 'transversetemporal R',...
    'isthmuscingulate R', 'posteriorcingulate R', 'caudalanteriorcingulate R', 'rostralanteriorcingulate R',...
    'precentral R', 'superiorfrontal R', 'caudalmiddlefrontal R', 'rostralmiddlefrontal R', 'parsopercularis R', 'parstriangularis R', 'parsorbitalis R', 'medialorbitofrontal R', 'lateralorbitofrontal R',...
    };


shiftLabels = 15;
newLabels = [newLabels(end-shiftLabels:end), newLabels(1:end-shiftLabels-1)];

% create symmetric adjacency matrix with zeros at diagonal
g = triu(g, 1) + triu(g, 1)';
g(isnan(g)) = 0;

[gNew, old2new] = matrixReorder(g, l.rois, newLabels);

g = gNew;
labels = newLabels;
m = m(old2new);

connMatrix = g; membership = m;

% % just for plotting we set low values to zero
% g(g<0.2) = 0;  



% get colors to modules
moduleIndices = unique(m);
colorMap = zeros(length(m), 3);
for i = 1:length(moduleIndices)
    colorMap(m==moduleIndices(i), :) = repmat(colorTriplets(i, :), [sum(m==moduleIndices(i)), 1]);
end


% % draw circular network figure
% circularGraph(g,'Colormap', colorMap,'Label',labels);

% create built-in graph object
G = graph(g, labels, 'upper');
% edge weights into vector
weights = G.Edges.Weight;
% edge ending nodes into cell array
nodesPerEdge = G.Edges.EndNodes;

% go through all modules, set different edge properties per module
edgeColors = repmat([0.5, 0.5, 0.5], [size(weights, 1), 1]);  % preallocate variable for edge color, filled with "gray"
edgeWidth = weights;
moduleEdges = zeros(size(weights, 1), length(moduleIndices));
for i = 1:length(moduleIndices)
    moduleNodes = labels(m == moduleIndices(i));
    moduleEdges(:, i) = ismember(nodesPerEdge(:, 1), moduleNodes) & ismember(nodesPerEdge(:, 2), moduleNodes);
    edgeColors(logical(moduleEdges(:, i)), :) = repmat(colorTriplets(i, :), [sum(moduleEdges(:, i)), 1]);
    edgeWidth(logical(moduleEdges(:, i))) = weights(logical(moduleEdges(:, i)))*10;
end

outModEdgeIdx = ismember(edgeColors, [0.5, 0.5, 0.5], 'rows');
edgeStyle = repmat({'-'}, [size(weights, 1), 1]);
edgeStyle(outModEdgeIdx) = repmat({'none'}, [sum(outModEdgeIdx, 1), 1]);

% create base graph plot
H = G.plot('Layout', 'circle',... 
    'LineWidth', edgeWidth,... 
    'EdgeColor', edgeColors,... 
    'EdgeAlpha', 0.3,... 
    'NodeColor', colorMap,...
    'MarkerSize', 10,...
    'LineStyle', edgeStyle);

% highlight a module
% nodes
modNodes = labels(m==0);
highlight(H, modNodes, 'MarkerSize', 10, 'NodeColor', colorTriplets(2,:));
% edges
% get edge indices for connections inside the module
modG = subgraph(G, modNodes);
subgNodes = modG.Edges.EndNodes;
subgWeights = modG.Edges.Weight;
for subedge = 1:size(subgNodes, 1)
    highlight(H, subgNodes(subedge, 1), subgNodes(subedge, 2), 'LineWidth', subgWeights(subedge)*12, 'EdgeColor', colorTriplets(2,:));
end



