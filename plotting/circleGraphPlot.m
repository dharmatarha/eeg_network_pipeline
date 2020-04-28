function graphPlot = circleGraphPlot(connMatrix, membership, colorTriplets, varargin)
%% Plotting network connectivity with module-structure in a circle layout
%
% USAGE: graphPlot = circleGraphPlot(connMatrix, 
%                                   membership,
%                                   colorTriplets,
%                                   trimmingThr=0.2, 
%                                   labels=[], 
%                                   silent=0)
% 
% Create a circle plot of the graph contructed from the supplied
% connectivity (adjacency) matrix and module membership vector. 
% Its details are fine-tuned for undirected EEG connectivity data. 
% Returns the plot object and optionally displays it in a figure window  
% 
% FOR UNDIRECTED GRAPH!
%
% Mandatory inputs:
% connMatrix      - Numeric matrix containing connectivity (adjacency) 
%               values (square matrix). Only upper triangle is used for 
%               graph construction. 
% membership      - Numeric vector containing the module membership of each
%               node in the graph. Values must be in range 0:1:20 (zero is
%               a valid module identifier).
% colorTriplets   - Numeric matrix with three columns, each row defines an
%               RGB triplet for module (and within-module edge) coloring. 
% 
% Optional inputs:
% trimmingThr     - One- or two-element vector containing threshold(s) 
%               for trimming (deleting) weak connections before plotting.
%               If trimmingThr is only one value, the same threshold is
%               applied to all connections. If two values, the first one is
%               applied to within-module, the second to between-module
%               connections. Defaults to [0.2].
% labels          - Cell array of node labels / names.  

% basic data
meanConnFile = '/home/adamb/eeg_network_pipeline/dev/edgePruningResults/alpha/avg_alpha_edgePruningInfo.mat';
membershipFile = '/home/adamb/eeg_network_pipeline/dev/edgePruningResults/alpha/avg_alpha_edgePruningInfo_modules.mat';
load(meanConnFile, 'prunedConnAvg');
load(membershipFile, 'memberships');

% select graph (and corresponding memberships)
epochIdx = 1;
condIdx = 1;
g = prunedConnAvg(:, :, epochIdx, condIdx);
m = memberships(:, epochIdx, condIdx);

% create symmetric adjacency matrix with zeros at diagonal
g = triu(g, 1) + triu(g, 1)';
g(isnan(g)) = 0;
% % just for plotting we set low values to zero
% g(g<0.2) = 0;  

% load RGB triplets
load('colorTriplets.mat');

% get colors to modules
moduleIndices = unique(m);
colorMap = zeros(length(m), 3);
for i = 1:length(moduleIndices)
    colorMap(m==moduleIndices(i), :) = repmat(colorTriplets(i, :), [sum(m==moduleIndices(i)), 1]);
end

% load labels / ROI names
rois = '/home/adamb/eeg_network_pipeline/utils/roiNamesInOrder.mat';
l = load(rois); labels = l.roisShort;

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



