function graphPlot = circleGraphPlot(connMatrix, membership, colorTriplets, varargin)
%% Plotting network connectivity with module-structure in a circle layout
%
% USAGE: graphPlot = circleGraphPlot(connMatrix, 
%                                   membership,
%                                   colorTriplets,
%                                   trimmingThr=0.2, 
%                                   labels={}, 
%                                   drawFlag=1)
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
%               connections. Defaults to [0.2], must be in range
%               [0:0.01:0.9].
% labels          - Cell array of node labels / names. Defaults to 
%               empty array (no labels). 
% drawFlag        - String, one of {'draw', 'nodraw'}. Flag for 
%               displaying the plot (=1) or only returning the plot object 
%               handle (0). Defaults to 1 (display plot).
%
% Outputs:
% graphPlot       - Graph plot object handle
%
%


%% Input checks

% check number of args
if nargin < 3 || nargin > 6
    error(['Function circleGraphPlot requires mandatory input args "connMatrix", '... 
        '"membership" and "colorTriplets", while input args "trimmingThr", ',...
        '"labels" and "drawFlag" are optional!']);
end

% check mandatory inputs
if ~ismatrix(connMatrix) || size(connMatrix, 1) ~=size (connMatrix, 2)
    error('Input arg "connMatrix" should be a square matrix!');
end
if ~isvector(membership) || length(membership) ~= size(connMatrix, 1)
    error(['Input arg "membership" should be a vector with the same ',...
        'length as either dimension of "connMatrix"!']);
end
if ~ismatrix(colorTriplets) || size(colorTriplets, 2) ~= 3
    error(['Input arg "colorTriplets" should be a matrix with three ',...
        'columns, with each row specifying an RGB color!']);
end

% check optional arguments, parse them
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && ismember(length(varargin{v}), [1 2])
            trimmingThr = varargin{v};
            for t = 1: length(trimmingThr)
                if ~ismember(trimmingThr(t), 0:0.01:0.9)
                    error('Optional input arg "trimmingThr" has value(s) outside 0:0.01:0.9!');
                end
            end
        elseif iscell(varargin{v}) && length(varargin{v}) == length(membership)
            labels = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'draw', 'nodraw'})
            drawFlag = varargin{v};
        else
            error(['An input arg could not be parsed as any of "trimmingThr", ',...
                '"labels" or "drawFlag"!']);
        end
    end
end

% defaults
if ~exist('trimmingThr', 'var')
    trimmingThr = 0.2;
end
if ~ exist('labels', 'var')
    labels = {};
end
if ~ exist('drawFlag', 'var')
    drawFlag = 'draw';
end

% further checks, transformations
if strcmp(drawFlag, 'draw')
    drawFlag = 1;
elseif strcmp(drawFlag, 'nodraw')
    drawFlag = 0;
end
if length(trimmingThr) == 1
    doubleTrim = 0;
elseif length(trimmingThr) == 2
    doubleTrim = 1;
end
if isrow(membership)
    membership = membership';
end
if isrow(labels)
    labels = labels';
end

% get number of modules
modNo = length(unique(membership));

% user message
disp([char(10), 'Function circleGraphPlot is called with inputs: ',...
    char(10), 'Connectivity matrix with size ', num2str(size(connMatrix)),...
    char(10), 'Membership (module) vector with size ', num2str(size(membership)),...
    char(10), 'Number of modules: ', num2str(modNo),...
    char(10), 'Edge trimming theshold: ', num2str(trimmingThr),...
    char(10), 'Node label array with size ', num2str(size(labels)),...
    char(10), 'Figure display flag: ', num2str(drawFlag)]);



%% Hard-coded params

% RGB color for between-module edges
baseEdgeColor = [0.5, 0.5, 0.5];
% multiplier for the width of within-module edges
edgeWidthMultip = 10;
% edge line styles for within- and between-module edges
edgeTypes = {'-', 'none'};
% general transparency setting for edges
edgeAlpha = 0.3;
% general node size setting
nodeSize = 12;
% graph plot layout
graphLayout = 'circle';


%% Prepare connectivity, sort colors to modules

% create symmetric adjacency matrix with zeros at diagonal
connMatrix = triu(connMatrix, 1) + triu(connMatrix, 1)';
connMatrix(isnan(connMatrix)) = 0;  % in many cases connectivity matrices contain NaN values

% indices in "membership" for unique modules
moduleIndices = unique(membership);
% sanity check number of supplied colors
if size(colorTriplets, 1) < length(moduleIndices)
    error('There are less colors specified in "colorTriplets" than modules!');
end

% sort colors to nodes
nodeColors = zeros(length(membership), 3);
for i = 1:length(moduleIndices)
    nodeColors(m==moduleIndices(i), :) = repmat(colorTriplets(i, :), [sum(m==moduleIndices(i)), 1]);
end


%% Init graph object, set properties of within- and between-module edges

% create built-in graph object
G = graph(connMatrix, labels, 'upper');
% edge weights in a vector
weights = G.Edges.Weight;
% edge ending nodes in a cell array
nodesPerEdge = G.Edges.EndNodes;

% go through all modules, set different edge properties per module
edgeColors = repmat(baseEdgeColor, [size(weights, 1), 1]);  % preallocate variable for edge colors, filled with base color
edgeWidth = weights;  % at start, edge width is deterimned by connectivity weight
moduleEdges = zeros(size(weights, 1), modNo);  % binary vectors identifying within-module edges (one column per module)
for i = 1:modNo
    moduleNodes = labels(membership == moduleIndices(i));  % node indices (as binary vector) for given module
    moduleEdges(:, i) = ismember(nodesPerEdge(:, 1), moduleNodes) & ismember(nodesPerEdge(:, 2), moduleNodes);  % edge indices (as binary vector) for edges within given module
    edgeColors(logical(moduleEdges(:, i)), :) = repmat(colorTriplets(i, :), [sum(moduleEdges(:, i)), 1]);  % set edge color for current module
    edgeWidth(logical(moduleEdges(:, i))) = weights(logical(moduleEdges(:, i)))*edgeWidthMultip;  % set edge width for current module
end

% identify between-module edges
betweenModEdgeIdx = ~logical(sum(moduleEdges, 2));
% set line styles for within- and between-module edges
edgeStyle = repmat(edgeTypes(1), [size(weights, 1), 1]);
edgeStyle(betweenModEdgeIdx) = repmat(edgeTypes(2), [sum(betweenModEdgeIdx, 1), 1]);


%% Define a subgraph for each module

% each subgraph is stored in a cell array
subGraphs = cell(modNo, 1);

for i = 1:modNo
    moduleNodes = labels(membership == moduleIndices(i));  % node indices (as binary vector) for given module
    subGraphs{i} = subgraph(G, moduleNodes);
end


%% Plot main graph

% enlarge figure to full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

% positioned on the upper half of the figure, centered
subplot('position', [0 0.5 0.5 0.5]);

% base graph plot
H = G.plot('Layout', graphLayout,... 
    'LineWidth', edgeWidth,... 
    'EdgeColor', edgeColors,... 
    'EdgeAlpha', edgeAlpha,... 
    'NodeColor', nodeColors,...
    'MarkerSize', nodeSize,...
    'LineStyle', edgeStyle);


%% Plot sub graphs

for s = 1:length(subGraphs)
    
    subplot('position', [(1/modNo)*(s-1), 0, 1/modNo, 1/modNo]);
    
    subGraphs{s}.plot('Layout', graphLayout,... 
    'LineWidth', subGraphs{s}.Edges.Weight*edgeWidthMultip,... 
    'EdgeColor', colorTriplets(s, :),... 
    'EdgeAlpha', edgeAlpha,... 
    'NodeColor', colorTriplets(s, :),...
    'MarkerSize', nodeSize,...
    'LineStyle', edgeStyle{1});

end


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



