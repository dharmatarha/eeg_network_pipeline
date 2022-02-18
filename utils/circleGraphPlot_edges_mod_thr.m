function mainFig = circleGraphPlot_edges_mod_thr(edgeContr, edgeStrength, edgeMembership, colorTriplets, varargin)
%% Plotting network connectivity with module-structure in a circle layout
% Version highlighting given edge sets
%
% USAGE: mainFig = circleGraphPlot_edges(connMatrix, 
%                                       edgeMembership,
%                                       colorTriplets,
%                                       group2color=[],
%                                       trimmingThr=0.2, 
%                                       labels={}, 
%                                       figTitle=[];
%                                       drawFlag='draw')
% 
% Creates a circle plot for the supplied network (graph) highlighting 
% edge groupings (memberships) with different colors.
%
% Figure details are fine-tuned for undirected EEG connectivity data. 
% 
% Some finer details (e.g. lobule labels, node grouping into lobules) only
% work with a specific set of labels (62 anatomic ROIs with specific
% order), otherwise they are avoided.
% 
% Specific size and position settings are for 23" screen with a resolution
% 1920 x 1080, no guarantee that figures look decent with anything else.
%
% FOR UNDIRECTED GRAPH!
%
% Mandatory inputs:
% connMatrix      - Numeric matrix containing connectivity (adjacency) 
%               values (square matrix). Only upper triangle is used for 
%               graph construction. 
% edgeMembership  - Numeric matrix containing the membership of each
%               edge in the graph. Only upper triangle is used for 
%               graph construction. Upper triangle values must be positive 
%               integers or zero. 
%               The maximum number of unique values is determined by the 
%               number of colors specified by input arg "colorTriplets". 
%               Zero is treated as identifier for "background", i.e., 
%               not-highlighted edges which are left gray.
% colorTriplets   - Numeric matrix with three columns, each row defines an
%               RGB triplet for edge membership coloring. 
% 
% Optional inputs:
% group2color     - Numeric matrix with two columns. Contains
%               edge group-to-color assignments. The first color contains
%               edge group identifiers, the second corresponding row 
%               numbers of the input arg "colorTriplets". Defaults to
%               empty, in which case edge group-to-color assignment is 
%               automatic and based on ascending corresponding numbers 
%               (e.g., lowest edge group  index is assinged to first row 
%               of "colorTriplets", and so on).
% trimmingThr     - One- or two-element vector containing threshold(s) 
%               for trimming (deleting) weak connections before plotting.
%               If trimmingThr is only one value, the same threshold is
%               applied to all connections. If two values, the first one is
%               applied to highlighted edges, the second to all other edges. 
%               Defaults to [0.2], value(s) must be in range [0:0.001:0.9].
% labels          - Cell array of node labels / names. Defaults to a cell
%               array of numbers {'1', '2', ...}. 
% figTitle        - Char array, displayed as title on the figures. Defaults
%               to [] (= no title).
% drawFlag        - Char array, one of {'draw', 'nodraw'}. Flag for 
%               displaying the plot (='draw') or only returning the plot 
%               object handle ('nodraw'). Defaults to 'draw' (display plot).
%
% Outputs:
% mainFig       - Figure handle for the plot depicting the whole network
%               and highlighting modules with colors.
%
%


%% Input checks

% check number of args
if ~ismember(nargin, 3:8)
    error(['Function circleGraphPlot_edges requires mandatory input args "connMatrix", '... 
        '"edgeMembership" and "colorTriplets", while input args "group2color", "trimmingThr", ',...
        '"labels", "figTitle" and "drawFlag" are optional!']);
end

connMatrix = edgeContr;
% check mandatory inputs
if ~isnumeric(connMatrix) || ~ismatrix(connMatrix) || size(connMatrix, 1) ~=size (connMatrix, 2)
    error('Input arg "connMatrix" should be a numeric square matrix!');
end
if ~isnumeric(edgeMembership) || ~ismatrix(edgeMembership) || ~isequal(size(edgeMembership), size(connMatrix))
    error(['Input arg "edgeMembership" should be a numeric matrix with the same ',...
        'size as "connMatrix"!']);
end
if ~isnumeric(colorTriplets) || ~ismatrix(colorTriplets) || size(colorTriplets, 2) ~= 3
    error(['Input arg "colorTriplets" should be a matrix with three ',...
        'columns, with each row specifying an RGB color!']);
end

% check optional arguments, parse them
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && ismatrix(varargin{v}) && size(varargin{v}, 2)==2 && ~exist('group2color', 'var')
            group2color = varargin{v}; 
        elseif isnumeric(varargin{v}) && ismember(length(varargin{v}), [1 2]) && ~exist('trimmingThr', 'var')
            trimmingThr = varargin{v};
            for t = 1: length(trimmingThr)
                if ~ismember(trimmingThr(t), 0:0.0001:0.9)
                    error('Optional input arg "trimmingThr" has value(s) outside 0:0.001:0.9!');
                end
            end
        elseif iscell(varargin{v}) && length(varargin{v}) == length(edgeMembership) && ~exist('labels', 'var')
            labels = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'draw', 'nodraw'}) && ~exist('drawFlag', 'var')
            drawFlag = varargin{v};
        elseif ischar(varargin{v}) && ~ismember(varargin{v}, {'draw', 'nodraw'}) && ~exist('figTitle', 'var')
            figTitle = varargin{v};            
        else
            error(['An input arg could not be parsed as any of "group2color", "trimmingThr", ',...
                '"labels", "drawFlag" or "figTitle"!']);
        end
    end
end
   
% defaults
if ~exist('node2color', 'var')
    node2color = [];
end
if ~exist('group2color', 'var')
    group2color = [];
end
if ~exist('trimmingThr', 'var')
    trimmingThr = 0.2;
end
if ~ exist('labels', 'var')
    labels = cellstr(num2str([1:length(edgeMembership)]'));  % {'1', '2', '3', ...}
end
if ~ exist('figTitle', 'var')
    figTitle = [];
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
if isrow(labels)
    labels = labels';
end

% get number of edge groupings
groupIndices = unique(edgeMembership(triu(true(size(edgeMembership, 1)), 1)));  % take only values from upper triangle
groupIndices(groupIndices==0) = [];  % zero identifies background / not-highlighted edges, does not count as a specific group of edges
groupNo = length(groupIndices);

% sanity checks on number of supplied colors
if size(colorTriplets, 1) < length(groupIndices)
    error('There are less colors specified in "colorTriplets" than edge groups!');
end
if ~isempty(group2color)
    if ~isempty(setxor(groupIndices, group2color(:,1)))
        error('Edge group indices in "group2color" do not match completely the ones in "edgeMembership"!');
    end
    if any(~ismember(group2color(:, 2), 1:size(colorTriplets, 1)))
        error('At least one color reference in "group2color" is out of bounds!');
    end
end

% user message
disp([char(10), 'Function circleGraphPlot is called with inputs: ',...
    char(10), 'Connectivity matrix with size ', num2str(size(connMatrix)),...
    char(10), 'Edge grouping (membership) matrix with size ', num2str(size(edgeMembership)),...
    char(10), '(Number of groups: ', num2str(groupNo), ')',...
    char(10), 'Edge trimming theshold: ', num2str(trimmingThr),...
    char(10), 'Node label array with size ', num2str(size(labels)),...
    char(10), 'Figure display flag: ', num2str(drawFlag)]);


%% Compare labels to the specific anatomic ROI set 
% if they match, the function can highlight the ROI structure automatically

% matching is performed by roiLabelMatching
[equalFlag, ~, ~] = roiLabelMatching(labels);

% set flag if supplied ROI/node labels is equal to either one of the expected sets
if equalFlag
    lobuleFlag = 1;
    disp([char(10), 'Supplied labels let us highlight the lobules with additional lines and annotations, will do so']);
else
    lobuleFlag = 0;
    disp([char(10), 'Cannot apply lobule-highlighting for given ROI/node label set']);
end


%% Hard-coded params

% RGB color for not-highlighted (background) edges
baseEdgeColor = [0.5, 0.5, 0.5];
% base edge width range - we map the supplied data to this range 
% irrespective of actual weights 
baseEdgeWidthRange = [0.1, 0.3];
% multiplier for the width of highlighted edges
highlEdgeWidthMultip = 2;
% edge line styles for highlighted and not-highlighted edges
%edgeTypes = {'-', 'none'};
edgeTypes = {'-', '-'};
% general transparency setting for edges
edgeAlpha = 0.3;
% general node size setting
nodeSize = 10;
% graph plot layout
graphMainLayout = 'circle';
% figure (gcf) background color
gcfColor = [1 1 1];
% axes (gca) color in subplots
gcaLinesColor = [1 1 1];
% figure position and size in normalized units 
gcfMainPos = [0.25, 0, 0.5, 1];
% axes position relative to figure for main figure
gcaPosInFig = [0.05, 0.05, 0.9, 0.9];
% figure title texts
mainFigTitle = ['Graph with highlighted edge groupings. ', figTitle];

% properties for text box displaying trimming info
if ~doubleTrim
    trimmingText = ['Edges with weight > ', num2str(trimmingThr), ' are depicted'];
elseif doubleTrim
    trimmingText = ['Edges with weight > ', num2str(trimmingThr(1)),...
        ' and > ', num2str(trimmingThr(2)),...
        ' (for highlighted and not-highlighted edges, respectively) are depicted'];
end
trimmingBoxPos = [0.01, 0.01, 0.4, 0.03];

% if ROIs are grouped into lobules, we draw lines and annotations, set
% their params
if lobuleFlag
    % params for drawing lines between lobules
    % angles of lines in radians
    lineAngles = [0.24, 0.665, pi/2, pi-0.665, pi-0.24, pi+0.55,... 
        pi+1.15, pi*3/2, 2*pi-0.55, 2*pi-1.15];
    lineL = 1.75;  % line length
    lineRatio = 0.60;  % lines are not drawn completely, only their ends outside the main graph, lineRatio controls how much is invisible
    xL = cos(lineAngles)*lineL;  % get end point coordinates from angles and length
    yL = sin(lineAngles)*lineL;
    xL = [xL*lineRatio; xL];  % add start points based on lineRatio
    yL = [yL*lineRatio; yL];
    % line properties
    lineWidth = 1.5;
    lineColor = [0 0 0];
    lineStyle = '--';

    % annotation / textbox properties
    typeA = 'textbox';
    % positions: one row per lobule label
    posA = [0.39, 0.06, 0.1, 0.1;
        0.20, 0.15, 0.1, 0.1;
        0.09, 0.35, 0.1, 0.1;
        0.08, 0.57, 0.1, 0.1;
        0.30, 0.77, 0.1, 0.1; 
        0.62, 0.77, 0.1, 0.1;
        0.85, 0.57, 0.1, 0.1;
        0.84, 0.35, 0.1, 0.1;
        0.73, 0.15, 0.1, 0.1;
        0.53, 0.06, 0.1, 0.1];
    % lobule labels
    textA = {'L Occipital', 'L Parietal', 'L Temporal', 'L Cingulate', 'L Frontal',... 
        'R Frontal', 'R Cingulate', 'R Temporal', 'R Parietal', 'R Occipital'};
    
end


%% Prepare connectivity, sort colors to edges

% create symmetric adjacency matrix with zeros at diagonal
connMatrix = triu(connMatrix, 1) + triu(connMatrix, 1)';
connMatrix(isnan(connMatrix)) = 0;  % in many cases connectivity matrices contain NaN values

% sort colors to nodes
nodeColors = repmat(baseEdgeColor, [size(connMatrix, 1), 1]);  % preallocate with base color
if ~isempty(node2color)
    % assignment is based on "node2color", go through all colors used
    tmpColors = unique(node2color(:, 2));
    for i = 1:length(tmpColors)
        currentColorIdx = tmpColors(i);
        nodeColors(node2color(: ,2)==currentColorIdx, :) = repmat(colorTriplets(currentColorIdx, :), [sum(node2color(: ,2)==currentColorIdx), 1]);
    end
end


%% Init graph object, set properties of within- and between-module edges

% create built-in graph object
if isempty(labels)
    G = graph(connMatrix, 'upper');
else
    G = graph(connMatrix, labels, 'upper');
end
% edge weights in a vector
weights = G.Edges.Weight;
% edge ending nodes in a cell array
nodesPerEdge = G.Edges.EndNodes;

% go through all edge groups, set different edge properties per group
edgeColors = repmat(baseEdgeColor, [size(weights, 1), 1]);  % preallocate variable for edge colors, filled with base color

% map connectivity values to basic edge width range specified earlier
% edgeWidth = (weights-min(weights))./(max(weights)-min(weights))*(baseEdgeWidthRange(2)-baseEdgeWidthRange(1))+baseEdgeWidthRange(1); 

% set edge width values based on weights
edgeWidth = weights./mean(weights).*mean(baseEdgeWidthRange);
edgeColor = weights;
edgeStrength = edgeStrength(triu(true(size(edgeStrength, 1)), 1));

% prepare binary vectors identifying edges per group (one column per group)
groupEdges = zeros(size(weights, 1), groupNo); 

% go through group edges
for i = 1:groupNo
    % find label pairs for edges in group "groupIndices(i)"
    currentGroupIdx = groupIndices(i);  % current edge group whose edges we work with
    [iRow, iCol] = ind2sub(size(edgeMembership), find(edgeMembership==currentGroupIdx));  % find row, col indices of edges in edgeMembership who belong to "currentGroupIdx"
    groupLabels = [labels(iRow), labels(iCol)];  % get cell array of node labels for each edge in current edge group
    % compare label pairs of edges in group "currentGroupIdx" to all edges, get binary
    % vector 
    groupEdges(:, i) = ismember(string(nodesPerEdge), string(groupLabels), 'rows');
    % assign RGB color for edges in group "currentGroupIdx"
    edgeGroupColor = colorTriplets(group2color(group2color(:, 1) == groupIndices(i), 2), :);  % get RGB color for current module
    edgeColors(logical(groupEdges(:, i)), :) = repmat(edgeGroupColor, [sum(groupEdges(:, i)), 1]);  % set edge color for edges within current group
    % for highlighted edge groups, apply the corresponding width multiplier
    edgeWidth(logical(groupEdges(:, i))) = edgeWidth(logical(groupEdges(:, i)))*highlEdgeWidthMultip;  % set edge width for current module
    edgeColor(logical(groupEdges(:, i))) = edgeColor(logical(groupEdges(:, i)));
    edgeStrength(logical(groupEdges(:, i))) = edgeStrength(logical(groupEdges(:, i)))*20;
end

% identify not-highlighted edges
backgroundEdgeIdx = ~logical(sum(groupEdges, 2));
% set line styles for highlighted and not-highlighted edges
edgeStyle = repmat(edgeTypes(1), [size(weights, 1), 1]);
edgeStyle(backgroundEdgeIdx) = repmat(edgeTypes(2), [sum(backgroundEdgeIdx, 1), 1]);


%% Trimming edges

% preallocate variable to collect edges to be deleted
edgesToTrim = [];
edgesBelowThr = abs(weights) < trimmingThr(2);
prunedEdges = (edgeStrength == 0);
edgesBelowThr = or(edgesBelowThr, prunedEdges);
edgesToTrim = sort(find(edgesBelowThr));
G = G.rmedge(edgesToTrim);

% delete corresponding rows from edge attribute arrays
if ~isempty(edgesToTrim)
    edgeColors(edgesToTrim, :) = [];
    edgeWidth(edgesToTrim) = [];
    edgeColor(edgesToTrim) = [];
    edgeStrength(edgesToTrim) = [];
    edgeStyle(edgesToTrim) = [];
end


%% Plot main graph

% main graph plot figure
mainFig = figure;

% set figure size and background color
set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);
set(gcf, 'Color', gcfColor);
% figure shown or not?
if ~drawFlag
    set(mainFig, 'Visible', 'off');
end

% graph plot
h = G.plot('Layout', graphMainLayout,... 
    'LineWidth', edgeStrength,... 
    'EdgeColor', edgeColors,... 
    'EdgeAlpha', edgeAlpha,... 
    'NodeColor', nodeColors,...
    'MarkerSize', nodeSize,...
    'LineStyle', edgeStyle);

h.EdgeCData = edgeColor;
colormap('jet');
colorbar;

% title
title(mainFigTitle, 'Interpreter', 'none');

% lines and text boxes highlight the ROIs in each lobule in case of a
% specific ROI set (labels)
if lobuleFlag
    %  draw separating lines between lobules
    line(xL, yL, 'Color', lineColor, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
    % lobule labels as annotations (text boxes)
    for a = 1:length(textA)
        annotation(typeA, posA(a,:), 'String', textA{a}, 'EdgeColor', gcaLinesColor); 
    end
end

% extra annotation displaying trimming info
annotation('textbox', trimmingBoxPos, 'String', trimmingText, 'EdgeColor', gcaLinesColor);

% set axes boundary line colors 
set(gca,'XColor', gcaLinesColor,'YColor', gcaLinesColor);
% set axes position relative to figure
set(gca, 'Position', gcaPosInFig);



return



