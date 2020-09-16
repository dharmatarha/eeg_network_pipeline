function [mainFig, subFig] = circleGraphPlot(connMatrix, membership, colorTriplets, varargin)
%% Plotting network connectivity with module-structure in a circle layout
%
% USAGE: [mainFig, subFig] = circleGraphPlot(connMatrix, 
%                                       membership,
%                                       colorTriplets,
%                                       mod2color=[],
%                                       trimmingThr=0.2, 
%                                       labels={}, 
%                                       figTitle=[];
%                                       drawFlag='draw')
% 
% Creates two figures for the supplied network and its modules.
%
% The first figure is a circle plot of the whole graph, highlighting
% both the module membership of each node and within-module edges with 
% different colors. 
% The second figure depicts the separate modules as subplots. 
% Module layout calculation happens after edge trimming!! Switch the
% positions of the trimming and subgraph extraction blocks if want to
% preserve all edges for module layouts (and subplots).
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
% membership      - Numeric vector containing the module membership of each
%               node in the graph. Values must be positive integers or zero. 
%               The maximum number of unique values is determined by the 
%               number of colors specified by input arg "colorTriplets". 
%               Zero is treated as a valid module identifier.
% colorTriplets   - Numeric matrix with three columns, each row defines an
%               RGB triplet for module (and within-module edge) coloring. 
% 
% Optional inputs:
% mod2color       - Numeric matrix with two columns. Contains
%               module-to-color assignments. The first color contains
%               module identifiers, the second corresponding row numbers of
%               the input arg "colorTriplets". Defaults to empty, in which
%               case module-to-color assignment is automatic and based on
%               ascending corresponding numbers (e.g., lowest module indice is
%               assinged to first row of "colorTriplets", and so on).
% trimmingThr     - One- or two-element vector containing threshold(s) 
%               for trimming (deleting) weak connections before plotting.
%               If trimmingThr is only one value, the same threshold is
%               applied to all connections. If two values, the first one is
%               applied to within-module, the second to between-module
%               connections. Defaults to [0.2], value(s) must be in range
%               [0:0.001:0.9].
% labels          - Cell array of node labels / names. Defaults to a cell
%               array of numbers {'1', '2', ...}. 
% figTitle        - String, displayed as title on the figures. Defaults to
%               [] (none).
% drawFlag        - String, one of {'draw', 'nodraw'}. Flag for 
%               displaying the plot (='draw') or only returning the plot object 
%               handle ('nodraw'). Defaults to 'draw' (display plot).
%
% Outputs:
% mainFig       - Figure handle for the plot depicting the whole network
%               and highlighting modules with colors.
% subFig        - Figure handle for plot with the modules as subplots.
%
%


%% Input checks

% check number of args
if ~ismember(nargin, 3:8)
    error(['Function circleGraphPlot requires mandatory input args "connMatrix", '... 
        '"membership" and "colorTriplets", while input args "mod2color", "trimmingThr", ',...
        '"labels", "figTitle" and "drawFlag" are optional!']);
end

% check mandatory inputs
if ~ismatrix(connMatrix) || size(connMatrix, 1) ~=size (connMatrix, 2)
    error('Input arg "connMatrix" should be a square matrix!');
end
if ~isvector(membership) || length(membership) ~= size(connMatrix, 1)
    error(['Input arg "membership" should be a vector with the same ',...
        'length as either dimension of "connMatrix"!']);
end
if length(unique(membership)) > size(colorTriplets, 1)
    error(['Input arg "membership" contains more unique module ids than ',...
        'the number of colors coded by "colorTriplets"!']);
end
if ~ismatrix(colorTriplets) || size(colorTriplets, 2) ~= 3
    error(['Input arg "colorTriplets" should be a matrix with three ',...
        'columns, with each row specifying an RGB color!']);
end

% check optional arguments, parse them
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && ismatrix(varargin{v}) && size(varargin{v}, 2)==2 && ~exist('mod2color', 'var')
            mod2color = varargin{v}; 
        elseif isnumeric(varargin{v}) && ismember(length(varargin{v}), [1 2]) && ~exist('trimmingThr', 'var')
            trimmingThr = varargin{v};
            for t = 1: length(trimmingThr)
                if ~ismember(trimmingThr(t), 0:0.001:0.9)
                    error('Optional input arg "trimmingThr" has value(s) outside 0:0.001:0.9!');
                end
            end
        elseif iscell(varargin{v}) && length(varargin{v}) == length(membership) && ~exist('labels', 'var')
            labels = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'draw', 'nodraw'}) && ~exist('drawFlag', 'var')
            drawFlag = varargin{v};
        elseif ischar(varargin{v}) && ~ismember(varargin{v}, {'draw', 'nodraw'}) && ~exist('figTitle', 'var')
            figTitle = varargin{v};            
        else
            error(['An input arg could not be parsed as any of "mod2color", "trimmingThr", ',...
                '"labels", "drawFlag" or "figTitle"!']);
        end
    end
end
   
% defaults
if ~exist('mod2color', 'var')
    mod2color = [];
end
if ~exist('trimmingThr', 'var')
    trimmingThr = 0.2;
end
if ~ exist('labels', 'var')
    labels = cellstr(num2str([1:length(membership)]'));  % {'1', '2', '3', ...}
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
if isrow(membership)
    membership = membership';
end
if isrow(labels)
    labels = labels';
end
if ~isempty(mod2color)
    if ~isempty(setxor(unique(membership), mod2color(:,1)))
        error('Module indices in "mod2color" do not match completely the ones in "membership"!');
    end
    if any(~ismember(mod2color(:, 2), 1:size(colorTriplets, 1)))
        error('At least one color reference in "mod2color" is out of bounds!');
    end
end

% get number of modules
modNo = length(unique(membership));

% user message
disp([char(10), 'Function circleGraphPlot is called with inputs: ',...
    char(10), 'Connectivity matrix with size ', num2str(size(connMatrix)),...
    char(10), 'Membership (module) vector with size ', num2str(size(membership)),...
    char(10), '(Number of modules: ', num2str(modNo), ')',...
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

% RGB color for between-module edges
baseEdgeColor = [0.5, 0.5, 0.5];
% base edge width range - we map the supplied data to this range 
% irrespective of actual weights 
baseEdgeWidthRange = [0.1, 4];
% multiplier for the width of within-module edges
withinEdgeWidthMultip = 2;
% edge line styles for within- and between-module edges
%edgeTypes = {'-', 'none'};
edgeTypes = {'-', '-'};
% general transparency setting for edges
edgeAlpha = 0.3;
% general node size setting
nodeSize = 10;
% graph plot layout
graphMainLayout = 'circle';
% graphSubLayout = 'subspace';  % consider 'force' with 'WeightEffect' set to 'inverse'
graphSubLayout = 'force';
% figure (gcf) background color
gcfColor = [1 1 1];
% axes (gca) color in subplots
gcaLinesColor = [1 1 1];
% figure position and size in normalized units 
gcfMainPos = [0.25, 0, 0.5, 1];
gcfSubPos = [0, 0, 1, 1];
% axes position relative to figure for main figure
gcaPosInFig = [0.05, 0.05, 0.9, 0.9];
% figure title texts
mainFigTitle = ['Full module structure. ', figTitle];
subFigTitle = ['Modules separately. ', figTitle];

% properties for text box displaying trimming info
if ~doubleTrim
    trimmingText = ['Edges with weight > ', num2str(trimmingThr), ' are depicted'];
elseif doubleTrim
    trimmingText = ['Edges with weight > ', num2str(trimmingThr(1)),...
        ' and > ', num2str(trimmingThr(2)),...
        ' (for within- and between-module edges, respectively) are depicted'];
end
trimmingBoxPos = [0.01, 0.01, 0.4, 0.03];

% for subFig that contains the module-level graph plots, we define subplot
% positions depending on module number
% modules are plotted into subplots, the positions of subplots depend on
% the number of modules
subPlotPos = subPlotPositions(modNo);

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
nodeColors = zeros(length(membership), 3);  % preallocate
% if the input arg "mod2color" was not supplied, we assign them in
% ascending order
if isempty(mod2color)
    for i = 1:length(moduleIndices)
        nodeColors(membership==moduleIndices(i), :) = repmat(colorTriplets(i, :), [sum(membership==moduleIndices(i)), 1]);
    end
% else assignment is based on "mod2color"
else
    for i = 1:size(mod2color, 1)
        nodeColors(membership==mod2color(i,1), :) = repmat(colorTriplets(mod2color(i, 2), :), [sum(membership==mod2color(i,1)), 1]);
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

% go through all modules, set different edge properties per module
edgeColors = repmat(baseEdgeColor, [size(weights, 1), 1]);  % preallocate variable for edge colors, filled with base color
% map connectivity values to basic edge width range specified earlier
% edgeWidth = (weights-min(weights))./(max(weights)-min(weights))*(baseEdgeWidthRange(2)-baseEdgeWidthRange(1))+baseEdgeWidthRange(1); 
edgeWidth = weights./mean(weights).*mean(baseEdgeWidthRange);
moduleEdges = zeros(size(weights, 1), modNo);  % binary vectors identifying within-module edges (one column per module)
for i = 1:modNo
    moduleNodes = labels(membership == moduleIndices(i));  % node indices (as binary vector) for given module
    moduleEdges(:, i) = ismember(nodesPerEdge(:, 1), moduleNodes) & ismember(nodesPerEdge(:, 2), moduleNodes);  % edge indices (as binary vector) for edges within given module
    moduleEdgeColor = colorTriplets(mod2color(mod2color(:, 1) == moduleIndices(i), 2), :);  % get RGB color for current module
    edgeColors(logical(moduleEdges(:, i)), :) = repmat(moduleEdgeColor, [sum(moduleEdges(:, i)), 1]);  % set edge color for edges within current module
%     edgeColors(logical(moduleEdges(:, i)), :) = repmat(colorTriplets(i, :), [sum(moduleEdges(:, i)), 1]);  % set edge color for current module
%     edgeWidth(logical(moduleEdges(:, i))) = weights(logical(moduleEdges(:, i)))*withinEdgeWidthMultip;  % set edge width for current module
    edgeWidth(logical(moduleEdges(:, i))) = edgeWidth(logical(moduleEdges(:, i)))*withinEdgeWidthMultip;  % set edge width for current module
end

% identify between-module edges
betweenModEdgeIdx = ~logical(sum(moduleEdges, 2));
% set line styles for within- and between-module edges
edgeStyle = repmat(edgeTypes(1), [size(weights, 1), 1]);
edgeStyle(betweenModEdgeIdx) = repmat(edgeTypes(2), [sum(betweenModEdgeIdx, 1), 1]);


%% Trimming edges

% preallocate variable to collect edges to be deleted
edgesToTrim = [];
% if same threshold applies to within- and between-module edges
if ~doubleTrim && trimmingThr ~= 0
    edgesToTrim = find(weights < trimmingThr);  % graph.rmedge does not work with logical indexing, requires numeric edge indices
    G = G.rmedge(edgesToTrim);
% if there are different threshold values for within- and between-module
% edges
elseif doubleTrim && any(trimmingThr ~= 0)
    % within-module edges
    withinEdges = logical(sum(moduleEdges, 2));
    edgesBelowThr = weights < trimmingThr(1);
    withinEdgesToTrim = find(withinEdges & edgesBelowThr);
    % between-module edges
    edgesBelowThr = weights < trimmingThr(2);
    betweenEdgesToTrim = find(~withinEdges & edgesBelowThr);
    % summarize into one vector of indices and remove edges
    edgesToTrim = sort([withinEdgesToTrim; betweenEdgesToTrim]);
    G = G.rmedge(edgesToTrim);
end

% delete corresponding rows from edge attribute arrays
if ~isempty(edgesToTrim)
    edgeColors(edgesToTrim, :) = [];
    edgeWidth(edgesToTrim) = [];
    edgeStyle(edgesToTrim) = [];
end


%% Define a subgraph for each module
   
% subgraphs are defined before removing any edges

% subgraphs are stored in a cell array
subGraphs = cell(modNo, 1);

for i = 1:modNo
    moduleNodes = labels(membership == moduleIndices(i));  % node indices (as binary vector) for given module
    subGraphs{i} = subgraph(G, moduleNodes);
end


%% Plot main graph

% main graph plot figure
mainFig = figure;

% set figure size and background color
set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);
set(gcf, 'Color', gcfColor);

% graph plot
G.plot('Layout', graphMainLayout,... 
    'LineWidth', edgeWidth,... 
    'EdgeColor', edgeColors,... 
    'EdgeAlpha', edgeAlpha,... 
    'NodeColor', nodeColors,...
    'MarkerSize', nodeSize,...
    'LineStyle', edgeStyle);

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


%% Plot sub graphs

% module-level graphs figure
subFig = figure;

% set figure size and background color
set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfSubPos);
set(gcf, 'Color', gcfColor);

for s = 1:modNo
    
    % positioned in lower half, one next to the other
    subplot('position', subPlotPos(s, :));
    
    % edge and node color for given module
    subModColor = colorTriplets(mod2color(mod2color(:,1)==moduleIndices(s), 2), :);
    
    % module subplot
    subGraphs{s}.plot('Layout', graphSubLayout,... 
        'LineWidth', subGraphs{s}.Edges.Weight*withinEdgeWidthMultip,... 
        'EdgeColor', subModColor,... 
        'EdgeAlpha', edgeAlpha,... 
        'NodeColor', subModColor,...
        'MarkerSize', nodeSize,...
        'LineStyle', edgeTypes{1});

    % set axes boundary line colors
    set(gca,'XColor', gcaLinesColor,'YColor', gcaLinesColor);

end

% title
suptitle(subFigTitle);
% extra annotation displaying trimming info
annotation('textbox', trimmingBoxPos, 'String', trimmingText, 'EdgeColor', gcaLinesColor);


return



%% Helper function storing details for subplot positions for subFig

function subPlotPos = subPlotPositions(modNo)

switch modNo 
    case 2  
        subPlotPos = [0.05, 0.05, 0.35, 0.7;
             0.55, 0.05, 0.35, 0.7];       
    case 3
         subPlotPos = [0.05, 0.55, 0.35, 0.35;
             0.55, 0.55, 0.35, 0.35;
             0.25, 0.05, 0.45, 0.35];          
    case 4 
         subPlotPos = [0.05, 0.55, 0.35, 0.35;
             0.55, 0.55, 0.35, 0.35;
             0.05, 0.05, 0.35, 0.35;
             0.55, 0.05, 0.35, 0.35];       
    case 5
         subPlotPos = [0.05, 0.55, 0.25, 0.25;
             0.32, 0.55, 0.25, 0.25;
             0.65, 0.55, 0.25, 0.25;
             0.05, 0.05, 0.35, 0.25;
             0.55, 0.05, 0.35, 0.25];     
    case 6
         subPlotPos = [0.05, 0.55, 0.25, 0.25;
             0.32, 0.55, 0.25, 0.25;
             0.65, 0.55, 0.25, 0.25;
             0.05, 0.05, 0.25, 0.25;
             0.32, 0.05, 0.25, 0.25;
             0.65, 0.05, 0.25, 0.25;];   
    case 7
         subPlotPos = [0.05, 0.66, 0.25, 0.25;
             0.32, 0.66, 0.25, 0.25;
             0.65, 0.66, 0.25, 0.25;
             0.05, 0.33, 0.25, 0.25;
             0.32, 0.33, 0.25, 0.25;
             0.65, 0.33, 0.25, 0.25;
             0.35, 0.05, 0.25, 0.25];      
    case 8
         subPlotPos = [0.05, 0.66, 0.25, 0.25;
             0.32, 0.66, 0.25, 0.25;
             0.65, 0.66, 0.25, 0.25;
             0.05, 0.33, 0.25, 0.25;
             0.32, 0.33, 0.25, 0.25;
             0.65, 0.33, 0.25, 0.25;
             0.20, 0.05, 0.25, 0.25;
             0.60, 0.05, 0.25, 0.25];       
    case 9
         subPlotPos = [0.05, 0.66, 0.25, 0.25;
             0.32, 0.66, 0.25, 0.25;
             0.65, 0.66, 0.25, 0.25;
             0.05, 0.33, 0.25, 0.25;
             0.32, 0.33, 0.25, 0.25;
             0.65, 0.33, 0.25, 0.25;
             0.05, 0.05, 0.25, 0.25;
             0.32, 0.05, 0.25, 0.25;
             0.65, 0.05, 0.25, 0.25];  
    case 10
         subPlotPos = [0.02, 0.66, 0.20, 0.20;
             0.25, 0.66, 0.20, 0.20;
             0.50, 0.66, 0.20, 0.20;
             0.75, 0.66, 0.20, 0.20;
             0.02, 0.33, 0.20, 0.20;
             0.25, 0.33, 0.20, 0.20;
             0.50, 0.33, 0.20, 0.20;
             0.75, 0.33, 0.20, 0.20;
             0.20, 0.05, 0.20, 0.20;
             0.60, 0.05, 0.20, 0.20];
    case 11
         subPlotPos = [0.02, 0.66, 0.20, 0.20;
             0.25, 0.66, 0.20, 0.20;
             0.50, 0.66, 0.20, 0.20;
             0.75, 0.66, 0.20, 0.20;
             0.02, 0.33, 0.20, 0.20;
             0.25, 0.33, 0.20, 0.20;
             0.50, 0.33, 0.20, 0.20;
             0.75, 0.33, 0.20, 0.20;
             0.10, 0.05, 0.20, 0.20;
             0.40, 0.05, 0.20, 0.20;
             0.70, 0.05, 0.20, 0.20];
    case 12
         subPlotPos = [0.02, 0.66, 0.20, 0.20;
             0.25, 0.66, 0.20, 0.20;
             0.50, 0.66, 0.20, 0.20;
             0.75, 0.66, 0.20, 0.20;
             0.02, 0.33, 0.20, 0.20;
             0.25, 0.33, 0.20, 0.20;
             0.50, 0.33, 0.20, 0.20;
             0.75, 0.33, 0.20, 0.20;
             0.02, 0.05, 0.20, 0.20;
             0.25, 0.05, 0.20, 0.20;
             0.50, 0.05, 0.20, 0.20;
             0.75, 0.05, 0.20, 0.20];     
    case 13
         subPlotPos = [0.02, 0.66, 0.20, 0.20;
             0.25, 0.66, 0.20, 0.20;
             0.50, 0.66, 0.20, 0.20;
             0.75, 0.66, 0.20, 0.20;
             0.02, 0.33, 0.20, 0.20;
             0.25, 0.33, 0.20, 0.20;
             0.50, 0.33, 0.20, 0.20;
             0.75, 0.33, 0.20, 0.20;
             0.02, 0.05, 0.17, 0.17;
             0.20, 0.05, 0.17, 0.17;
             0.40, 0.05, 0.17, 0.17;
             0.60, 0.05, 0.17, 0.17;
             0.80, 0.05, 0.17, 0.17];   
     case 14
        subPlotPos = [0.02, 0.66, 0.20, 0.20;
             0.25, 0.66, 0.20, 0.20;
             0.50, 0.66, 0.20, 0.20;
             0.75, 0.66, 0.20, 0.20;
             0.02, 0.33, 0.17, 0.17;
             0.20, 0.33, 0.17, 0.17;
             0.40, 0.33, 0.17, 0.17;
             0.60, 0.33, 0.17, 0.17;
             0.80, 0.33, 0.17, 0.17;
             0.02, 0.05, 0.17, 0.17;
             0.20, 0.05, 0.17, 0.17;
             0.40, 0.05, 0.17, 0.17;
             0.60, 0.05, 0.17, 0.17;
             0.80, 0.05, 0.17, 0.17]; 
     case 15
        subPlotPos = [0.02, 0.66, 0.17, 0.17;
             0.20, 0.66, 0.17, 0.17;
             0.40, 0.66, 0.17, 0.17;
             0.60, 0.66, 0.17, 0.17;
             0.80, 0.66, 0.17, 0.17;
             0.02, 0.33, 0.17, 0.17;
             0.20, 0.33, 0.17, 0.17;
             0.40, 0.33, 0.17, 0.17;
             0.60, 0.33, 0.17, 0.17;
             0.80, 0.33, 0.17, 0.17;
             0.02, 0.05, 0.17, 0.17;
             0.20, 0.05, 0.17, 0.17;
             0.40, 0.05, 0.17, 0.17;
             0.60, 0.05, 0.17, 0.17;
             0.80, 0.05, 0.17, 0.17];      
     case 16
        subPlotPos = [0.02, 0.73, 0.20, 0.20;
             0.25, 0.73, 0.20, 0.20;
             0.50, 0.73, 0.20, 0.20;
             0.75, 0.73, 0.20, 0.20;
             0.02, 0.48, 0.20, 0.20;
             0.25, 0.48, 0.20, 0.20;
             0.50, 0.48, 0.20, 0.20;
             0.75, 0.48, 0.20, 0.20;
             0.02, 0.24, 0.20, 0.20;
             0.25, 0.24, 0.20, 0.20;
             0.50, 0.24, 0.20, 0.20;
             0.75, 0.24, 0.20, 0.20;
             0.02, 0.04, 0.20, 0.20;
             0.25, 0.04, 0.20, 0.20;
             0.50, 0.04, 0.20, 0.20;
             0.75, 0.04, 0.20, 0.20];      
     case 17
        subPlotPos = [0.02, 0.73, 0.20, 0.20;
             0.25, 0.73, 0.20, 0.20;
             0.50, 0.73, 0.20, 0.20;
             0.75, 0.73, 0.20, 0.20;
             0.02, 0.48, 0.20, 0.20;
             0.25, 0.48, 0.20, 0.20;
             0.50, 0.48, 0.20, 0.20;
             0.75, 0.48, 0.20, 0.20;
             0.02, 0.24, 0.20, 0.20;
             0.25, 0.24, 0.20, 0.20;
             0.50, 0.24, 0.20, 0.20;
             0.75, 0.24, 0.20, 0.20;
             0.02, 0.04, 0.17, 0.17;
             0.20, 0.04, 0.17, 0.17;
             0.40, 0.04, 0.17, 0.17;
             0.60, 0.04, 0.17, 0.17;
             0.80, 0.04, 0.17, 0.17];    
     case 18
        subPlotPos = [0.02, 0.73, 0.20, 0.20;
             0.25, 0.73, 0.20, 0.20;
             0.50, 0.73, 0.20, 0.20;
             0.75, 0.73, 0.20, 0.20;
             0.02, 0.48, 0.20, 0.20;
             0.25, 0.48, 0.20, 0.20;
             0.50, 0.48, 0.20, 0.20;
             0.75, 0.48, 0.20, 0.20;
             0.02, 0.24, 0.17, 0.17;
             0.20, 0.24, 0.17, 0.17;
             0.40, 0.24, 0.17, 0.17;
             0.60, 0.24, 0.17, 0.17;
             0.80, 0.24, 0.17, 0.17;
             0.02, 0.04, 0.17, 0.17;
             0.20, 0.04, 0.17, 0.17;
             0.40, 0.04, 0.17, 0.17;
             0.60, 0.04, 0.17, 0.17;
             0.80, 0.04, 0.17, 0.17];  
     case 19
        subPlotPos = [0.02, 0.73, 0.20, 0.20;
             0.25, 0.73, 0.20, 0.20;
             0.50, 0.73, 0.20, 0.20;
             0.75, 0.73, 0.20, 0.20;
             0.02, 0.48, 0.17, 0.17;
             0.20, 0.48, 0.17, 0.17;
             0.40, 0.48, 0.17, 0.17;
             0.60, 0.48, 0.17, 0.17;
             0.80, 0.48, 0.17, 0.17;
             0.02, 0.24, 0.17, 0.17;
             0.20, 0.24, 0.17, 0.17;
             0.40, 0.24, 0.17, 0.17;
             0.60, 0.24, 0.17, 0.17;
             0.80, 0.24, 0.17, 0.17;
             0.02, 0.04, 0.17, 0.17;
             0.20, 0.04, 0.17, 0.17;
             0.40, 0.04, 0.17, 0.17;
             0.60, 0.04, 0.17, 0.17;
             0.80, 0.04, 0.17, 0.17];      
     case 20
        subPlotPos = [0.02, 0.73, 0.17, 0.17;
             0.20, 0.73, 0.17, 0.17;
             0.40, 0.73, 0.17, 0.17;
             0.60, 0.73, 0.17, 0.17;
             0.80, 0.73, 0.17, 0.17;
             0.02, 0.48, 0.17, 0.17;
             0.20, 0.48, 0.17, 0.17;
             0.40, 0.48, 0.17, 0.17;
             0.60, 0.48, 0.17, 0.17;
             0.80, 0.48, 0.17, 0.17;
             0.02, 0.24, 0.17, 0.17;
             0.20, 0.24, 0.17, 0.17;
             0.40, 0.24, 0.17, 0.17;
             0.60, 0.24, 0.17, 0.17;
             0.80, 0.24, 0.17, 0.17;
             0.02, 0.04, 0.17, 0.17;
             0.20, 0.04, 0.17, 0.17;
             0.40, 0.04, 0.17, 0.17;
             0.60, 0.04, 0.17, 0.17;
             0.80, 0.04, 0.17, 0.17];      
     case 21
        subPlotPos = [0.02, 0.73, 0.17, 0.17;
             0.20, 0.73, 0.17, 0.17;
             0.40, 0.73, 0.17, 0.17;
             0.60, 0.73, 0.17, 0.17;
             0.80, 0.73, 0.17, 0.17;
             0.02, 0.48, 0.17, 0.17;
             0.20, 0.48, 0.17, 0.17;
             0.40, 0.48, 0.17, 0.17;
             0.60, 0.48, 0.17, 0.17;
             0.80, 0.48, 0.17, 0.17;
             0.02, 0.24, 0.17, 0.17;
             0.20, 0.24, 0.17, 0.17;
             0.40, 0.24, 0.17, 0.17;
             0.60, 0.24, 0.17, 0.17;
             0.80, 0.24, 0.17, 0.17;
             0.02, 0.04, 0.15, 0.17;
             0.17, 0.04, 0.15, 0.17;
             0.33, 0.04, 0.15, 0.17;
             0.49, 0.04, 0.15, 0.17;
             0.65, 0.04, 0.15, 0.17;
             0.81, 0.04, 0.15, 0.17]; 
     case 22
        subPlotPos = [0.02, 0.73, 0.17, 0.17;
             0.20, 0.73, 0.17, 0.17;
             0.40, 0.73, 0.17, 0.17;
             0.60, 0.73, 0.17, 0.17;
             0.80, 0.73, 0.17, 0.17;
             0.02, 0.48, 0.17, 0.17;
             0.20, 0.48, 0.17, 0.17;
             0.40, 0.48, 0.17, 0.17;
             0.60, 0.48, 0.17, 0.17;
             0.80, 0.48, 0.17, 0.17;
             0.02, 0.24, 0.15, 0.17;
             0.17, 0.24, 0.15, 0.17;
             0.33, 0.24, 0.15, 0.17;
             0.49, 0.24, 0.15, 0.17;
             0.65, 0.24, 0.15, 0.17;
             0.81, 0.24, 0.15, 0.17;
             0.02, 0.04, 0.15, 0.17;
             0.17, 0.04, 0.15, 0.17;
             0.33, 0.04, 0.15, 0.17;
             0.49, 0.04, 0.15, 0.17;
             0.65, 0.04, 0.15, 0.17;
             0.81, 0.04, 0.15, 0.17];   
     case 23
        subPlotPos = [0.02, 0.73, 0.17, 0.17;
             0.20, 0.73, 0.17, 0.17;
             0.40, 0.73, 0.17, 0.17;
             0.60, 0.73, 0.17, 0.17;
             0.80, 0.73, 0.17, 0.17;
             0.02, 0.48, 0.15, 0.17;
             0.17, 0.48, 0.15, 0.17;
             0.33, 0.48, 0.15, 0.17;
             0.49, 0.48, 0.15, 0.17;
             0.65, 0.48, 0.15, 0.17;
             0.81, 0.48, 0.15, 0.17;
             0.02, 0.24, 0.15, 0.17;
             0.17, 0.24, 0.15, 0.17;
             0.33, 0.24, 0.15, 0.17;
             0.49, 0.24, 0.15, 0.17;
             0.65, 0.24, 0.15, 0.17;
             0.81, 0.24, 0.15, 0.17;
             0.02, 0.04, 0.15, 0.17;
             0.17, 0.04, 0.15, 0.17;
             0.33, 0.04, 0.15, 0.17;
             0.49, 0.04, 0.15, 0.17;
             0.65, 0.04, 0.15, 0.17;
             0.81, 0.04, 0.15, 0.17];   
     case 24
        subPlotPos = [0.02, 0.73, 0.15, 0.17;
             0.17, 0.73, 0.15, 0.17;
             0.33, 0.73, 0.15, 0.17;
             0.49, 0.73, 0.15, 0.17;
             0.65, 0.73, 0.15, 0.17;
             0.81, 0.73, 0.15, 0.17;
             0.02, 0.48, 0.15, 0.17;
             0.17, 0.48, 0.15, 0.17;
             0.33, 0.48, 0.15, 0.17;
             0.49, 0.48, 0.15, 0.17;
             0.65, 0.48, 0.15, 0.17;
             0.81, 0.48, 0.15, 0.17;
             0.02, 0.24, 0.15, 0.17;
             0.17, 0.24, 0.15, 0.17;
             0.33, 0.24, 0.15, 0.17;
             0.49, 0.24, 0.15, 0.17;
             0.65, 0.24, 0.15, 0.17;
             0.81, 0.24, 0.15, 0.17;
             0.02, 0.04, 0.15, 0.17;
             0.17, 0.04, 0.15, 0.17;
             0.33, 0.04, 0.15, 0.17;
             0.49, 0.04, 0.15, 0.17;
             0.65, 0.04, 0.15, 0.17;
             0.81, 0.04, 0.15, 0.17];           
end

return
