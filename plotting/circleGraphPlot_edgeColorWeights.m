function [mainFig, G, cmapColors] = circleGraphPlot_edgeColorWeights(connMatrix, edgeColorWeights, colorMap, varargin)
%% Plotting network connectivity with module-structure in a circle layout
% Version highlighting given edge sets
%
% USAGE: mainFig = circleGraphPlot_edgeColorWeights(connMatrix, 
%                                           edgeColorWeights,
%                                           colorMap,
%                                           trimmingThr=0.2, 
%                                           labels={}, 
%                                           figTitle=[];
%                                           drawFlag='draw')
% 
% Creates a circle plot for the supplied network (graph) highlighting 
% edges based on their weights ("edgeColorWeights") with colors from the 
% defined colormap ("colorMap").
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
% connMatrix    - Numeric matrix containing connectivity (adjacency) 
%               values (square matrix). Only upper triangle is used for 
%               graph construction. 
% edgeColorWeights   - Numeric matrix containing the weight of each
%               edge in the graph used for coloring. Only upper triangle 
%               is used for graph construction. 
%               The values in the matrix will be mapped to the colormap 
%               in "colorMap" for edge-based coloring. Zero is valid value 
%               and is mapped to the colormap as well, use NaN to avoid 
%               coloring of edges.
% colorMap      - Char array, name of valid Matlab colormap (e.g. "jet"). 
% 
% Optional inputs:
% trimmingThr   - Numeric value. threshold for trimming (deleting) weak 
%               connections before plotting. The same threshold is 
%               applied to all connections.
%               Defaults to [0.2], value(s) must be in range [0:0.001:0.9].
% labels        - Cell array of node labels / names. Defaults to a cell
%               array of numbers {'1', '2', ...}. 
% figTitle      - Char array, displayed as title on the figures. Defaults
%               to [] (= no title).
% drawFlag      - Char array, one of {'draw', 'nodraw'}. Flag for 
%               displaying the plot (='draw') or only returning the plot 
%               object handle ('nodraw'). Defaults to 'draw' (display plot).
%
% Outputs:
% mainFig       - Figure handle for the plot depicting the network.
% G             - Matlab's graph object for the network depicted
%


%% Input checks

% check number of args
if ~ismember(nargin, 3:7)
    error(['Function circleGraphPlot_edgeColorWeights requires mandatory input args "connMatrix", '... 
        '"edgeColorWeights" and "colorMap", while input args "trimmingThr", ',...
        '"labels", "figTitle" and "drawFlag" are optional!']);
end

% check mandatory inputs
if ~isnumeric(connMatrix) || ~ismatrix(connMatrix) || size(connMatrix, 1) ~=size (connMatrix, 2)
    error('Input arg "connMatrix" should be a numeric square matrix!');
end
if ~isnumeric(edgeColorWeights) || ~ismatrix(edgeColorWeights) || ~isequal(size(edgeColorWeights), size(connMatrix))
    error(['Input arg "edgeColorWeights" should be a numeric matrix with the same ',...
        'size as "connMatrix"!']);
end
if ~ischar(colorMap)
    error(['Input arg "colorMap" should be a character array, the name of a valid matlab colormap (e.g. "jet")!']);
% check if "colorMap" can be treated as a colormap function
else 
    mycmap = flipud(str2func(colorMap));
    try
        tmp = mycmap(1);
        if numel(tmp)~=3 || any(tmp<0) || any(tmp>1)
            error('Oops, wrong colormap name');
        end
    catch ME
        error('Supplied char array in "colorMap" is not a valid colormap name!');
    end 
end

% check optional arguments, parse them
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && ismember(length(varargin{v}), [1 2]) && ~exist('trimmingThr', 'var')
            trimmingThr = varargin{v};
            for t = 1: length(trimmingThr)
                if ~ismembertol(trimmingThr(t), 0:0.001:0.9)
                    error('Optional input arg "trimmingThr" has value(s) outside 0:0.001:0.9!');
                end
            end
        elseif iscell(varargin{v}) && length(varargin{v}) == size(edgeColorWeights, 1) && ~exist('labels', 'var')
            labels = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'draw', 'nodraw'}) && ~exist('drawFlag', 'var')
            drawFlag = varargin{v};
        elseif ischar(varargin{v}) && ~ismember(varargin{v}, {'draw', 'nodraw'}) && ~exist('figTitle', 'var')
            figTitle = varargin{v};            
        else
            error(['An input arg could not be parsed as any of "trimmingThr", ',...
                '"labels", "drawFlag" or "figTitle"!']);
        end
    end
end
   
% defaults
if ~exist('trimmingThr', 'var')
    trimmingThr = 0.2;
end
if ~ exist('labels', 'var')
    labels = cellstr(num2str([1:size(edgeColorWeights, 1)]'));  % {'1', '2', '3', ...}
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

% user message
disp([char(10), 'Function circleGraphPlot_edgeColorWeights is called with inputs: ',...
    char(10), 'Connectivity matrix with size ', num2str(size(connMatrix)),...
    char(10), 'Edge coloring matrix with size ', num2str(size(edgeColorWeights)),...
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
baseEdgeWidthRange = [2, 6];
% % multiplier for the width of highlighted edges
% highlEdgeWidthMultip = 2;
% edge line styles for highlighted and not-highlighted edges
%edgeTypes = {'-', 'none'};
edgeTypes = {'-', '-'};
% general transparency setting for edges
edgeAlpha = 0.5;
% general node size setting
nodeSize = 10;
% graph plot layout
graphMainLayout = 'circle';
% figure (gcf) background color
gcfColor = [1 1 1];
% axes (gca) color in subplots
gcaLinesColor = [1 1 1];
% figure position and size in normalized units 
gcfMainPos = [0, 0, 0.5, 1];
% axes position relative to figure for main figure
gcaPosInFig = [0.05, 0.05, 0.95, 0.95];
% figure title texts
mainFigTitle = ['Graph with colored edges. ', figTitle];

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
    lineWidth = 2;
    lineColor = [0 0 0];
    lineStyle = '--';

    % annotation / textbox properties
    typeA = 'textbox';
    % positions: one row per lobule label
    posA = [0.39, 0.06, 0.1, 0.1;
        0.20, 0.15, 0.1, 0.1;
        0.09, 0.35, 0.1, 0.1;
        0.08, 0.59, 0.1, 0.1;
        0.32, 0.77, 0.1, 0.1; 
        0.66, 0.77, 0.1, 0.1;
        0.85, 0.59, 0.1, 0.1;
        0.87, 0.35, 0.1, 0.1;
        0.77, 0.15, 0.1, 0.1;
        0.57, 0.06, 0.1, 0.1];
    % lobule labels
    textA = {'L Occipital', 'L Parietal', 'L Temporal', 'L Cingulate', 'L Frontal',... 
        'R Frontal', 'R Cingulate', 'R Temporal', 'R Parietal', 'R Occipital'};
    
end


%% Prepare connectivity, sort colors to edges

% create symmetric adjacency matrix with zeros at diagonal
connMatrix = triu(connMatrix, 1) + triu(connMatrix, 1)';
connMatrix(isnan(connMatrix)) = 0;  % in many cases connectivity matrices contain NaN values

% sort colors to edges
nodeNo = size(edgeColorWeights, 1);
colorWeights = edgeColorWeights;
colorWeights = colorWeights(triu(true(nodeNo),1));  % extract only upper triangle
colorWeights(isnan(colorWeights)) = [];  % delete NaN values
% get colors for the number of values we have
cmapColors = mycmap(numel(colorWeights));  % get an RGB color for each value
% values ordered according to colors in "edgeColors"
colorWeightsSorted = sort(colorWeights, 'descend');

% % match color vectors to edges, in the same format as the connectivity
% % matrix but with three values per entry (matrix nodeNo X nodeNo X RGB)
% edgeColorsMatrix = nan(nodeNo, nodeNo, 3);
% for y = 1:numel(colorWeightsSorted)
%     [i, j] = find(edgeColorWeights==colorWeightsSorted(y));
%     if ~isempty(i)
%         for idxNo = 1: numel(i)
%             edgeColorsMatrix(i(idxNo), j(idxNo), :) = cmapColors(y);
%         end
%     end
% end
        

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

% map connectivity values to basic edge width range specified earlier
edgeWidth = (weights-min(weights))./(max(weights)-min(weights))*(baseEdgeWidthRange(2)-baseEdgeWidthRange(1))+baseEdgeWidthRange(1); 

% % set edge width values based on weights
% edgeWidth = weights./mean(weights).*mean(baseEdgeWidthRange);

% set line styles for highlighted and not-highlighted edges
edgeStyle = repmat(edgeTypes(1), [size(weights, 1), 1]);

% Set edge colors
edgeColors = nan(size(weights, 1), 3);
for i = 1:numel(weights)
    % find the corresponding edgeColorWeights value for current weight
    idxM = connMatrix==weights(i);
    tmpColorW = edgeColorWeights(idxM);
    tmpColorW(isnan(tmpColorW)) = [];
    if numel(tmpColorW)>1
        tmpColorW = tmpColorW(1);
    end
    % get corresponding color
    tmpColor = cmapColors(colorWeightsSorted==tmpColorW, :);
    if size(tmpColor, 1)>1
        tmpColor = tmpColor(1,:);
    end
    % to edgeColors
    edgeColors(i, :) = tmpColor;
end
    
% Set node colors
nodeColors = repmat(baseEdgeColor, [nodeNo, 1]);


%% Trimming edges

% preallocate variable to collect edges to be deleted
edgesToTrim = [];
% if same threshold applies to highlighted and not-highlighted edges
if ~doubleTrim && trimmingThr ~= 0
    edgesToTrim = find(weights < trimmingThr);  % graph.rmedge does not work with logical indexing, requires numeric edge indices
    G = G.rmedge(edgesToTrim);
end
% % if there are different threshold values for highlighted and not-highlighted edges
% elseif doubleTrim && any(trimmingThr ~= 0)
%     % highlighted edges
%     highlEdges = logical(sum(groupEdges, 2));
%     edgesBelowThr = weights < trimmingThr(1);
%     highlEdgesToTrim = find(highlEdges & edgesBelowThr);
%     %  not-highlighted edges
%     edgesBelowThr = weights < trimmingThr(2);
%     backgrEdgesToTrim = find(~highlEdges & edgesBelowThr);
%     % summarize into one vector of indices and remove edges
%     edgesToTrim = sort([highlEdgesToTrim; backgrEdgesToTrim]);
%     G = G.rmedge(edgesToTrim);
% end

% delete corresponding rows from edge attribute arrays
if ~isempty(edgesToTrim)
    edgeColors(edgesToTrim, :) = [];
    edgeWidth(edgesToTrim) = [];
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
G.plot('Layout', graphMainLayout,... 
    'LineWidth', edgeWidth,... 
    'EdgeColor', edgeColors,... 
    'EdgeAlpha', edgeAlpha,... 
    'NodeColor', nodeColors,...
    'MarkerSize', nodeSize,...
    'LineStyle', edgeStyle);

% % title
% title(mainFigTitle, 'Interpreter', 'none');

% lines and text boxes highlight the ROIs in each lobule in case of a
% specific ROI set (labels)
if lobuleFlag
    %  draw separating lines between lobules
    line(xL, yL, 'Color', lineColor, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
    % lobule labels as annotations (text boxes)
    for a = 1:length(textA)
        annotation(typeA, posA(a,:), 'String', textA{a}, 'EdgeColor', gcaLinesColor, 'FontSize', 11); 
    end
end

% % extra annotation displaying trimming info
% annotation('textbox', trimmingBoxPos, 'String', trimmingText, 'EdgeColor', gcaLinesColor);

% set axes boundary line colors 
set(gca,'XColor', gcaLinesColor,'YColor', gcaLinesColor);
% set axes position relative to figure
set(gca, 'Position', gcaPosInFig);
% % set axes font sizes
set(gca, 'FontSize', 18);


return



