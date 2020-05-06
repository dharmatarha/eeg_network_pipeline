function [mainFig, subFig] = circleGraphPlotSimple(connMatrix, varargin)
%% Plotting network connectivity with module-structure in a circle layout
%
% USAGE: [mainFig, subFig] = circleGraphPlot(connMatrix, 
%                                       colorTriplet=[0, 0.447, 0.741],
%                                       trimmingThr=0.2, 
%                                       labels={}, 
%                                       figTitle=[];
%                                       drawFlag=1)
% 
% Creates a circle graph figure of the supplied network.
%
% Figure details are fine-tuned for undirected EEG connectivity data.
% 
% Some finer details (e.g. lobule labels, node grouping into lobules) only
% work with a specific set of labels (62 anatomic ROIs with specific
% order), otherwise they are avoided.
% 
% Specific size and position settings are for 23" screen with a resolution
% 1920 x 1080, no guarantee that figure looks decent with anything else.
%
% FOR UNDIRECTED GRAPH!
%
% Mandatory inputs:
% connMatrix      - Numeric matrix containing connectivity (adjacency) 
%               values (square matrix). Only upper triangle is used for 
%               graph construction. 
% 
% Optional inputs:
% colorTriplet    - Numeric vector of length = 3, defining an RGB color.
%               Defaults to [0, 0.447, 0.741] (default matlab blue). 
% trimmingThr     - Numeric value defining the threshold 
%               for trimming (deleting) weak connections before plotting. 
%               Defaults to [0.2], value(s) must be in range [0:0.01:0.9].
% labels          - Cell array of node labels / names. Defaults to a cell
%               array of numbers {'1', '2', ...}. 
% figTitle        - String, displayed as title on the figures. Defaults to
%               [] (none).
% drawFlag        - String, one of {'draw', 'nodraw'}. Flag for 
%               displaying the plot (='draw') or only returning the plot object 
%               handle ('nodraw'). Defaults to 'draw' (display plot).
%
% Outputs:
% figHandle      - Figure handle for the the whole network circle graph.
%
%



%% Input checks

% check number of args
if nargin < 1 || nargin > 6
    error(['Function circleGraphPlotSimple requires mandatory input arg "connMatrix", '... 
        ', while input args "colorTriplet", "trimmingThr", ',...
        '"labels", "figTitle" and "drawFlag" are optional!']);
end

% check mandatory inputs
if ~ismatrix(connMatrix) || size(connMatrix, 1) ~=size (connMatrix, 2)
    error('Input arg "connMatrix" should be a square matrix!');
end
if ~ismatrix(colorTriplets) || size(colorTriplets, 2) ~= 3
    error(['Input arg "colorTriplets" should be a matrix with three ',...
        'columns, with each row specifying an RGB color!']);
end

% check optional arguments, parse them
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && length(varargin{v})==1 && ismember(trimmingThr(t), 0:0.01:0.9)
            trimmingThr = varargin{v};
        elseif isnumeric(varargin{v}) && length(varargin{v})==3
            colorTriplet = varargin{v};
        elseif iscell(varargin{v}) && length(varargin{v}) == size(connmatrix, 1)
            labels = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'draw', 'nodraw'})
            drawFlag = varargin{v};
        elseif ischar(varargin{v}) && ~ismember(varargin{v}, {'draw', 'nodraw'})
            figTitle = varargin{v};            
        else
            error(['An input arg could not be parsed as any of "colorTriplet", ""trimmingThr", ',...
                '"labels", "drawFlag" or "figTitle"!']);
        end
    end
end
   
% defaults
if ~exist('colorTriplet', 'var')
    colorTriplet = [0, 0.447, 0.741];
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
if isrow(labels)
    labels = labels';
end
if iscolumn(colorTriplet)
    colorTriplet = colorTriplet';
end

% user message
disp([char(10), 'Function circleGraphPlotSimple is called with inputs: ',...
    char(10), 'Connectivity matrix with size ', num2str(size(connMatrix)),...
    char(10), 'Main graph color: ', num2str(colorTriplet),...
    char(10), 'Edge trimming theshold: ', num2str(trimmingThr),...
    char(10), 'Node label array with size ', num2str(size(labels)),...
    char(10), 'Figure display flag: ', num2str(drawFlag)]);


%% Compare labels to the specific anatomic ROI set 
% if they match, the function can highlight the ROI structure automatically

% expected set and order
expectedLabels = {'lateralorbitofrontal L', 'medialorbitofrontal L', 'parsorbitalis L', 'parstriangularis L', 'parsopercularis L', 'rostralmiddlefrontal L', 'caudalmiddlefrontal L', 'superiorfrontal L', 'precentral L',...
    'rostralanteriorcingulate L', 'caudalanteriorcingulate L', 'posteriorcingulate L', 'isthmuscingulate L',...
    'transversetemporal L', 'superiortemporal L', 'middletemporal L', 'inferiortemporal L', 'entorhinal L', 'parahippocampal L', 'fusiform L', 'insula L',...
    'supramarginal L', 'inferiorparietal L', 'superiorparietal L', 'postcentral L', 'paracentral L', 'precuneus L',...
    'cuneus L', 'lingual L', 'pericalcarine L', 'lateraloccipital L',...
    'lateraloccipital R', 'pericalcarine R', 'lingual R', 'cuneus R',...
    'precuneus R', 'paracentral R', 'postcentral R', 'superiorparietal R', 'inferiorparietal R', 'supramarginal R',...
    'insula R', 'fusiform R', 'parahippocampal R', 'entorhinal R', 'inferiortemporal R', 'middletemporal R', 'superiortemporal R', 'transversetemporal R',...
    'isthmuscingulate R', 'posteriorcingulate R', 'caudalanteriorcingulate R', 'rostralanteriorcingulate R',...
    'precentral R', 'superiorfrontal R', 'caudalmiddlefrontal R', 'rostralmiddlefrontal R', 'parsopercularis R', 'parstriangularis R', 'parsorbitalis R', 'medialorbitofrontal R', 'lateralorbitofrontal R',...
    };
shiftLabels = 15;
expectedLabels = [expectedLabels(end-shiftLabels:end), expectedLabels(1:end-shiftLabels-1)]';

% set flag if supplied arg "labels" match expected set of labels
if isequal(expectedLabels, labels)
    lobuleFlag = 1;
    disp([char(10), 'Supplied labels let us highlight the lobules with additional lines and annotations, will do so']);
else
    lobuleFlag = 0;
    disp([char(10), 'Cannot apply lobule-highlighting for given ROI/node label set']);
end


%% Hard-coded params

% base multiplier for the width of all edges (they are based on
% connectivity strength which is < 1)
edgeWidthMultip = 5;
% edge line styles for within- and between-module edges
%edgeTypes = {'-', 'none'};
edgeType = '-';
% general transparency setting for edges
edgeAlpha = 0.3;
% general node size setting
nodeSize = 10;
% graph plot layout
graphLayout = 'circle';
% figure (gcf) background color
gcfColor = [1 1 1];
% axes (gca) color in subplots
gcaLinesColor = [1 1 1];
% figure position and size in normalized units 
gcfMainPos = [0.25, 0, 0.5, 1];
% axes position relative to figure for main figure
gcaPosInFig = [0.05, 0.05, 0.9, 0.9];
% figure title texts
mainFigTitle = ['Full graph. ', figTitle];

% properties for text box displaying trimming info
trimmingText = ['Edges with weight > ', num2str(trimmingThr), ' are depicted'];
trimmingBoxPos = [0.01, 0.01, 0.4, 0.03];

% if ROIs are grouped into lobules, we draw lines and annotations, set
% their params
% ONLY WORKS PROPERLY ON 23" SCREEN WITH 1920*1080
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


%% Prepare connectivity matrix, init graph object, trim edges

% create symmetric adjacency matrix with zeros at diagonal
connMatrix = triu(connMatrix, 1) + triu(connMatrix, 1)';
connMatrix(isnan(connMatrix)) = 0;  % in many cases connectivity matrices contain NaN values

% create built-in graph object
if isempty(labels)
    G = graph(connMatrix, 'upper');
else
    G = graph(connMatrix, labels, 'upper');
end

% edge trimming
edgesToTrim = find(G.Edges.Weight < trimmingThr);  % graph.rmedge does not work with logical indexing, requires numeric edge indices
G = G.rmedge(edgesToTrim);

% set edge width values
edgeWidth = G.Edges.Weight*edgeWidthMultip;


%% Plot figure

% main graph plot figure
mainFig = figure;

% set figure size and background color
set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);
set(gcf, 'Color', gcfColor);

% graph plot
G.plot('Layout', graphLayout,... 
    'LineWidth', edgeWidth,... 
    'EdgeColor', colorTriplet,... 
    'EdgeAlpha', edgeAlpha,... 
    'NodeColor', colorTriplet,...
    'MarkerSize', nodeSize,...
    'LineStyle', edgeType);

% title
title(mainFigTitle);

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

