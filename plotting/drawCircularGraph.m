function drawCircularGraph(adjacencyMatrix, varargin)

%% Draw undirected circular graph specified by its adjacency matrix
%
% USAGE: drawCircularGraph(adjacencyMatrix, colorMatrix = zeros(N, N), textLabels = {'', '', ...}, lobuleFlag = 0)
%
% ATTENTION: Reorders rows and columns of the adjacency matrix and color
%            matrix as follows:
%                   permutationVector = [(32:62) (1:31)];
%                   adjacencyMatrix = adjacencyMatrix(permutationVector, :);
%                   adjacencyMatrix = adjacencyMatrix(:, permutationVector);
%                   colorMatrix = colorMatrix(permutationVector, :);
%                   colorMatrix = colorMatrix(:, permutationVector);
%
% Draws an undirected circular graph specified by its adjacency matrix,
% optionally with specified colors, textlabels and lobule labels.
%
% Mandatory input:
% adjacencyMatrix  - N by N numeric square adjacency matrix of the graph.
%               Only the upper triangular is considered.
% 
% Optional inputs:
% colorMatrix      - N by N numeric square matrix specifying color for each
%               edge by one numeric value. Must have the same size as the
%               adjacency matrix. Edge colors will be proportional to the
%               values given in the colorMatrix. Only the upper triangular
%                is considered. Defaults to deep blue color for each edge.
% textLabels      - Cell array listing text labels for each node. Must have
%               the same number of elements as the number of rows in the
%               adjacency matrix. Defaults to empty string for each node.
% lobuleFlag      - Flag indicating whether lobules are shown by separating
%               lines and labels. Defaults to 0 (lobules not shown).
% 
% Based on the great implementation of Paul Kassebaum:
% https://github.com/paul-kassebaum-mathworks/circularGraph


%% Input checks

% check for mandatory argument
if nargin < 1
    error('Input arg adjacencyMatrix is required!');
end

% check mandatory inputs
if ~isnumeric(adjacencyMatrix) || ~ismatrix(adjacencyMatrix) || size(adjacencyMatrix, 1) ~=size (adjacencyMatrix, 2)
    error('Input arg "adjacencyMatrix" should be a numeric symmetric matrix!');
end

% check optional arguments, parse them
if ~isempty(varargin)
    for v = 1:length(varargin)
        if isnumeric(varargin{v}) && ismatrix(varargin{v}) && size(varargin{v}, 1)>1 && ~exist('colorMatrix', 'var')
            colorMatrix = varargin{v};
        elseif iscell(varargin{v}) && length(varargin{v}) == size(adjacencyMatrix, 1) && ~exist('textLabels', 'var')
            textLabels = varargin{v};
        elseif isnumeric(varargin{v}) && ismember(varargin{v}, [0 1]) && ~exist('lobuleFlag', 'var')
            lobuleFlag = varargin{v};          
        else
            error(['An input arg could not be parsed as any of "colorMatrix", ',...
                '"textLabels", or "lobuleFlag"!']);
        end
    end
end

% Check if the color matrix has the same size as the adjacency matrix
% (if specified).
if exist('colorMatrix', 'var')
    if size(colorMatrix, 1) ~= size(adjacencyMatrix, 1) || ...
            size(colorMatrix, 2) ~= size(adjacencyMatrix, 2)
        error('The color matrix has to have the same size as the adjacency matrix');
    end
end

% defaults
if ~exist('colorMatrix', 'var')
    colorMatrix = zeros(size(adjacencyMatrix));
end
if ~ exist('textLabels', 'var')
    for i = 1 : size(adjacencyMatrix, 1)
        textLabels{i} = blanks(1);
    end
end
if ~ exist('lobuleFlag', 'var')
    lobuleFlag = 0;
end

% user message
disp([char(10), 'Function drawCircularGraph is called with inputs: ',...
    char(10), 'Connectivity matrix with size ', num2str(size(adjacencyMatrix)),...
    char(10), 'Edge coloring matrix with size ', num2str(size(adjacencyMatrix)),...
    char(10), 'Node label array with size ', num2str(size(textLabels)),...
    char(10), 'Show lobules flag: ', num2str(lobuleFlag)]);


%% Basics

% Symmetrize adjacency and color matrix
adjacencyMatrix = (triu(adjacencyMatrix, 1) + triu(adjacencyMatrix, 1)');
colorMatrix = (triu(colorMatrix, 1) + triu(colorMatrix, 1)');

% Permute rows and columns according to our node ordering rule
permutationVector = [(32:62) (1:31)];
adjacencyMatrix = adjacencyMatrix(permutationVector, :);
adjacencyMatrix = adjacencyMatrix(:, permutationVector);
colorMatrix = colorMatrix(permutationVector, :);
colorMatrix = colorMatrix(:, permutationVector);

% textLabels = {'midTemp L', 'infTemp L', 'entorhinal L', 'paraHippoc L', 'fusiform L', 'insula L',...
%     'supraMarg L', 'infPar L', 'supPar L', 'postcentral L', 'paracentral L', 'precuneus L',...
%     'cuneus L', 'lingual L', 'periCalc L', 'latOcc L',...
%     'latOcc R', 'periCalc R', 'lingual R', 'cuneus R',...
%     'precuneus R', 'paracentral R', 'postcentral R', 'supPar R', 'infPar R', 'supraMarg R',...
%     'insula R', 'fusiform R', 'paraHippoc R', 'entorhinal R', 'infTemp R', 'midTemp R', 'supTemp R', 'transvTemp R',...
%     'isthmusCing R', 'postCing R', 'caudAntCing R', 'rostrAntCing R',...
%     'precentral R', 'supFront R', 'caudMidFront R', 'rostrMidFront R', 'parsOpercul R', 'parsTriang R', 'parsOrb R', 'medOrbFront R', 'latOrbFront R',...
%     'latOrbFront L', 'medOrbFront L', 'parsOrb L', 'parsTriang L', 'parsOpercul L', 'rostrMidFront L', 'caudMidFront L', 'supFront L', 'precentral L',...
%     'rostrAntCing L', 'caudAntCing L', 'postCing L', 'isthmusCing L',...
%     'transvTemp L', 'supTemp L'
%     };


%% Draw nodes

t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
nodeLines = [];
extents = [];
labelOffsetFactor = 1.1;
for i = 1:length(adjacencyMatrix)
    x = cos(t(i));
    y = sin(t(i));
    nodeLines(end+1) = line(...
        cos(t(i)),...
        sin(t(i)),...
        2,...
        'Color', 'k',...
        'Marker','o',...
        'LineStyle','none',...
        'PickableParts','all');
    tau = atan2(y,x);
    textLabel = text(0,0,textLabels{i});
    textLabel.Position = labelOffsetFactor*([x, y]);
    if abs(tau) > pi/2
        textLabel.Rotation = 180*(tau/pi + 1);
        textLabel.HorizontalAlignment = 'right';
    else
        textLabel.Rotation = tau*180/pi;
    end
    extents(i) = textLabel.Extent(3);
end


%% Draw connections

% Find non-zero values of the adjacencyMatrix and their indices
[row,col,v] = find(adjacencyMatrix);

% Select the same values from the colorMatrix too
relevantColorValues = diag(colorMatrix(row, col));

% Rescale color values
if max(relevantColorValues) - min(relevantColorValues) > 0
    scaledColorValues = (relevantColorValues - min(relevantColorValues) ) ./ ...
    (max(relevantColorValues) - min(relevantColorValues));
else
    scaledColorValues = relevantColorValues - min(relevantColorValues);
end

% Scaled color values are used to select color for each edge from 1000
% parula values. The color is proportional to the scaled color value.
colorMap = parula(1000);

% Calculate line widths based on values of s (stored in v).
minLineWidth  = 0.5;
lineWidthCoef = 5;
lineWidth = v./max(v);
if sum(lineWidth) == numel(lineWidth) % all lines are the same width.
    lineWidth = repmat(minLineWidth,numel(lineWidth),1);
else % lines of variable width.
    lineWidth = lineWidthCoef*lineWidth + minLineWidth;
end
      
% Draw connections on the Poincare hyperbolic disk.
%
% Equation of the circles on the disk:
% x^2 + y^2
% + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x
% - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
% where u and v are points on the boundary.
%
% Standard form of equation of a circle
% (x - x0)^2 + (y - y0)^2 = r^2
%
% Therefore we can identify
% x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
% y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
% r^2 = x0^2 + y0^2 - 1

connections = [];

for i = 1:length(v)
    if row(i) ~= col(i)
        if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0
            % points are diametric, so draw a straight line
            u = [cos(t(row(i)));sin(t(row(i)))];
            v = [cos(t(col(i)));sin(t(col(i)))];
            connections(end+1) = line(...
                [u(1);v(1)],...
                [u(2);v(2)],...
                'LineWidth', lineWidth(i),...
                'Color', colorMap(round(999*scaledColorValues(i))+1, :),...
                'PickableParts','none');
        else % points are not diametric, so draw an arc
            u  = [cos(t(row(i)));sin(t(row(i)))];
            v  = [cos(t(col(i)));sin(t(col(i)))];
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
            
            if u(1) >= 0 && v(1) >= 0
                % ensure the arc is within the unit disk
                theta = [linspace(max(thetaLim),pi,50),...
                    linspace(-pi,min(thetaLim),50)].';
            else
                theta = linspace(thetaLim(1),thetaLim(2)).';
            end
            
            connections(end+1) = line(...
                r*cos(theta)+x0,...
                r*sin(theta)+y0,...
                'LineWidth', lineWidth(i),...
                'Color', colorMap(round(999*scaledColorValues(i))+1, :),...
                'PickableParts','none');
            
        end % if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0
        
    end % if row(i) ~= col(i)
    
end % for i = 1:length(v)


%% Set axes properties

axis image;
ax = gca;
extent = max(extents(:));
fudgeFactor_x = 2; % Not sure why this is necessary. Eyeballed it.
fudgeFactor_y = 2.2; % Not sure why this is necessary. Eyeballed it.
ax.XLim = ax.XLim + fudgeFactor_x*extent*[-1 1];
ax.YLim = ax.YLim + fudgeFactor_y*extent*[-1 1];
ax.Visible = 'off';
ax.SortMethod = 'depth';

fig = gcf;
fig.Color = [1 1 1];


%% Draw lobules if needed

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
        0.10, 0.59, 0.1, 0.1;
        0.32, 0.78, 0.1, 0.1; 
        0.66, 0.78, 0.1, 0.1;
        0.85, 0.59, 0.1, 0.1;
        0.87, 0.35, 0.1, 0.1;
        0.77, 0.15, 0.1, 0.1;
        0.57, 0.06, 0.1, 0.1];
    % lobule labels
    textA = {'L Occipital', 'L Parietal', 'L Temporal', 'L Cingulate', 'L Frontal',... 
        'R Frontal', 'R Cingulate', 'R Temporal', 'R Parietal', 'R Occipital'};
    
end

% axes (gca) color in subplots
gcaLinesColor = [1 1 1];

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

% figure position and size in normalized units 
gcfMainPos = [0, 0, 0.5, 1];
set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);

end