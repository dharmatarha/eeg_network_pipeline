function [connMatrixRo, labelsRo, h, hs] = heatmapReordered(connMatrix, labels)
%% Simple heatmap with reordered label structure and lobule highlighting
%
% USAGE: [connMatrixRo, labelsRo, h, hs] = heatmapReordered(connMatrix, labels, cmap)
%
% A connectivity matrix is routinely plotted as a heatmap. The current
% function does that, but reorders first the matrix to a predefined
% lobule-grouped ROI structure. 
%
% IMPORTANT! Heatmap objects have changed over different Matlab releases,
% this function is intended for release 2017a. 
%
% Inputs:
% connMatrix    - Numeric square matrix, contains connectivity values.
%               Only upper triangle is taken into account. Might contain
%               NaN or zero values for missing connections, or a mix of
%               both (NaN and zero values are treated the same way).
% labels        - Cell vector, contains ROI/channel names as char arrays. 
%               Its length equals the number of rows/columns of "connMatrix".
%
% Outputs:
% connMatrixRo  - Numeric matrix, reordered version of "connMatrix". 
% labelsRo      - Cell vector, reordered version of "labels"
% h             - Figure handle generated by calling "heatmap".
% hs            - Output of struct(h).
%
% 
% TODO:
% - Allow for arbitrary reordering via new optional input arg "labelsRo"
% - Allow for arbitrary colormap via new optional input arg
% - Get lobule and laterality highlighting


%% Input checks

if nargin~=2
    error('Function heatmapReordered requires input args "connMatrix" and "labels"!');
end
if ~isnumeric(connMatrix) || ~ismatrix(connMatrix) || ~isequal(size(connMatrix, 1), size(connMatrix, 2))
    error('Input arg "connMatrix" should be a numeric square matrix!');
end
if ~iscell(labels) || ~isvector(labels) || ~isequal(length(labels), size(connMatrix, 1))
    error('Input arg "labels" should be a cell vector with its length equal to the rows/columns of "connMatrix"!');
end
if ~all(cellfun('isclass', labels, 'char'))
    error('Input arg "labels" should be a cell array of char arrays!');
end

% warning if release is not 2017a
if ~strcmp(version('-release'), '2017a')
    warning([char(10), 'The way this function treats heatmap objects is ',...
        'compatible with Matlab 2017a ', char(10),...
        'and might be incompatible with other releases. ',...
        char(10), 'Current version is not 2017a!']);
end

% user message
disp([char(10), 'Called heatmapReordered with input args: ',...
    char(10), 'Connectivity matrix with size: ', num2str(size(connMatrix)), ...
    char(10), 'ROI/channel labels in cell array with size: ', num2str(size(labels))]);


%% Reordering

% turn labels into row vector if column
if iscolumn(labels)
    labels = labels';
end

% Make sure connMatrix is symmetric before reordering it. Diagonal is
% treated as missing data (self-connectivity is meaningless)
connMatrix = triu(connMatrix, 1) + triu(connMatrix, 1)';
% Zeros mark missing connections, set them to NaN
connMatrix(connMatrix==0) = NaN;

% target labels (ROI structure with ROIs grouped according to lobules and
% laterality)
labelsRo = {'latOrbFront L', 'medOrbFront L', 'parsOrb L', 'parsTriang L', 'parsOpercul L', 'rostrMidFront L', 'caudMidFront L', 'supFront L', 'precentral L',...
    'rostrAntCing L', 'caudAntCing L', 'postCing L', 'isthmusCing L',...
    'transvTemp L', 'supTemp L', 'midTemp L', 'infTemp L', 'entorhinal L', 'paraHippoc L', 'fusiform L', 'insula L',...
    'supraMarg L', 'infPar L', 'supPar L', 'postcentral L', 'paracentral L', 'precuneus L',...
    'cuneus L', 'lingual L', 'periCalc L', 'latOcc L',...
    'latOcc R', 'periCalc R', 'lingual R', 'cuneus R',...
    'precuneus R', 'paracentral R', 'postcentral R', 'supPar R', 'infPar R', 'supraMarg R',...
    'insula R', 'fusiform R', 'paraHippoc R', 'entorhinal R', 'infTemp R', 'midTemp R', 'supTemp R', 'transvTemp R',...
    'isthmusCing R', 'postCing R', 'caudAntCing R', 'rostrAntCing R',...
    'precentral R', 'supFront R', 'caudMidFront R', 'rostrMidFront R', 'parsOpercul R', 'parsTriang R', 'parsOrb R', 'medOrbFront R', 'latOrbFront R',...
    };

% check is supplied labels correspond to target labels, reorder if true
if all(ismember(labelsRo, labels))
    [connMatrixRo, ~] = matrixReorder(connMatrix, labels, labelsRo);
else
    error('Could not reorder "connMatrix" as "labels" is an unexpected set!');
end


%% Heatmap

% plotting constants
gcfMainPos = [0, 0, 0.58, 96];
gcfBackground = [1 1 1];
missingEdgeLabel = 'No edge';
missingDataColor = [1 1 1];
fontSize = 14;

% heatmap 
h = heatmap(connMatrixRo, 'ColorMap', flipud(autumn),... 
    'MissingDataColor', missingDataColor,... 
    'MissingDataLabel', missingEdgeLabel,... 
    'FontSize', fontSize);

% background to white
set(gcf,'Color', gcfBackground);
% set size
set(gcf, 'Units', 'Normalized', 'OuterPosition', gcfMainPos);

% remove tick labels - release specific method accessing private properties:
old_warning_state = warning('off', 'MATLAB:structOnObject');
hs = struct(h);
warning(old_warning_state);
hs.XAxis.TickValues = [];
hs.YAxis.TickValues = [];


return

