function circleGraphPlot_genLouvainMod(realConn, modRes, gammaOmegaTargets, varargin)
%% Wrapper for plotting modularity results from genLouvain 
% 
% USAGE: circleGraphPlot_genLouvainMod(realConn, 
%                                       modRes, 
%                                       gammaOmegaTargets, 
%                                       roiLabels=load('utils/roinamesInOrder.mat', 'roisShort'),
%                                       colorTriplets=load('utils/colorTriplets.mat', 'colorTriplets'),
%                                       trimmingThr=0)
%
% Wrapper function for plotting the results of modularity detection on a
% set of epochs using genLouvain. Call circleGraphPlotting internally,
% after selecting and rearranging module memberships.
%
% Mandatory inputs:
% realConn          - 3D numeric array, contains the connectivity / adjacency
%           matrices for plotting. Its size is 
%           [node no. X node no. X epoch no.], that is, the first two
%           dimensions define the connectivity square matrix for a given
%           epoch.
% modRes            - Struct, output of multiCommDetectWrapper. While it has
%           multiple fields, this function relies on fields "gammaValues",
%           "omegaValues" and "res.consSim". res.consSim is (usually)
%           a 3D numeric array sized [gamma param no. X omega param
%           no. X modularityVector size] that contains the consensus module
%           membership of each node in each layer / epoch, in a
%           vectorized form. 
% gammaOmegaTargets - Numeric vector with length==2. Contains the gamma and
%           omega parameter values for the modularity detection algorithm. 
%           Usually we run multiCommDetectWrapper with a set of gamma and 
%           omega values, and modRes.consSim contains the module
%           memberships for each [gamma, omega] pairing. This input arg
%           controls which module membership array we work with, selecting
%           it on the basis of corresponding [gamma, omega] param values.
%
% Optional inputs:
% roiLabels        - Cell array holding the labels for nodes / ROIs. Be
%           Defaults to loading var 'roisShort' from
%           'utils/roiNamesInOrder.mat'. Be careful as this 
% colorTriplets    - Numeric matrix with three columns, each row defines an
%           RGB triplet for module (and within-module edge) coloring.
%           Defaults to loading the var 'colorTriplets' from
%           utils/colorTriplets.mat.
% trimmingThr      - One- or two-element numeric vector containing threshold(s) 
%           for trimming (deleting) weak connections before plotting.
%           If trimmingThr is only one value, the same threshold is
%           applied to all connections. If two values, the first one is
%           applied to within-module, the second to between-module
%           connections. Defaults to [0.2], value(s) must be in range
%           [0:0.001:0.9].
%


%% Input checks

% no. of inputs
if ~ismember(nargin, 3:6)
    error(['Function circleGraphPlot_genLouvainMod requires input args ',...
        '"realConn", "modRes" and "gammaOmegaTargets" while args "roiLabels", ',...
        '"colorTriplets" and "trimmingThr" are optional!']);
end
% mandatory inputs
if ~isnumeric(realConn) || length(size(realConn)) ~= 3 || size(realConn, 1) ~= size(realConn, 2)
    error(['Input arg "realConn" should be a 3D numeric array with equal ',...
        'sizes on first two dimensions (one square matrix per layer / epoch)!']);
end
if ~isstruct(modRes) || ~isfield(modRes, 'res') || ~isfield(modRes.res, 'consSim')
    error('Input arg "modRes" should be a struct containing a 3D numeric field "res.consSim"!');
end
if size(modRes.res.consSim, 3) ~= size(realConn, 1)*size(realConn, 3)
    error(['The field "res.consSim" of input arg "modRes" has unexpected size ',...
        '- its third dimension should be (node no X layer/epoch no), ',...
        'matching the node no. and layer/epoch no. inferred from input arg "realConn"!']);
end
if ~isnumeric(gammaOmegaTargets) || ~ismember(numel(gammaOmegaTargets), 1:2)
    error('Input arg "gammaOmegaTargets" should be a one- or two-element numeric vector!');
end
% optional inputs
if ~isempty(varargin)
    for v = 1:length(varargin)
        if iscell(varargin{v}) && isvector(varargin{v}) && ~exist('roiLabels', 'var')
            roiLabels = varargin{v};
        elseif isnumeric(varargin{v}) && ismatrix(varargin{v}) && size(varargin{v}, 2)==3 && ~exist('colorTriplets', 'var')
            colorTriplets = varargin{v};
        elseif isnumeric(varargin{v}) && ismember(numel(varargin{v}), 1:2) && ~exist('trimmingThr', 'var')
            trimmingThr = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to "roiLabels", "colorTriplets" or "trimmingThr"!');
        end  % if
    end  % for v
end  % if
% default values
if ~exist('roiLabels', 'var')
    tmp = which('roiNamesInOrder.mat', '-ALL');
    if size(tmp, 1) ~= 1
        error('Found either zero or multiple versions of "roiNamesInOrder.mat"! See the help on input arg "roiLabels"!');
    else
        tmp = load('roiNamesInOrder.mat');
        roiLabels = tmp.roisShort;
    end
end
if ~exist('colorTriplets', 'var')
    tmp = which('colorTriplets.mat', '-ALL');
    if size(tmp, 1) ~= 1
        error('Found either zero or multiple versions of "colorTriplets.mat"! See the help on input arg "colorTriplets"!');
    else
        tmp = load('colorTriplets.mat');
        colorTriplets = tmp.colorTriplets;
    end
end
if ~exist('trimmingThr', 'var')
    trimmingThr = 0;
end
% further checks
if length(roiLabels) ~= size(realConn, 1)
    error('Length of ROI / node labels cell array does not match the number of nodes in "realConn"!');
end

% user message
disp([char(10), 'Called circleGraphPlot_genLouvainMod with input args: ',...
    char(10), 'Connectivity array sized ', num2str(size(realConn)),...
    char(10), 'Modularity struct with fields ',...
    char(10), '      res.consSim sized ', num2str(size(modRes.res.consSim)),...
    char(10), '      gammaValues sized ', num2str(size(modRes.gammaValues)),...
    char(10), '      omegaValues sized ', num2str(size(modRes.omegaValues)),...
    char(10), 'ROI/node labels as cell array sized ', num2str(size(roiLabels)),...
    char(10), 'Coloring RGB triplets array sized ', num2str(size(colorTriplets)),...
    char(10), 'Trimming threshold(s): ', num2str(trimmingThr)]);


%% Select and rearrange module memberships

% get node and layer numbers
[nodeNo, ~, layerNo] = size(realConn);

% get modularity around target gamma, omega pairing
[~, gIdx] = min(abs(modRes.gammaValues-gammaOmegaTargets(1)));
[~, oIdx] = min(abs(modRes.omegaValues-gammaOmegaTargets(2)));
mod = modRes.res.consSim(gIdx, oIdx, :);
mod = double(squeeze(mod));  % res.consSim is uint16 by default if coming from multiCommDetectWrapper
% reshape modularity memberships into one module vector per layer/epoch
mod = reshape(mod, [nodeNo, layerNo]);

% rearrange label names and modularity indices
% expected set of labels
newLabels = {'lateralorbitofrontal L', 'medialorbitofrontal L', 'parsorbitalis L', 'parstriangularis L', 'parsopercularis L', 'rostralmiddlefrontal L', 'caudalmiddlefrontal L', 'superiorfrontal L', 'precentral L',...
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
newLabels = [newLabels(end-shiftLabels:end), newLabels(1:end-shiftLabels-1)];

% get intersect of supplied labels and expected ones
c = intersect(labels, newLabels);

% if the two sets match, set data transformation flag
if isequal(c, labels)
    transformFlag = 1;  % transformation flag
else
    error('Unexpected labels');
end


for layerIdx = 1:40
    
    % title
    figTitle = ['Alpha, layer ', num2str(layerIdx), ', stim 1'];
    
    % connectivity
    connMatrix = squeeze(connData(:,:,layerIdx));
    % modularity indices
    modIndicesVector = modAll(:,layerIdx);

    % transform connectivity and module data in case of special label set
    if transformFlag   
    
        % rearrange connectivity matrix based on new ROI/node label
        % order
        [connMatrix, old2new] = matrixReorder(connMatrix, labels, newLabels);
        % apply the same re-ordering to ROI/node module indices
        modIndicesVector = modIndicesVector(old2new);

        % plot
        [mainFig, subFig] = circleGraphPlot(connMatrix, modIndicesVector, colorTriplets, trimmingThr,... 
                                              newLabels, figTitle);
                                          
    end
    
    % save out plots
    mainFile = ['main_alpha_layer', num2str(layerIdx), '_stim1.png'];
    subFile = ['sub_alpha_layer', num2str(layerIdx), '_stim1.png'];
    saveas(mainFig, mainFile);
    saveas(subFig, subFile);
    
    % close figs
    close(mainFig);
    close(subFig);
                                          
end

                                  
                                  
                                  
                                  
                                  
                                  