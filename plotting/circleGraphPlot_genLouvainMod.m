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
% roiLabels        - Cell array holding the labels for nodes / ROIs.
%           Defaults to loading var 'roisShort' from
%           'utils/roiNamesInOrder.mat'. Be careful as this 
% colorTriplets    - Numeric matrix with three columns, each row defines an
%           RGB triplet for module (and within-module edge) coloring.
%           Defaults to loading the var 'colorTriplets20' from
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
        colorTriplets = tmp.colorTriplets20;
    end
end
if ~exist('trimmingThr', 'var')
    trimmingThr = 0;
end
% further checks
if length(roiLabels) ~= size(realConn, 1)
    error('Length of ROI / node labels cell array does not match the number of nodes in "realConn"!');
end
if ~iscolumn(roiLabels)
    roiLabels = roiLabels';
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
% provide feedback about selected gamma & omega values
disp([char(10), 'Requested [gamma, omega] values were: ', num2str(gammaOmegaTargets),...
    char(10), 'Closest match in the data was: ',... 
    num2str([modRes.gammaValues(gIdx), modRes.omegaValues(oIdx)])]);

mod = modRes.res.consSim(gIdx, oIdx, :);
mod = double(squeeze(mod));  % res.consSim is uint16 by default if coming from multiCommDetectWrapper
% reshape modularity memberships into one module vector per layer/epoch
mod = reshape(mod, [nodeNo, layerNo]);


%% Check if labels match either of the expected label sets

% see if the supplied labels match any of the standard sets
[equalFlag, matchingSetsFlag, roiLabelsPlotting] = roiLabelMatching(roiLabels);


%% Check if we have enough unique colors for all modules

% no. of colors
colorNo = size(colorTriplets, 1);

% if there are more modules overall than colors, set flag for only keeping 
% color-module pairings in successive layers/epochs for overlapping modules
if sum(unique(mod(:))) > colorNo
    moduleRecolorFlag = true;
else
    moduleRecolorFlag = false;
end


%% Loop through layers/epochs, generate plots

% for layerIdx = 1:layerNo
for layerIdx = 1:10
    
    % title
    figTitle = ['Alpha, layer ', num2str(layerIdx), ', stim 1'];
    
    % connectivity
    connMatrix = squeeze(realConn(:,:,layerIdx));
    % modularity indices
    modIndicesVector = mod(:,layerIdx);    
    
    % reassign module numbers / colors if we would not have enough colors
    % to cover all modules otherwise
    if layerIdx > 1 && moduleRecolorFlag
        
        % there is already a mod2color, as it is defined for the first
        % layer
        mod2colorOld = mod2color;
        mod2color = [modIndicesVector, zeros(size(modIndicesVector, 1), 1)];
        % check if there is any module in the current layer/epoch not there
        % in the previous one
        newModules = setdiff(modIndicesVector, mod2colorOld(:, 1));
        % only go on if this set is not empty
        if ~isempty(newModules)
            % define the set of colors not used in the previous layer/epoch
            availColors = setdiff(1:colorNo, mod2colorOld(:, 2));
            % assign all new modules a color
            for newModIdx = 1:length(newModules)
                mod2color(modIndicesVector==newModules(newModIdx), 2) = repmat(availColors(newModIdx), [sum(modIndicesVector==newModules(newModIdx)), 1]);
            end  % for newModIdx
        end  % if ~isempty
        
    else
        mod2color = [modIndicesVector, modIndicesVector]; 
        
    end
    
    % if the labels were not an exact match but overlapped with standard
    % sets, rearrange labels + data
    if ~equalFlag && matchingSetsFlag
        % rearrange connectivity matrix based on new ROI label order
        [connMatrix, old2new] = matrixReorder(connMatrix, roiLabels, roiLabelsPlotting);
        % apply the same re-ordering to ROI/node module indices
        modIndicesVector = modIndicesVector(old2new);
        
        % plotting
        [mainFig, subFig] = circleGraphPlot(connMatrix, modIndicesVector,... 
                                colorTriplets, mod2color, trimmingThr,... 
                                roiLabelsPlotting, figTitle);
        
    % in any other case, simply call the plotting function with the original labels    
    else
        [mainFig, subFig] = circleGraphPlot(connMatrix, modIndicesVector,... 
                                colorTriplets, mod2color, trimmingThr,... 
                                roiLabels, figTitle); 
                            
    end  % if
    
    % save out plots
    mainFile = ['main_alpha_layer', num2str(layerIdx), '_stim1.png'];
    subFile = ['sub_alpha_layer', num2str(layerIdx), '_stim1.png'];
    saveas(mainFig, mainFile);
    saveas(subFig, subFile);
    
    % close figs
    close(mainFig);
    close(subFig);      
    
end  % for


return

                                  
                                  
                                  
                                  
                                  
                                  