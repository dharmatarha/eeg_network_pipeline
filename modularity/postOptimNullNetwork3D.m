function postOptimConn = postOptimNullNetwork3D(multiLayerConn, varargin)
%% Multilayer post-optimization null network calculation
%
% USAGE: postOptimConn = postOptimNullNetwork3D(multiLayerConn, intraLayerEdgeRand = 1, temporalRand = 0)
%
% Calculates the post-optimization null network of multi-layer connectivity
% data (e.g. multi-epoch brain connectvitiy data). Optional input args
% (randomization flags) control the type of null network we calculate.
% Internally, the function calls functions edgeRandomization and
% temporalRand. 
% 
% Mandatory inputs:
% multiLayerConn        - Numeric tensor, 3D. Multi-layer connectvitiy
%                         network with dimensions 
%                         (node no., node no., layer no.). Contains all
%                         intra-layer connection values. 
%
% Optional inputs:
% intraLayerEdgeRand    - Numeric value, one of 0:2. Flag indicating if 
%                         intra-layer edge weight randomization (random 
%                         reassignment of edge weights) should be performed.
%                         0 means no edge randomization, 1 means edge
%                         randomization across existing edges (constrained 
%                         version) and 2 means randomization across all
%                         possible edges (i.e., 0/NaN connectivity 
%                         values are also shuffled around). Defaults to 1.
% temporalRand          - Numeric value, one of 0:1. Flag indicating if 
%                         temporal randomization (random reshuffling of 
%                         network layers) should be performed. 0 means no
%                         reshuffling, 1 means reshuffling. Defaults to 0.
%
% Output:
% postOptimConn         - Multi-layer post-optimization null network with 
%                         the same dimensions as the original network (node
%                         no., node no., layer no.).
% 


%% Input checks

% Check for mandatory arguments
if ~ismembertol(nargin, 1:3)
    error(['Function postOptimNullNetwork3D requires input arg "multiLayerConn" ',...
        'while input args "intraLayerEdgeRand" and "temporalRand" are optional!']);
end
if (numel(size(multiLayerConn)) ~= 3) || size(multiLayerConn, 1) ~= size(multiLayerConn, 2)
    error('Input arg "multiLayerConn" should be 3D tensor with the first two dimension having equal size!');
end

% Check optional arguments
if ~isempty(varargin)
    % If too many, raise error
    if length(varargin) > 2
        error('Too many variable inputs. Only "intraLayerEdgeRand" and "temporalRand" are allowed!');
    % If only one, it is intraLayerEdgeRand
    elseif length(varargin) == 1
            intraLayerEdgeRand = varargin{1};
            temporalRand = 0;
    % If two, the first one is intraLayerEdgeRand, the second one is temporalRand
    elseif length(varargin) == 2
        intraLayerEdgeRand = varargin{1};
        temporalRand = varargin{2};
    end
    % check assinged values
    if ~ismembertol(intraLayerEdgeRand, 0:2)
        error('Optional input arg "intraLayerEdgeRand" should be one of 0:2!');
    end
    if ~ismembertol(temporalRand, 0:1)
        error('Optional input arg "temporalRand" should be one of 0:1!');
    end
else
    intraLayerEdgeRand = 1;
    temporalRand = 0;
end


%% Determine data dimensions  

[numberOfChannels, ~, numberOfLayers] = size(multiLayerConn);
% preallocate the tensor for the post-optimization null network
postOptimConn = zeros(numberOfChannels, numberOfChannels, numberOfLayers);


%% Calculate post-optimization null network

% if intra-layer edge randomization was requested
if intraLayerEdgeRand ~= 0
    % loop across layers
    for layerIndex = 1:numberOfLayers
        % if constrained edge randomization (only across existing edges)
        % was requested
        if intraLayerEdgeRand == 1
            postOptimConn(:, :, layerIndex) = edgeRandomization(squeeze(multiLayerConn(:, :, layerIndex)), true); 
        % else if randomization across all possible edges (including
        % non-existing edges) was requested
        elseif intraLayerEdgeRand == 2
            postOptimConn(:, :, layerIndex) = edgeRandomization(squeeze(multiLayerConn(:, :, layerIndex)), false);
        end
    end
else
    postOptimConn = multiLayerConn;
end

% temporal randomization (layer shuffling) as well, if requested
if temporalRand == 1
    postOptimConn = temporalRandomization(postOptimConn);
end 
    

return