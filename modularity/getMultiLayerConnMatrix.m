function multiLayerConn = getMultiLayerConnMatrix(realConn, nullConn, varargin)
%% Get multilayer connectivity matrix from a set of 2D connectivity matrices
%
% USAGE: multiLayerConn = getMultiLayerConnMatrix(realConn, nullConn, gamma = 1, omega = 1, normalization = 'oneStep')
%
% Calculates the multilayer modularity matrix of 3D connectivity data 
% (realConn - usually a set of connectivity matrices, each corresponding 
% to one epoch).
% The output can be fed to the genlouvain.m function of the GenLouvain
% toolbox.
% 
% The function subtracts a null model from the connectivity data before 
% generating the multilayer matrix. As null model it expects an array with 
% the same size as the real connectivity data. Alternatively it accepts 
% one connectivity matrix that can be used as null model for all real 
% connectivity matrices. 
% The spatial resolution parameter (gamma) and the inter-layer edge weight 
% (omega) might also be specified. They are set to 1 by default but that is
% virtually guaranteed not to fit any real use case.
%
% In order to keep the gamma / omega values in interpretable ranges there
% is a normalization/rescaling step applied before subtracting the null
% model from the real connectivity data. The exact behavior might be
% controlled with the "normalization" arg, but by default, the function
% rescales both "realConn" and "nulConn" with their overall mean values
% (i.e., realConn/normalizeTensor(realConn, 'mean')).
%
% There is also a step where the null model gets masked according to existing
% (non-NaN) edges in the real connectivity data. This is included because 
% in our use cases the real data is already thresholded (based on
% permutation tests against distributions of null model values), resulting
% in NaN values, whereas the null model is always full (in the upper
% triangle, that is).
% 
% Mandatory inputs:
% realConn      - 3D numeric array, contains a set of connectivity matrices.
%                 Its dimensions are (nodeNo, nodeNo, epochNo).
% nullConn      - 2D or 3D numeric array, contains null (e.g. surrogate) 
%                 connectivity matrix / matrices. If 2D, the same null
%                 matrix is used for all real connectivity matrices.
%                 Its dimensions are (nodeNo, nodeNo (, epochNo)).
%
% Optional inputs:
% gamma         - Numeric value, spatial resolution parameter.
%                 Defaults to 1.
% omega         - Numeric value, inter-layer edge weight. Uniform
%                 within-node inter-layer connections are assumed (all 
%                 inter-layer edges have the same weight). Defaults to 1.
% normalization - Char array, one of {'perLayer', 'oneStep'}. Determines
%                 the type of normalization / rescaling procedure applied 
%                 to the connectivity and null arrays. 'perLayer' means 
%                 that normalization is applied separately for each layer 
%                 (using normalizeMatrix), while 'oneStep' means 
%                 normalization of 3D arrays in one step (with
%                 normalizeTensor).
%
% Output:
% B             - 2D numeric matrix, sized
%                 (nodeNo*epochNo, nodeNo*epochNo). 
%                 Multilayer connectivity matrix, where connectivity 
%                 value between node "n1" from epoch "e1" and node 
%                 "n2" from epoch "e2" is stored at (numberOfChannels*(e1-1)+n1,
%                 numberOfChannels*(e2-1)+m2). Off-diagonal submatrices
%                 represent inter-layer connectivity values (see
%                 omega).
% 
%

%% Input checks

% Check for mandatory arguments
if ~ismembertol(nargin, 2:5)
    error(['Function getmultiLayerConnMatrix requires ipnut args "realConn" ',...
        'and "nullConn", while input args "gamma", "omega" and "normalization" are optional!']);
end
if (numel(size(realConn)) ~= 3) || (size(realConn, 1) ~= size(realConn, 2))
    error('Input arg "realConn" has unexpected dimensions (should be (nodeNo, nodeNo, epochNo))!');
end
if ~isequal(size(nullConn), size(realConn)) && ~isequal(size(nullConn), [size(realConn, 1), size(realConn, 2)])
    error('Input arg "nullConn" has unexpected dimensions (should be (nodeNo, nodeNo (, epochNo)))!');
end
% Check optional arguments
if ~isempty(varargin)
    % If too many, raise error
    if length(varargin) == 1
            gamma = varargin{1};
            omega = 1;
            normalization = 'oneStep';
    % If two, the first one is gamma, the second one is omega
    elseif length(varargin) == 2
        gamma = varargin{1};
        omega = varargin{2};
        normalization = 'oneStep';
    % If three, the first one is gamma, the second one is omega    
    elseif length(varargin) == 3
        gamma = varargin{1};
        omega = varargin{2};        
        normalization = varargin{3};
    end
else
    % default values
    gamma = 1;
    omega = 1;
    normalization = 'oneStep';
end
% check for empty values
if isempty(gamma); gamma = 1; end
if isempty(omega); omega = 1; end
if isempty(normalization); normalization = 'oneStep'; end
% check optional input arg values
if ~isnumeric(gamma) || numel(gamma)~=1
    error('Input arg "gamma" should be a numeric value!');
end
if ~isnumeric(omega) || numel(omega)~=1
    error('Input arg "omega" should be a numeric value!');
end
if ~ismember(normalization, {'oneStep', 'perLayer'})
    error('Input arg "normalization" should be one of {''oneStep'', ''perLayer''}!');
end
% normalization flag
if strcmp(normalization, 'oneStep'); oneStepFlag = true; else; oneStepFlag = false; end


%% Determine data dimensions  

[nodeNo, ~, layerNo] = size(realConn);


%% Normalize real and null connectivity data

% if normalization is in one step
if oneStepFlag
    
    % real connectivity data
    realConnNorm = normalizeTensor(realConn, 'mean', false);  % suppress NaN warnings
    % null model data
    if length(size(nullConn)) == 2
        nullConn = repmat(nullConn, [1, 1, layerNo]);
    end
    % selecting edges existing (non-NaN) in real data 
    nullConn(isnan(realConnNorm)) = NaN;
    % normalization
    nullConnNorm = normalizeTensor(nullConn, 'mean', false);  % suppress NaN warnings
    
% if normalization is per layer    
else
    % preallocate for normalized arrays
    realConnNorm = zeros(nodeNo, nodeNo, layerNo);
    nullConnNorm = realConnNorm;
    % loop through layers / epochs
    for layerIdx = 1 : layerNo
        % Real and null connectivity matrices for a given epoch
        realConnMatrix = squeeze(realConn(:, :, layerIdx));
        if length(size(nullConn)) == 2
            nullConnMatrix = nullConn;
        else
            nullConnMatrix = squeeze(nullConn(:, :, layerIdx));
        end    
        
        % Keep only connections in surrogate matrices where measured
        % connectivity exists in real data
        nullConnMatrix(isnan(realConnMatrix)) = NaN;

        % Normalize measured and surrogate connectivity matrices
        realConnNorm(:, :, layerIdx) = normalizeMatrix(realConnMatrix);
        nullConnNorm(:, :, layerIdx) = normalizeMatrix(nullConnMatrix);       
    end  % for layerIdx
    
end  % if oneStepFlag


%% Apply gamma

intraLayerConnArr = realConnNorm-gamma*nullConnNorm;
        

%% Symmetrization of all layers

for layerIdx = 1 : layerNo
    
    % if layer is not symmetric, upper triangle is copied to lower triangle
    if nnz(intraLayerConnArr(:, :, layerIdx)-intraLayerConnArr(:, :, layerIdx)')
        intraLayerConnArr(:, :, layerIdx) = triu(intraLayerConnArr(:, :, layerIdx)) + triu(intraLayerConnArr(:, :, layerIdx), 1)';
    end
    
end  % for layerIdx


%% Turn NaN values to zero

intraLayerConnArr(isnan(intraLayerConnArr)) = 0;


%% Get multilayer connectivity as one 2D matrix, add omega

% Allocate sparse matrix for the multilayer (output) matrix
multiLayerConn = spalloc(nodeNo*layerNo, nodeNo*layerNo, nodeNo*nodeNo*layerNo+2*nodeNo*layerNo);

% Iterate over all layers
for layerIdx = 1:layerNo
    indx = [1:nodeNo] + (layerIdx-1)*nodeNo;
    % Diagonal blocks are the intra-layer modularity matrices
    multiLayerConn(indx, indx) = squeeze(intraLayerConnArr(:, :, layerIdx));
end

% Off-diagonal blocks contain inter-layer connections
multiLayerConn = multiLayerConn + omega*spdiags(ones(nodeNo*layerNo,2), [-nodeNo, nodeNo], nodeNo*layerNo, nodeNo*layerNo);



return



% % Earlier version:
%
% % Preallocate 3D array for intra-layer connectivity matrices after
% % null-model extraction
% intraLayerConnArr = zeros(layerNo, nodeNo, nodeNo);
%     
% 
% % normalize connectivity values then extract null model for each epoch
% for layerIdx = 1 : layerNo
% 
%     % Real and null connectivity matrices for a given epoch
%     realConnMatrix = squeeze(realConn(:, :, layerIdx));
%     if length(size(nullConn)) == 2
%         nullConnMatrix = nullConn;
%     else
%         nullConnMatrix = squeeze(nullConn(:, :, layerIdx));
%     end
% 
%     % Keep only connections in both measured and surrogate matrices, where measured connectivity exists
%     nullConnMatrix(isnan(realConnMatrix)) = 0;
%     realConnMatrix(isnan(realConnMatrix)) = 0;
% 
%     % Normalize measured and surrogate connectivity matrices
%     normedConnMatrix = normalizeMatrix(realConnMatrix);
%     normedNullMatrix = normalizeMatrix(nullConnMatrix);
% 
%     % Intra-layer modularity matrix
%     intraLayerConnMatrix = normedConnMatrix - gamma * normedNullMatrix;
% 
%     % If the intra-layer modularity matrix is not symmetric, symmetrization is forced
%     if nnz(intraLayerConnMatrix-intraLayerConnMatrix')
%         intraLayerConnMatrix = (intraLayerConnMatrix + intraLayerConnMatrix')/2; %disp('WARNING: Forced symmetric intra-layer modularity matrix ')
%     end
% 
%     intraLayerConnArr(layerIdx, :, :) = intraLayerConnMatrix;
% end