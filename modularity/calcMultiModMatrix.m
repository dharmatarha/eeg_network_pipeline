function B = calcMultiModMatrix(realConn, nullConn, varargin)
%% Get multilayer connectivity matrix from a set of 2D connectivity matrices
%
% USAGE: B = calcMultiModMatrix(realConn, nullConn, gamma = 1, omega = 1)
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
% Mandatory inputs:
% realConn  -       3D numeric array, contains a set of connectivity matrices.
%                   Its dimensions are 
%                   (nodeNo, nodeNo, epochNo).
% nullConn   -      2D or 3D numeric array, contains null (e.g. surrogate) 
%                   connectivity matrix / matrices. If 2D, the same null
%                   matrix is used for all real connectivity matrices.
%                   Its dimensions are 
%                   (nodeNo, nodeNo (, epochNo)).
%
% Optional inputs:
% gamma               - Numeric value, spatial resolution parameter.
%                     Defaults to 1.
% omega               - Numeric value, inter-layer edge weight. Uniform
%                     within-node inter-layer connections are assumed (all 
%                     inter-layer edges have the same weight). Defaults to 1.
%
% Output:
% B                   - 2D numeric matrix, sized
%                     (nodeNo*epochNo,
%                     nodeNo*epochNo). Multilayer
%                     connectivity matrix, where connectivity value between node "n1" from
%                     epoch "e1" and node "n2" from epoch "e2" is stored at 
%                     (numberOfChannels*(e1-1)+n1,
%                     numberOfChannels*(e2-1)+m2). Off-diagonal submatrices
%                     represent inter-layer connectivity values (see
%                     omega).
% 
%

%% Input checks

% Check for mandatory arguments
if nargin < 2
    error('Input args "realConn" and "nullConn" are required!');
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
    if length(varargin) > 2
        error('Too many variable inputs. Only "gamma" and "omega" are allowed!');
    % If only one, it is gamma
    elseif length(varargin) == 1
            gamma = varargin{1};
            omega = 1;
    % If two, the first one is gamma, the second one is omega
    elseif length(varargin) == 2
        gamma = varargin{1};
        omega = varargin{2};
    end
else
    % default values
    gamma = 1;
    omega = 1;
end


%% Determine data dimensions  

[nodeNo, ~, epochNo] = size(realConn);


%% Normalize real and null connectivity data, apply gamma

% Preallocate 3 array for intra-layer connectivity matrices after
% null-model extraction
intraLayerConnArr = zeros(epochNo, nodeNo, nodeNo);

% normalize connectivity values then extract null model for each epoch
for epochIdx = 1 : epochNo

    % Real and null connectivity matrices for a given epoch
    realConnMatrix = squeeze(realConn(:, :, epochIdx));
    if length(size(nullConn)) == 2
        nullConnMatrix = nullConn;
    else
        nullConnMatrix = squeeze(nullConn(:, :, epochIdx));
    end

    % Keep only connections in both measured and surrogate matrices, where measured connectivity exists
    nullConnMatrix(isnan(realConnMatrix)) = 0;
    realConnMatrix(isnan(realConnMatrix)) = 0;

    % Normalize measured and surrogate connectivity matrices
    normedConnMatrix = normalizeMatrix(realConnMatrix);
    normedNullMatrix = normalizeMatrix(nullConnMatrix);

    % Intra-layer modularity matrix
    intraLayerConnMatrix = normedConnMatrix - gamma * normedNullMatrix;

    % If the intra-layer modularity matrix is not symmetric, symmetrization is forced
    if nnz(intraLayerConnMatrix-intraLayerConnMatrix')
        intraLayerConnMatrix = (intraLayerConnMatrix + intraLayerConnMatrix')/2; %disp('WARNING: Forced symmetric intra-layer modularity matrix ')
    end

    intraLayerConnArr(epochIdx, :, :) = intraLayerConnMatrix;
end


%% Get multilayer connectivity as one 2D matrix, add omega

% Allocate sparse matrix for the multilayer (output) matrix
B = spalloc(nodeNo*epochNo, nodeNo*epochNo, nodeNo*nodeNo*epochNo+2*nodeNo*epochNo);

% Iterate over all layers
for s = 1:epochNo
    indx = [1:nodeNo] + (s-1)*nodeNo;
    % Diagonal blocks are the intra-layer modularity matrices
    B(indx, indx) = squeeze(intraLayerConnArr(s, :, :));
end
% Off-diagonal blocks contain inter-layer connections
B = B + omega*spdiags(ones(nodeNo*epochNo,2), [-nodeNo, nodeNo], nodeNo*epochNo, nodeNo*epochNo);


end